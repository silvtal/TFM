#!/usr/bin/env Rscript

start.time <- Sys.time() # para devolver al final el tiempo de ejecución


# -------------------------
#        FUNCIONES         
# -------------------------

library("parallel")
library("optparse")
home <- strsplit(paste0("./",getopt::get_Rscript_filename()),split="/")[[1]]
home <- paste(home[-length(home)],collapse="/")
source(paste0(home,"/utils.R")) # cargo las funciones del paquete

# -------------------------------
# 1 --> Definiciones preliminares
# -------------------------------

option_list <- list(
  make_option(c("-n", "--nodes"), type="character", default=NULL,
              help="Node input information text file name (e.g. 'my_node_data.txt'). The file needs to have the following format:\n\t\tNode35562	k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacteriales;f__Enterobacteriaceae;\n\t\tNode27828	k__Bacteria;p__Proteobacteria\nThe taxonomy must include at least two fields, with the second one being the phylum.", metavar="character"),
  make_option(c("-m", "--medium"), type="character", default="M9",
              help="medium (e.g. 'M9', 'M9[glc]'", metavar="character"),
  make_option(c("--mediadb"), type="character", default=paste0(home,"/my_media.tsv"),
              help="media database file name", metavar="character"),
  make_option(c("-c", "--checking"), type="logical", default=FALSE,
              help="TRUE or FALSE (default). If TRUE, generated and analysed models are limited to those pairs of species that are co-ocurring in pairs in one or more samples of a specified experiment.", metavar="logical"),
  make_option(c("-e", "--experiment"), type="character", default=NULL,
              help="Experiment input information text file name (e.g. 'table.from_biom.tsv'). The tab-separated file needs to have the following format:\n\t\t#OTU ID\tS1\tS2\tS3\t(...)\n\t\tO1\t000000\t000001\t000002\t(...) ", metavar="character"),
  make_option(c("--coupling"), type="logical", default=TRUE,
              help="If TRUE, smetana computes an aditional analysis without the --no-coupling option (see smetana help). Default is TRUE.", metavar="logical"),
  make_option(c("--nucmer"), type="character", default=paste0(home,"/MUMmer3.23/nucmer"),
              help="Path to the nucmer executable (e.g. './MUMmer3.23/nucmer', '~/my_apps/nucmer'.", metavar="character"),
  make_option(c("--showcoords"), type="character", default=paste0(home,"/MUMmer3.23/show-coords"),
              help="Path to the show-coords executable (e.g. './MUMmer3.23/show-coords', '~/my_apps/show_coords'.", metavar="character"),
  make_option(c("--tree"), type="character", default=paste0(home,"/99_otus_nodes.tree"),
              help="16S phylogenetic tree file name. The node names of the trees should be modified, and a genuine name should be given for all. This could be done for example using R with the function 'makeNodeLabel' from 'ape' package. Default is './99_otus_nodes.tree' (Greengenes gg_13_5).", metavar="character"),
  make_option(c("--fasta"), type="character", default=paste0(home,"/99_otus.fasta"),
              help="16S sequences to be analyzed in multifasta format. Default is './99_otus_nodes.tree' (Greengenes gg_13_5).", metavar="character"),
  make_option(c("--db16s"), type="character", default=paste0(home,"/bac120_ssu_reps_r95.fna"),
              help="16S sequences database. Nucmer will align the tree leaves' 16S sequences to this database. Default is './bac120_ssu_reps_r95.fna' (from GTDB https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/genomic_files_reps/).", metavar="character"),
  make_option(c("--dbproteins"), type="character", default=paste0(home,"/protein_faa_reps/bacteria/"),
              help="Aminoacid sequences database. CarveMe will take files from here that correspond to Nucmer hits and create SBML models from those files. Default is 'protein_faa_reps/bacteria/' (from GTDB https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/genomic_files_reps/).", metavar="character"),
  make_option(c("--run_smetana"), type="logical", default=TRUE,
              help="If TRUE, runs a Smetana analysis for inter-node pairs. If FALSE, skips the Smetana analysis. Default: TRUE.", metavar="logical"),
  make_option(c("--cores"), type="numeric", default=4,
              help="Number of cores to use in parallelization processes (mclapply). Default: 4.", metavar="numeric"))

parser <- OptionParser(option_list=option_list)
opt <- parse_args(parser)

t <- read.table(opt$nodes,sep = "\t")
nodos      <- levels(t[[1]])
taxonom    <- levels(t[[2]])
medium     <- opt$medium
mediadb    <- opt$mediadb
checking   <- opt$checking
run_smetana<- opt$run_smetana
coupling   <- opt$coupling

if (checking==TRUE & is.null(opt$experiment)) {
  stop("When --checking TRUE, a experiment file must be specified")
  } else if (checking == TRUE) {
  exp = opt$experiment
  }

nucmer_path     <- opt$nucmer
showcoords      <- opt$showcoords

tree_file         <- opt$tree
otus_fasta_file   <- opt$fasta
db_16S            <- opt$db16s
db_protein_folder <- opt$dbproteins

grampospath <- paste0(home,"/grampos.csv") # estos dos paths no son personalizables (no le veo sentido)
gramnegpath <- paste0(home,"/gramneg.csv")

cores <- opt$cores

# -----------------------
# 2 --> carga de archivos
# -----------------------
if (!require("ape", quietly=TRUE)) BiocManager::install("ape")
tn = ape::read.tree(tree_file)
nodos_hojas <- mclapply(nodos, function(nodo) {ape::extract.clade(tn, nodo)}, mc.cores=cores)

# divido el fasta en dos archivos para facilitar su parsing después:
system(paste0("cat ",otus_fasta_file," | grep '>' > 99_otus_col1"))
system(paste0('cat ',otus_fasta_file,' | grep ">" -v  > 99_otus_col2'))

df1 <- read.csv("99_otus_col1",header = F)
fasta <- read.csv("99_otus_col2",header = F)
rownames(fasta) <- sub(">", "", df1[,1]) # elimino ">"
fasta <- as.data.frame(t(fasta)) # Columna 1: > ID. Columna 2: secuencia

system("rm 99_otus_col*") # elimino archivos temporales

# ------------------------------------------------
# 3 --> alinear cada hoja, filtrar y hacer modelos
# ------------------------------------------------
# Asigno secuencias 16S a cada hoja
if (checking == TRUE) {
  # Pares de hojas de los dos nodos. INTER-NODO.
  filtered_pairs <- check(nodos=nodos_hojas, exp=exp,cores=cores) 
  # De cada hoja que pase el checking, tomo la secuencia de 16S del fasta original.
  # Son todas las hojas de cada nodo que pasan cualquiera de los dos checkings.
  checked_tipl <- list(levels(filtered_pairs[,1]),levels(filtered_pairs[,2]))
  nodos_16S    <- mclapply(checked_tipl, function(n) {fasta[n]},mc.cores=cores)
  } else {
    # De cada hoja (sin checking), cojo la secuencia de 16S del fasta original
    nodos_16S    <- mclapply(nodos_hojas, function(nodo) {fasta[nodo$tip.label]},mc.cores=cores)
  }


system("mkdir models")
# Para cada nodo creo una carpeta de resultados
for (i in c(1:length(nodos))) {
  filepath   <- paste("models/",nodos[i],"/",sep="")
  system(paste("mkdir", filepath)) # la carpeta tendrá el mismo nombre que el nodo

  # Alineo cada hoja de cada nodo con Nucmer y obtengo un genoma adecuado para cada una,
  # seleccionando solo las hojas cuyos hits pasan un filtro de calidad
  nucmer_res_final <- find_alignment_hits(filepath, nodos_16S[[i]], nucmer_path, db_16S, showcoords, cores)

  # Modelado con CarveMe de todas las hojas de cada nodo que pasan el filtro de Nucmer
  print(paste0("Creating models for ",nodos[i],"..."))
  dump <- simplify2array(mclapply(nucmer_res_final, 
                                  FUN = function(line) {carve(line, taxonom[i], filepath, db_protein_folder)},mc.cores=cores))
  print(paste0("Finished modelling for ",nodos[i],"."))
  
  }

# -------------------------------------
# 4 --> análisis metabólico con Smetana
# -------------------------------------

if (run_smetana == TRUE) {
  # Definimos la lista de parejas a analizar
  if (checking == TRUE) {
    pairs <- filtered_pairs
    } else {
      pairs <- expand.grid(nodos_hojas[[1]]$tip.label, nodos_hojas[[2]]$tip.label, KEEP.OUT.ATTRS = F)
    }
  
  # Ejecutamos Smetana para cada pareja inter-nodo de hojas para las que se ha creado
  # un modelo. Esta lista se parejas se guardará en smetana_results/generated_pairs.txt.
  # También guardamos la lista de parejas que no han sido analizadas por Smetana, en el
  # archivo smetana_results/filtered_out_pairs.txt. Incluimos el porcentaje de parejas 
  # que pasan y no pasan el filtro.
  #~~~~~~~~~~~~~~~~~~
  output = "smetana_results/"
  system(paste0("mkdir ",output))
  system(paste0("mkdir ",output,"global"))
  system(paste0("mkdir ",output,"detailed"))
  if (coupling==TRUE) {
    output_coupling = paste0(output,"coupling/")
    system(paste("mkdir",output_coupling))
    system(paste0("mkdir ",output_coupling,"global"))
    system(paste0("mkdir ",output_coupling,"detailed"))
  } else {
    output_coupling = NULL
  }
  #~~~~~~~~~~~~~~~~~
  generated_pairs_filename = "generated_pairs.txt"
  dump <- file.create(paste0(output,generated_pairs_filename)) # vaciamos el archivo, de existir, o lo creamos si no existe
  
  # Paralelización con parApply (solo para este paso)
  failed_pairs <- mcmapply(pairs, FUN=function(z) {
    smetana(z, modelfilepath = "models/", output=output, 
                        coupling=coupling, output_coupling=output_coupling, 
                        generated_pairs_filename=generated_pairs_filename)
    }, mc.cores=cores)
  
  failed_pairs[simplify2array(mclapply(failed_pairs, is.null, mc.cores=cores))] <- NULL
  failed_pairs <- t(failed_pairs) # preparamos la matriz para el siguiente paso
  colnames(failed_pairs) <-  seq(length(failed_pairs[1,]))
  
  # Anotamos las parejas que no han sido analizadas por Smetana (no pasaron el filtro de Nucmer)
  write(x=paste0("filtered out: ", 100*length(failed_pairs)/length(pairs[,1])/2,"%"), file=paste0(output,"filtered_out_pairs.txt"))
  dump <- mclapply(colnames(failed_pairs), FUN=function(col){
    cat(failed_pairs[,col],"\n", file=paste0(output,"filtered_out_pairs.txt"), append=T)
    }, mc.cores=cores)

  print(paste0("Finished SMETANA analysis."))
}
# 
end.time <- Sys.time()
time.taken <- end.time - start.time
print(paste0("Execution time: ",format(time.taken,format = "%H %M %S")))