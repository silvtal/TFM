#!/usr/bin/env Rscript

start.time <- Sys.time() # para devolver al final el tiempo de ejecución
# INPUT: annotate.R -g genomes -o outputname --outputdir outputdir
# INFO: THIS SCRIPT TAKES A LIST OF GENOMES, ANNOTATES ALL OF THEM WITH EGGNOG MAPPER AND RETURNS THE ANNOTATED FILES PLUS A CONSENSUS ONE

# =======================
#        FUNCIONES       
# =======================
library("optparse")
library("parallel")

home <- strsplit(paste0("./",getopt::get_Rscript_filename()),split="/")[[1]]
home <- paste(home[-length(home)],collapse="/")
source(paste0(home,"/utils.R")) # cargo las funciones del paquete

# ==============
#  Definiciones 
# ==============

option_list <- list(
  make_option(c("--skip_alignment"), type="logical", default=TRUE,
              help="If TRUE, a list of genomes to annotate is used as input (--genomes, default). If FALSE, the selected node or nodes' (--nodes) leaves of a given tree (--tree-file, /99_otus_nodes.tree by default) are aligned using Nucmer against a given database. The corresponding protein genome files are the genomes that will be annotated.", metavar="logical"),
  make_option(c("--skip_consensus"), type="logical", default=FALSE,
              help="If TRUE, eggNOG-mapper creates the annotation files but a consensus file is not made. Default: FALSE.", metavar="logical"),
  make_option(c("-g", "--genomes"), type="character", default=NULL,
              help="File containing the list of genomes to annotate. This input is only used when skip_alignment is TRUE (default).", metavar="character"),
  make_option(c("-n", "--nodes"), type="character", default=NULL,
              help="Node input information text file name (e.g. 'my_node_data.txt'). The file needs to have the following format:\n\t\tNode35562	k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacteriales;f__Enterobacteriaceae;\n\t\tNode27828	k__Bacteria;p__Proteobacteria\nThe taxonomy must include at least two fields, with the second one being the phylum.", metavar="character"),
  make_option(c("-m", "--medium"), type="character", default="M9",
              help="medium (e.g. 'M9', 'M9[glc]'", metavar="character"),
  make_option(c("--outputdir"), type="character", default="./annotate_results",
              help="output directory for annotation files and the consensus model(s)", metavar="character"),
  make_option(c("-o","--outputname"), type="character", default="node_consensus",
              help="output name for the consensus model(s)", metavar="character"),
  make_option(c("--mediadb"), type="character", default=paste0(home,"/my_media.tsv"),
              help="media database file name", metavar="character"),
  make_option(c("-c", "--checking"), type="logical", default=FALSE,
              help="TRUE or FALSE (default). If TRUE, generated and analysed models are limited to those pairs of species that are co-ocurring in pairs in one or more samples of a specified experiment.", metavar="logical"),
  make_option(c("-e", "--experiment"), type="character", default=NULL,
              help="Experiment input information text file name (e.g. 'table.from_biom.tsv'). The tab-separated file needs to have the following format:\n\t\t#OTU ID\tS1\tS2\tS3\t(...)\n\t\tO1\t000000\t000001\t000002\t(...) ", metavar="character"),
  make_option(c("--nucmer"), type="character", default=paste0(home,"/MUMmer3.23/nucmer"),
              help="Path to the nucmer executable (e.g. './MUMmer3.23/nucmer', '~/my_apps/nucmer'.", metavar="character"),
  make_option(c("--showcoords"), type="character", default=paste0(home,"/MUMmer3.23/show-coords"),
              help="Path to the show-coords executable (e.g. './MUMmer3.23/show-coords', '~/my_apps/show_coords'.", metavar="character"),
  make_option(c("--emapper_path"), type="character", default=paste0(home,"/eggnog-mapper-master/emapper.py"),
              help="Path to the emapper.py executable file (e.g. './eggnog-mapper-my_version/emapper.py').",metavar="character"),
  make_option(c("--tree"), type="character", default=paste0(home,"/99_otus_nodes.tree"),
              help="16S phylogenetic tree file name. The node names of the trees should be modified, and a genuine name should be given for all. This could be done for example using R with the function 'makeNodeLabel' from 'ape' package. Default is './99_otus_nodes.tree' (Greengenes gg_13_5).", metavar="character"),
  make_option(c("--fasta"), type="character", default=paste0(home,"/99_otus.fasta"),
              help="16S sequences to be analyzed in multifasta format. Default is './99_otus_nodes.tree' (Greengenes gg_13_5).", metavar="character"),
  make_option(c("--db16s"), type="character", default=paste0(home,"/bac120_ssu_reps_r95.fna"),
              help="16S sequences database. Nucmer will align the tree leaves' 16S sequences to this database. Default is './bac120_ssu_reps_r95.fna' (from GTDB https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/genomic_files_reps/).", metavar="character"),
  make_option(c("--dbproteins"), type="character", default=paste0(home,"/protein_faa_reps/bacteria/"),
              help="Aminoacid sequences database. CarveMe will take files from here that correspond to Nucmer hits and create SBML models from those files. Default is 'protein_faa_reps/bacteria/' (from GTDB https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/genomic_files_reps/).", metavar="character"),
  make_option(c("--dmnd_db"), type="character", default=NULL,
              help="Path to a custom Diamond database. Specially useful when using a Diamond version different from the eggNOG-mapper one, such as the one created with create_compatible_database.R.", metavar="character"),
  make_option(c("--cores"), type="numeric", default=4,
              help="Number of cores to use in parallelization processes (mclapply). Default: 4.", metavar="numeric"))

parser <- OptionParser(option_list=option_list)
opt <- parse_args(parser)

cores <- opt$cores
skip_alignment <- opt$skip_alignment 
if (skip_alignment == FALSE) {
  stop("The Nucmer alignment functionality is not yet available.")
  t <- read.table(opt$nodes,sep = "\t")
  nodos      <- levels(t[[1]])
  taxonom    <- levels(t[[2]])
  nucmer_path<- opt$nucmer
  showcoords <- opt$showcoords
  tree_file         <- opt$tree
  otus_fasta_file   <- opt$fasta
  db_16S            <- opt$db16s
  checking   <- opt$checking
  if (checking==TRUE & is.null(opt$experiment)) {
    stop("When --checking TRUE, a experiment file must be specified")
    } else if (checking == TRUE) {
      exp = opt$experiment
    }
} else {
  if (is.null(opt$genomes)) {
    stop("Please provide a file with the list of genomes to annotate")
  } else {
    genomes <- opt$genomes # TODO como puedo hacer que sea capaz de coger una lista, como por ej poniendo "Nodo*" y que coja varias carpetas?
    }
  }

outputdir         <- opt$outputdir
outputname        <- opt$outputname
db_protein_folder <- opt$dbproteins
emapper_path      <- opt$emapper_path

if (is.null(opt$dmnd_db)) {
  emapper_folder  <- system(paste0("echo $(dirname ",emapper_path,")"),intern=TRUE)
  dmnd_db         <- paste0(emapper_folder,"/data/eggnog_proteins.dmnd")
} else
  dmnd_db         <- opt$dmnd_db

# =======================
#   Anotar los genomas   
# =======================

if (skip_alignment == FALSE) {
  genomes <- scan(genomes,character(),sep="\n")
  
  #   Anotar el genoma asignado a cada hoja con eggNOG-mapper   
  # ------------------------------------------------------------
  annotate(genomes, outputdir, db_protein_folder, emapper_path)
  system("rm outputdir/*seed_orthologs")
  
  #   Crear genoma consenso   
  # --------------------------
  if (!skip_consensus) {
    system(paste("python3 consenso_EGG.py $(realpath",outputdir,") $(realpath",outputdir,")",outputname))
  }
  
  #  Fabricar modelo con CarveMe 
  # -----------------------------
  # TODO en este caso no podemos leer la taxonomía, pero podríamos pedirla al user o permitirle que indique gram u otros argumentos de carveme...
  system(paste0("carve --egg ",outputname,".tsv -o ",outputdir,"/",outputname,".xml"))
  
# TODO reconsider esto de abajo: meto esta funcionalidad?
} else {
  
  #   Obtener secuencias 16S para todas las hojas de ambos nodos  
  # --------------------------------------------------------------
  if (checking == TRUE) {
    filtered_pairs <- check(nodos=nodos_hojas, exp=exp)
    checked_tipl <- list(levels(filtered_pairs[,1]),levels(filtered_pairs[,2]))
    nodos_16S    <- mclapply(checked_tipl, function(n) {fasta[n]},mc.cores=cores)
  } else {
    nodos_16S    <- mclapply(nodos_hojas, function(nodo) {fasta[nodo$tip.label]},mc.cores=cores)
  }
  
  system("mkdir models")
  # Para cada nodo haremos alineamientos, anotaciones y un consenso
  
  genomes=list()
  for (i in c(1:length(nodos))) {
    filepath   <- paste("models/",nodos[i],"/",sep="")
    system(paste("mkdir", filepath))

    #    Alineamiento con Nucmer    
    # ------------------------------
    nucmer_res_final <- find_alignment_hits(filepath, nodos_16S[[i]], nucmer_path, db_16S, showcoords, cores)
    genomes[[i]] <- scan(file=paste0(filepath,"genomes"),what=character(),sep="\n")
    
    #   Anotar el genoma asignado a cada hoja con eggNOG-mapper   
    # ------------------------------------------------------------
    node_outputdir  <- paste0(outputdir,"/",nodos[i])
    node_outputname <- paste0(outputname,"_",nodos[i])
    system(paste0("mkdir ",node_outputdir))
    
    annotate(genomes[[i]], node_outputdir, db_protein_folder, emapper_path)
    system("rm outputdir/*seed_orthologs")
    
    #   Crear genoma consenso   
    # --------------------------
    if (!skip_consensus) {
      system(paste("python3 consenso_EGG.py $(realpath",node_outputdir,") $(realpath",node_outputdir,")",node_outputname))
    }
    
    #  Fabricar modelo con CarveMe 
    # -----------------------------
    system(paste0("carve --egg ",node_outputname,".tsv ",gram(taxonom[i])," -o ",node_outputdir,"/",node_outputname,".xml"))
  }
}


end.time <- Sys.time()
time.taken <- end.time - start.time
print(paste0("Execution time: ",format(time.taken,format = "%H %M %S")))