#!/usr/bin/env Rscript
library(parallel)

rndnum <- runif(n = 1)*10000000 # random number for tempfiles. Helps avoid superposition when there are multiple parallel processes

nucmer <- function(nucmer_path, db_16S, fasta_16S, prefix="out") {
  # Alinea cada secuencia del archivo fasta indicado ("fasta_16S") con las secuencias
  # de la base de datos ("db_16S") y devuelve todos los hits en un archivo llamado
  # "./out.delta".
  system(paste(nucmer_path, db_16S, fasta_16S, "-p", prefix))
}

check = function(nodos, abuntable, cutoff = 0, cores=1) {
  # Dada una pareja de nodos ("nodos", formato ape::phylo) y una tabla de 
  # abundancias, determina y devuelve las combinaciones de OTUs inter-nodo que
  # están presentes en una misma muestra experimental.
  
  if (length(nodos) > 2) {
    warning(paste("The node name list is longer than 2. ", length(nodos)," Only the first two nodes will be checked for coexistence"))
    nodos <- nodos[1:2]
  }
  # Cutoff: minimum abundance necessary for considering an OTU to be "present"
  # in a given sample
  
  # Leo abundancias
  abuntable <- read.csv(abuntable, sep="\t",skip = 1,row.names=1)
  abuntable <- abuntable[1:(ncol(abuntable)-1)] # exclude taxonomy column
  
  # Seleccionamos solo las hojas que estaban en el experimento
  lista_nodos <- mclapply(nodos, function(nodo) { nodo$tip.label[nodo$tip.label %in% rownames(abuntable)]},mc.cores=cores)
  
  # Obtengo todas las combinaciones de hojas para cada nodo
  pairs = expand.grid(lista_nodos[[1]], lista_nodos[[2]], KEEP.OUT.ATTRS = F)
  
  # Selecciono las parejas presentes en una misma muestra
  pairs[,3] = vector(mode="logical", length=length(pairs[,2])) # indicador de si coinciden o no
  for (muestra in colnames(abuntable)){
    # anoto qué otus hay en cada muestra
    otus=rownames(subset(abuntable[muestra], abuntable[muestra]>cutoff))
    # indico para cada pareja si coinciden en esa muestra o en otra ya comprobada
    pairs[,3] = pairs[,3] | pairs[,1] %in% otus & pairs[,2] %in% otus
  }
  filtered_pairs  = subset(pairs, pairs[,3]==TRUE)[,0:2]
  return(filtered_pairs)
}

check_intra = function(nodos, abuntable) {
  # Dada una pareja de nodos ("nodos", formato ape::phylo) y la ruta de un archivo con
  # datos experimentales ("abuntable"), determina y devuelve las combinaciones de OTUs 
  # intra-nodo que coinciden en una misma muestra experimental
  abuntable = read.csv(abuntable,sep="\t",skip = 1,row.names=1)
  abuntable <- abuntable[1:(ncol(abuntable)-1)] # exclude taxonomy column
  
  # Seleccionamos solo las hojas que estaban en el experimento, para agilizar el proceso
  lista_nodos <- lapply(nodos, function(nodo) { nodo$tip.label[nodo$tip.label %in% rownames(abuntable)]} )
  
  # Obtengo todas las combinaciones de hojas para cada nodo
  pairs1 = expand.grid(lista_nodos[[1]], lista_nodos[[1]], KEEP.OUT.ATTRS = F)
  pairs2 = expand.grid(lista_nodos[[2]], lista_nodos[[2]], KEEP.OUT.ATTRS = F)
  
  filtered_pairs  = list(pairs1,pairs2)
  return(filtered_pairs)
}



gram  = function(taxonom, mediadb, media, gramneg = "gramneg.csv", grampos = "grampos.csv") {
  # Dado un string con la taxonomía en formato GreenGenes, devuelve si es 
  # gram-positiva o gram-negativa. Si no se reconoce como G+ ni G-, devuelve
  # un string vacío.
  phylum  <- strsplit(taxonom, split=";")[[1]][2]
  gramneg <- levels(read.table(gramnegpath)[[1]])
  grampos <- levels(read.table(grampospath)[[1]])
  
  if (phylum %in% gramneg) { #FIXTHIS
    result <- paste0("--verbose -u gramneg -g ", media, " --mediadb ", mediadb)
  }
  else if (phylum %in% grampos) {
    result <- paste0("--verbose -u grampos -g ", media, " --mediadb ", mediadb)
  }
  else {
    result <- ""
  }
  return(result)
}




carve = function(line, taxonom, mediadb, media, outputpath, db_protein_folder) {
  # Dado un string del tipo "<ID de hoja> <ID de genoma anotado>", crea un
  # modelo metabólico de dicha hoja a partir del archivo correspondiente al
  # genoma dado.
  leaf = strsplit(line, split=" ")[[1]][1]
  file = strsplit(line, split=" ")[[1]][2]
  outf = paste0(outputpath,leaf,".xml")
  if (file.exists(outf)) {
    print(paste0(outf," model already exists. Moving to the next one..."))
  } else {
    system(paste0("carve ",db_protein_folder,file,"_protein.faa ",
                  gram(taxonom, mediadb, media)," -o ",outf)) # nombres de archivo = nombre de hoja
    print(paste0(outf," model created"))
  }
}



smetana = function(pair, nodos, modelfilepath="models/", output, coupling=TRUE, output_coupling=NULL, generated_pairs_filename="generated_
pairs.txt") {
  # Dado un string del tipo <hoja1> <hoja2>, lanza smetana para las hojas dadas,
  # si existen sus archivos. Si no existe alguno de los dos archivos o ninguno, 
  # devuelve sus nombres. Si coupling==TRUE, se computan también los análisis sin
  # la opción de smetana "--no-coupling".
  # File creation
  if (coupling==TRUE & is.null(output_coupling)) {
    stop("A output name for the coupling results must be specified.")
  }
  filepath1 = paste0(modelfilepath,nodos[1],"/") 
  filepath2 = paste0(modelfilepath,nodos[2],"/") 
  m1 = pair[[1]]
  m2 = pair[[2]]
  if (length(pair)>2) {warning(paste("The node name list is longer than 2. Only the first two nodes will be analysed with Smetana (", m1, ", ", m2, ")"))}
  if (file_test("-f",paste0(filepath1,m1,".xml")) & file_test("-f",paste0(filepath2,m2,".xml"))) {
    for (i in c("global","detailed")){
      output_filename=paste0(output,i,"/",m1,"_",m2,"_",medium)
      if (!file.exists(paste0(output_filename,"_",i,".tsv"))) { #solo crea el archivo si no existía ya
        system(paste0("smetana --",i," ",filepath1,m1,".xml"," ",filepath2,m2,".xml"," --flavor bigg -m ",medium,
                      " --mediadb ",mediadb," --molweight --no-coupling -o ",output_filename))
        write(paste(m1,m2),file=paste0(output, generated_pairs_filename),append=TRUE) # lista de parejas que se han analizado con Smetana
      } else {
        print(paste0(output_filename,"_",i,".tsv already exists. Moving to the next pair..."))
      }
    }
    if (coupling==TRUE) {  # se hace una ejecución más, pero solo para detailed.
      output_filename_c=paste0(output_coupling,i,"/",m1,"_",m2,"_",medium)
      if (!file.exists(paste0(output_filename,"_",i,".tsv"))) { #solo crea el archivo si no existía ya
        system(paste0("smetana --detailed"," ",filepath1,m1,".xml"," ",filepath2,m2,".xml"," --flavor bigg -m ",medium,
                      " --mediadb ",mediadb," --molweight -o ",output_filename_c))
      } else {
        print(paste0(output_filename_c,"_detailed.tsv already exists. Moving to the next pair..."))
      }
    }
  } else {
    print("estamos aqui!!!!!!!!!!!!!!!")
    return(pair) # se devuelve 
  }
}





find_alignment_hits = function(filepath, node_16S, nucmer_path, db_16S, showcoords, cores) {
  # Runs nucmer on each of the leaves of a given node (input is 16S sequences).        
  # Creates 4 different files:                                                         
  #    - nucmer_unfiltered: table including all the best hits for every leaf before the
  #                         97% identity and 90% coverage filter.                      
  #    - nucmer_filtered: table including the best hits with over 97% identity and 90% 
  #                       coverage. Some leaves may be filtered out at this point.     
  #    - queries_v_hits: list of every leaf that passed the filter vs its best database
  #                      hit according to Nucmer. It's a two-column table              
  #    - genomes: list of every database hit that passed the filter. One-column table. 
  #                                                                                    
  # Returns a character vector with the leaf names and hits (the same list of the      
  # queries_v_hits file).                                                              
  
  # Creo un archivo temporal en formato fasta que contenga las secuencias de ese nodo
  write(
    simplify2array(mclapply(colnames(node_16S), FUN=function(hoja){
      paste(paste(">",hoja,sep=""),as.character(node_16S[hoja][,]),sep="\n")}, mc.cores=cores)),
    file=paste0("sec_temp_", rndnum, ".fasta"))
  
  # Alineo cada una de sus hojas con la base de datos
  nucmer(nucmer_path = nucmer_path, 
         db_16S =db_16S,
         fasta_16S = paste0("sec_temp_", rndnum, ".fasta"),
         prefix = rndnum)
  
  system(paste0("mv ", rndnum, ".delta ", filepath))
  
  # Descarto hits con identidad por debajo de 97% con show-coords
  nucmer_res <- system(paste(showcoords," -c -l -I 97 ",filepath,
                             paste0(rndnum, ".delta"), sep=""), intern = TRUE)
  # -c  Include percent coverage information
  # -l  Include the sequence length information
  # -I float    Set minimum percent identity to display
  
  # Creo archivo de mejores resultados sin filtrar por coverage
  header = nucmer_res[4:5]
  nucmer_res = gsub("|", "", nucmer_res[-(0:5)], fixed=T) # elimino separador
  nucmer_res = gsub("(?<=[\\s])\\s*|^\\s+|\\s+$", "", nucmer_res, perl=TRUE) # fusiono espacios
  
  write(header, paste(filepath,"nucmer_unfiltered",sep=""))
  write(nucmer_res, paste(filepath,"nucmer_unfiltered",sep=""),append = TRUE) # guardo
  
  # Creo archivo de mejores resultados filtrados por %ID y coverage y lo guardo
  write(header, paste(filepath,"nucmer_filtered",sep=""))
  
  system(paste("awk '{ if ($10>=90) { print } }' ", # filtro por cobertura
               filepath,"nucmer_unfiltered > ",
               filepath,"nucmer_temp1",sep=""))
  system(paste("sort -r -k13,13 -k7,7 ", # ordeno por query y por %ID
               filepath,"nucmer_temp1 > ",
               filepath,"nucmer_temp2",sep=""))
  system(paste("awk -F ' ' -v q='' '{if ($13!=q) {q=$13; print}}' ", # escojo los mejores resultados
               filepath,"nucmer_temp2 >>",
               filepath,"nucmer_filtered",sep="")) #guardo
  
  # Lista filtrada de hits a modelar
  queries_v_hits = system(paste("awk -F ' ' -v q='' '{if ($13!=q) {q=$13; print $13,$12}}' ",
                                filepath,"nucmer_temp2",sep=""), intern = TRUE)
  
  # guardo parejas como archivo (a modo de índice informativo)
  write(queries_v_hits, file=paste0(filepath, "queries_v_hits"))
  
  # guardo solo los genomas en un archivo (útil para anotar posteriormente por ejemplo)
  system(paste("awk -F ' ' -v q='' '{if ($13!=q) {q=$13; print $12}}' ",
               filepath,"nucmer_temp2 > ",filepath,"genomes",sep=""), intern = TRUE)
  
  # Borrado de archivos temporales
  system(paste0("rm sec_temp_", rndnum, ".fasta ", filepath,"nucmer_temp*",sep=""))
  
  return(queries_v_hits)
} 



emapper = function(input_fa, db_protein_folder, outputname, outputdir, emapper_path, cores) {
  system(paste0(emapper_path," -m diamond --cpu ",cores," --no_annot --no_file_comments -i ",
                db_protein_folder, input_fa,"_protein.faa"," -o ",outputname," --output_dir ",
                outputdir," --temp_dir /dev/shm --dmnd_db $PWD/",dmnd_db," --override"))
  returned <- system(paste0(emapper_path," --annotate_hits_table ",outputdir,"/",outputname,
                            ".emapper.seed_orthologs -o ",outputname," --output_dir ",outputdir," --override"))
  if (returned != 0) {
    stop("eggNOG-mapper returned non-zero status. If your Diamond version is different from the eggNOG-mapper one 
         (e.g. the CarveMe version is in the PATH) please create a new database with create_compatible_database.R 
         and select it with --dmnd_db when next running annotate.R")
  }
}



annotate = function(genomes, outputdir, db_protein_folder, emapper_path, cores) {
  if (!file.exists(outputdir)){
    system(paste0("mkdir ",outputdir)) }
  
  N = length(genomes) # número total de genomas
  dump <- simplify2array(mclapply(genomes,FUN=function(genome) {
    print(paste0("Annotating genome ",genome," (total: ",N,")"))
    emapper(input_fa=genome,
            db_protein_folder = db_protein_folder,
            outputname=genome,
            outputdir=outputdir,
            emapper_path=emapper_path,
            cores=cores) 
  },mc.cores=cores))
}