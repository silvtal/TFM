#!/usr/bin/env Rscript

library("optparse")
library("parallel")

# ===================
#     Take input     
# ===================
### HELP ###
# parser.R -d/-g input_file -c cores --report_name reportname(can include path)

option_list <- list(
  make_option(c("-g", "--global"), type="character", default=NULL,
              help="Input folder with Smetana global results (e.g. ./my_results/NodeXXXXX/smetana_results/global).", metavar="character"),
  make_option(c("-d", "--detailed"), type="character", default=NULL,
              help="Input folder with Smetana global results (e.g. ./my_results/NodeXXXXX/smetana_results/global).", metavar="character"),
  make_option(c("-c","--cores"), type="numeric", default=4,
              help="Number of cores to use in parallelization", metavar = "character"),
  make_option(c("--report_name"), type="character", default="report.txt",
              help="Name of the output report file. Default: ./report.txt", metavar="character")
)

parser <- OptionParser(option_list=option_list)
opt <- parse_args(parser)

global      <- opt$global
detailed    <- opt$detailed
cores       <- opt$cores
report_name <- opt$report_name

if (is.null(global)) {
  if (is.null(detailed)) {
    stop("Please specify an input folder with --global or --detailed.")
  } else {
    is.detailed <- TRUE
    is.global   <- FALSE
  }
} else if (is.null(detailed)) {
  is.global   <- TRUE
  is.detailed <- FALSE
} else {
  stop("Only --global or --detailed input is accepted. Please choose only one.")
}

if (is.global){
  # Open the dataset
  # ================
  filenames  <- list.files(global, pattern="*.tsv", full.names=TRUE)
  
  dataset    <- do.call("rbind", mclapply(filenames, FUN=function(filename)  {
      df <- read.table(filename, sep="\t", header=TRUE, na.strings = "n/a")
      rownames(df) <- tail(strsplit(filename,split="/")[[1]],n=1) # TODO fix this
      return(df)
    } ,mc.cores=cores))
  
  # Report MRO=n/a files
  # ====================
  mro.is.na <- as.logical(is.na(dataset["mro"]))
  row_names  <- list.files(global, pattern="*.tsv", full.names=FALSE)
  not_growing <- row_names[mro.is.na]
  
  #  MRO values  
  # =============
  growing <- dataset[!mro.is.na,]
  
  if (length(growing)==0) {
    write("MRO was n/a for all files. The co-culture can't grow.",file=report_name)
    stop ("MRO was n/a for all files. The co-culture can't grow.")
  } else {
    mean_mro <- mean(growing[,"mro"])
    min_mro  <-  min(growing[,"mro"])
    max_mro  <-  max(growing[,"mro"])
  }
  
  # MIP report
  # ==========
  mip.is.na <- as.logical(is.na(growing[,"mip"]))
  mip.is.not.na <- growing[!mip.is.na,]
  
  mean_mip <- mean(mip.is.not.na[,"mip"])
  min_mip  <-  min(mip.is.not.na[,"mip"])
  max_mip  <-  max(mip.is.not.na[,"mip"])
  
  # Print report
  # ============
  # TODO open only once?
  write("The following files have n/a MRO and MIP, which means they can't grow by themselves:\n",file=report_name)
  write(not_growing,file=report_name,append=TRUE)
  
  write("\n\nMRO (metabolic resource overlap) calculates how much the species compete for the same metabolites", file=report_name, append=TRUE)
  write(paste("\nMean MRO:",   mean_mro, sep="\n"), file=report_name, append=TRUE)
  write(paste("\nMinimum MRO:", min_mro, sep="\n"), file=report_name, append=TRUE)
  write(paste("\nMaximum MRO:", max_mro, sep="\n"), file=report_name, append=TRUE)
  
  write("\n\nMIP (metabolic interaction potential) calculates how many metabolites the species can share to decrease their dependency on external resources", file=report_name, append=TRUE)
  write(paste("\nMean MIP:",   mean_mip, sep="\n"), file=report_name, append=TRUE)
  write(paste("\nMinimum MIP:", min_mip, sep="\n"), file=report_name, append=TRUE)
  write(paste("\nMaximum MIP:", max_mip, sep="\n"), file=report_name, append=TRUE) 
  
  write("\n\nFiles where MRO values are available:", file=report_name, append=TRUE)
  write.table(growing, file=report_name, sep="\t", quote=FALSE, append=TRUE)
  
  
} else {
  # Open the data files
  # ===================
  filenames   <- list.files(detailed, pattern="*.tsv", full.names=TRUE)
  
  files_text  <- mclapply(filenames, FUN=function(filename)  {
    read.csv(filename, sep="\t", header=TRUE)
  } ,mc.cores=cores)
  
  files_names <- list.files(detailed, pattern="*.tsv", full.names=FALSE)
  names(files_text) <- files_names
  
  
  # Report empty files
  # ==================
  is.empty  <- as.numeric(mclapply(files_text,FUN=function(x){length(x[,1])},mc.cores=cores))==0
  
  empty     <- files_names[is.empty]
  not.empty <- files_names[!is.empty]
  
  # Analyze metabolites
  # ===================
  merged.files    <- Reduce(f = function(...) merge(...,all=T), files_text)
  rm(files_text) # free memory
  
  # Create node columns
  check_node      <- function(x,node_1_index,tags=c("first","second")){
    if (x %in% node_1_index) {result <- tags[1]} else {result <- tags[2]}
    return(result)
    }
  
  node_1_index    <- unique(sapply(files_names,FUN=function(filename){head(strsplit(filename,split="_")[[1]],n=1)},USE.NAMES = FALSE))
  
  donor_node      <- sapply(merged.files[,"donor"], check_node, node_1_index, USE.NAMES=FALSE)
  receiver_node   <- sapply(merged.files[,"receiver"],check_node, node_1_index, USE.NAMES=FALSE)

  merged.files    <- cbind(merged.files, donor_node, receiver_node)
  
  # Exchanged in general
  all.metabolites <- setNames(aggregate(x=merged.files[,c(7,9)],by=list(merged.files$compound),mean),
                              nm=c("compound","mus","smetana"))
  
  # Exchanged by node
  all.metabolites.by.node <- setNames(aggregate(x=merged.files[,c(7,9)], 
                               by=list(merged.files$compound,merged.files$donor_node, merged.files$receiver_node),
                               mean,),nm=c("compound","donor_node","receiver_node","mus","smetana"))
  
  # Ordered
  all.metabolites <- all.metabolites[order(all.metabolites$smetana,decreasing=TRUE),]
  all.metabolites.by.node <- all.metabolites.by.node[order(all.metabolites.by.node$smetana,decreasing=TRUE),]
  
  # Print report
  # ============
  # TODO open only once?
  write("The following files are empty, which means no exchange between species:\n",file=report_name)
  write(empty,file=report_name,append=TRUE)
  write(paste0("\nThere are ",length(not.empty)," files which are not empty:"),file=report_name,append=TRUE)
  write(not.empty,file=report_name,append=TRUE)
  
  write("\n10 most exchanged metabolites:",file=report_name,append=TRUE)
  write.table(all.metabolites[1:10,],file=report_name,append=TRUE,quote=FALSE,sep="\t",row.names=FALSE)
  write("\n10 most exchanged metabolites by donor node:",file=report_name,append=TRUE)
  write.table(all.metabolites.by.node[1:10,],file=report_name,append=TRUE,quote=FALSE,sep="\t",row.names=FALSE)
  write("\n'first' node includes ",file=report_name,append=TRUE)
  write(node_1_index, file=report_name,append=TRUE)
  write("\nThe whole dataset is as follows:",file=report_name,append=TRUE)
  write.table(merged.files,file=report_name,append=TRUE,quote=FALSE,sep="\t",row.names=FALSE)
}
