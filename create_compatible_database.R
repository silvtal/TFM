# ----------------------------------------------------------
# Crear una base de datos compatible con el diamond del PATH
# ----------------------------------------------------------
library("optparse")

option_list <- list(
  make_option(c("-o","--outputdir"), type="character", default=".",
              help="Output directory for the new database.", metavar="character"),
  make_option(c("-d","--diamond"), type="character", default="diamond",
              help="New Diamond version location. Default is 'diamond' (PATH).", metavar="character"))

parser <- OptionParser(option_list=option_list)
opt <- parse_args(parser)

outputdir <- opt$outputdir
diamond   <- opt$diamond

home         <- strsplit(paste0("./",getopt::get_Rscript_filename()),split="/")[[1]]
emapper_path <- "eggnog-mapper-master" # TODO quizá cambiar esto cuando 
old_db       <- paste0(home,"/",emapper_path,"/data/eggnog_proteins.dmnd")
old_dmnd     <- paste0(home,"/",emapper_path,"/bin/diamond")


make_compatible_database <- function(emapper_path) {
  system(paste0("$(echo ",old_dmnd," getseq --db ",old_db," > ",home,"/",emapper_path,"/data/eggnog_proteins.faa;
              ",diamond,"makedb --in ",home,"/",emapper_path,"/data/eggnog_proteins.faa --db ",outputdir,"/eggnog_proteins_compatible"))
}

print(paste0("New Diamond database saved as: ",outputdir,"/eggnog_proteins_compatible.dmnd"))