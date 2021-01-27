#!/usr/bin/env Rscript

# ----------------------------------------------------------
# Crear una base de datos compatible con el diamond del PATH
# ----------------------------------------------------------
library("optparse")

option_list <- list(
  make_option(c("-o","--outputdir"), type="character", default=".",
              help="Output directory for the new database.", metavar="character"),  
  make_option(c("-e","--emapper_path"), type="character", default=NULL,
              help="EggNOG-mapper directory. Contains the old database and the old Diamond version.", metavar="character"),  
  make_option(c("-d","--diamond"), type="character", default="diamond",
              help="New Diamond version location. Default is 'diamond' (PATH).", metavar="character"))

parser <- OptionParser(option_list=option_list)
opt <- parse_args(parser)

outputdir <- opt$outputdir
diamond   <- opt$diamond

if (is.null(opt$emapper_path)) {
  stop("Please specify the eggNOG-mapper path")
} else {
  emapper_path <- opt$emapper_path
}
old_db       <- paste0(emapper_path,"/data/eggnog_proteins.dmnd")
old_dmnd     <- paste0(emapper_path,"/bin/diamond")

make_compatible_database <- function(emapper_path) {
  print(paste0(old_dmnd," getseq --db ",old_db," > '",outputdir,"/eggnog_proteins.faa';
              ",diamond," makedb --in '",outputdir,"/eggnog_proteins.faa' --db ",outputdir,"/eggnog_proteins_compatible"))
}

make_compatible_database(emapper_path)

print(paste0("New Diamond database saved as: ",outputdir,"/eggnog_proteins_compatible.dmnd"))