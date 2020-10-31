suppressMessages(library(tidyverse))
suppressMessages(library(foreach))
suppressMessages(library(doParallel))
suppressMessages(library(optparse))
suppressMessages(library(R.utils))

option_list <- list(
  make_option(c("-w", "--workDir"), type = "character", default=FALSE,
              help="Root location of all data and codes (workDir/anz)"),
  make_option(c("-n", "--starIndexThread"), type = "integer", default=40,
              help="Threads of STAR building index")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

root_dir <- opt$workDir
starIndexThread <- opt$starIndexThread
if(root_dir[length(root_dir)] != "/")
  root_dir <- paste(root_dir,"/",sep = "")



STAR_1stIndex <- function(root_dir,starIndexThread){
  STAR_Firstindex <- paste(root_dir,"anz/STAR_Firstindex",sep = "")
  gtf <- paste(root_dir,"anz/ensembl_v97/Homo_sapiens.GRCh38.97.chr.gtf",sep = "")
  reference_fasta <- paste(root_dir,"anz/ensembl_v97/Homo_sapiens.GRCh38.dna.primary_assembly.fa",sep = "")
  
  if(!file.exists(STAR_Firstindex)){
    dir.create(STAR_Firstindex,recursive = T)
    ####build STAR First index
    print("*******Starting building STAR first Index...............\n")
    stopifnot(system(paste("STAR --runThreadN ",starIndexThread,
                           " --runMode genomeGenerate", 
                           " --genomeDir ",STAR_Firstindex,
                           " --genomeFastaFiles ",reference_fasta,
                           " --sjdbGTFfile ",gtf,
                           sep = "")
                     ,wait = T) == 0)
    print("*******Finish building STAR first Index!!!\n")
  }else{
    print("STAR 1st index has already been indexed!!!! Thank you, next.....")
  }
}



STAR_1stIndex(root_dir,starIndexThread)