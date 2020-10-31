
suppressMessages(library(optparse))
suppressMessages(library(R.utils))

option_list <- list(
  make_option(c("-w", "--workDir"), type = "character", default=FALSE,
              help="Root location of all data and codes (workDir/anz)")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
root_dir <- opt$workDir
if(root_dir[length(root_dir)] != "/")
  root_dir <- paste(root_dir,"/",sep = "")

##download ensembl v97 reference genome and annotation files
download_Ensembl <- function(root_dir){
  gtf <- paste(root_dir,"anz/ensembl_v97/Homo_sapiens.GRCh38.97.chr.gtf",sep = "")
  reference_fasta <- paste(root_dir,"anz/ensembl_v97/Homo_sapiens.GRCh38.dna.primary_assembly.fa",sep = "")
  print(root_dir)
  if(!file.exists(gtf) & !file.exists(reference_fasta)){
    print("*******Starting downloading Ensembl v97 reference genome.........\n")
    ###download
    stopifnot(system(paste("wget -P ",root_dir,"anz/ensembl_v97/ ftp://ftp.ensembl.org/pub/release-97/gtf/homo_sapiens/Homo_sapiens.GRCh38.97.chr.gtf.gz",sep = ""),wait = T) == 0)
    gunzip(paste(gtf,".gz",sep = ""),remove = TRUE)
    stopifnot(system(paste("wget -P ",root_dir,"anz/ensembl_v97/ ftp://ftp.ensembl.org/pub/release-97/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",sep = ""),wait = T) == 0)
    gunzip(paste(reference_fasta,".gz",sep = ""),remove = TRUE)
    print("*******Finish downloading Ensembl v97 reference genome!!!!\n")
  }else{
    print("Ensembl v97 has already been downloaded!!!!")
  }
}

download_Ensembl(root_dir)


