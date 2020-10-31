suppressMessages(library(tidyverse))
suppressMessages(library(optparse))
suppressMessages(library(R.utils))

option_list <- list(
  make_option(c("-w", "--workDir"), type = "character", default=FALSE,
              help="Root location of all data and codes (workDir/anz)"),
  make_option(c("-c", "--cancer"), type = "character", default=FALSE,
              help="cancer type"),
  make_option(c("-d", "--WhippletIndexjls"), type = "character", default=FALSE,
              help="Whipplet_index_jls location"),
  make_option(c("-i","--index"),type="integer") # lizx: added this option
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

root_dir <- opt$workDir
if(root_dir[length(root_dir)] != "/")
  root_dir <- paste(root_dir,"/",sep = "")
cancertype <- opt$cancer
Whipplet_index_jls <- opt$WhippletIndexjls

manifest_all <- read.table(paste(root_dir,"anz/data/manifest/TCGA_all_sample_info.txt",sep = ""),
                           sep = "\t",header = T,check.names = F,quote = "",stringsAsFactors = F)
manifest <- manifest_all[manifest_all$cases == cancertype,]

sample_number <- nrow(manifest)

stopifnot(1<=opt$index,opt$index<=sample_number)


bulid_WhippletIndex <- function(root_dir,cancertype,Whipplet_index_jls,bam_file,Whippet_index_outdir){
  
  gtf <- paste(root_dir,"anz/ensembl_v97/Homo_sapiens.GRCh38.97.chr.gtf",sep = "")
  reference_fasta <- paste(root_dir,"anz/ensembl_v97/Homo_sapiens.GRCh38.dna.primary_assembly.fa",sep = "")
  

  print("Start building Whipplet Index..............\n")

  ##index
  stopifnot(system(paste("samtools index -@ 8 ",bam_file,sep = ""),wait = T) == 0)
  
  ###########
  ##Build Whipplet Index for current cancertype
  stopifnot(system(paste("julia ",Whipplet_index_jls, 
                         " --fasta ",reference_fasta,
                         " --gtf ",gtf,
                         " --bam ",bam_file,
                         " -x ", Whippet_index_outdir,
                         " --bam-both-novel ",
                         " --bam-min-reads 5",
                         sep = ""),wait = TRUE) == 0)
  
  ##delect bam_file
  # lizx: don't remove files, I can delete myself
  # file.remove(bam_file) 
  
  print("Finish building Whipplet Index!!!!!!!!\n")
  
}


# lizx: removed parallel for loop. changed every index j to opt$index

bam_file_name <- paste(root_dir,"anz/result/",cancertype,"/05.STAR_2ndAlign/",manifest$associated_entities[opt$index],".rmdup.bam",sep = "")
whippet_index_rootdir <- paste(root_dir,"anz/result/",cancertype,"/06.buildWhippletIndex/",sep = "")
Whippet_index_outdir <- paste(whippet_index_rootdir,manifest$associated_entities[opt$index],sep = "")
if(!file.exists(whippet_index_rootdir)){
  dir.create(whippet_index_rootdir,recursive = T)
}
bulid_WhippletIndex(root_dir,cancertype,Whipplet_index_jls,bam_file_name,Whippet_index_outdir)
        
