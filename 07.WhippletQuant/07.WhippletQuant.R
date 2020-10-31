suppressMessages(library(tidyverse))
suppressMessages(library(optparse))
suppressMessages(library(R.utils))

option_list <- list(
  make_option(c("-w", "--workDir"), type = "character", default=FALSE,
              help="Root location of all data and codes (workDir/anz)"),
  make_option(c("-c", "--cancer"), type = "character", default=FALSE,
              help="cancer type"),
  make_option(c("-d", "--WhippletQuantjls"), type = "character", default=FALSE,
              help="Whipplet_quant_jls location"),
  make_option(c("-i", "--index"),type="integer")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

root_dir <- opt$workDir
if(root_dir[length(root_dir)] != "/")
  root_dir <- paste(root_dir,"/",sep = "")
cancertype <- opt$cancer
Whipplet_quant_jls <- opt$WhippletQuantjls
manifest_all <- read.table(paste(root_dir,"anz/data/manifest/TCGA_all_sample_info.txt",sep = ""),
                           sep = "\t",header = T,check.names = F,quote = "",stringsAsFactors = F)
manifest <- manifest_all[manifest_all$cases == cancertype,]

sample_number <- nrow(manifest)
stopifnot(1<=opt$index,opt$index<=sample_number)

runWhipplet <- function(root_dir,Whipplet_index_jls,Whipplet_quant_jls,fastq1_dir, fastq2_dir = "",prefix){
  ##create dir `
  
  ##run whipplet
  julia_command <- paste("julia ",Whipplet_quant_jls, " ", fastq1_dir, " ", fastq2_dir,
                         " -o ",prefix,
                         " -x ",Whipplet_index_jls,
                         " --circ",
                         " --biascorrect",sep = "")

  # lizx: do not delete files

  # delete_file1_command <- paste("rm -f ",fastq1_dir,sep = "")
  # if(fastq2_dir != ""){
  #   delete_file2_command <- paste("rm -f ",fastq2_dir,sep = "")
  #   all_command <- paste(julia_command,delete_file1_command,delete_file2_command,sep =" && ")
  # }
  # else{
  #   all_command <- paste(julia_command,delete_file1_command,sep =" && ")
  # }
  stopifnot(system(julia_command,wait = T) == 0)
}

data_dir <- paste(root_dir,"anz/result/",cancertype,"/02.downTCGA/",sep = "")
log_dir <- paste(root_dir,"anz/log/",cancertype,"/07.WhippletQuant/",sep = "")
out_dir <- paste(root_dir,"anz/result/",cancertype,"/07.WhippletQuant/",sep = "")
if(!file.exists(log_dir)){
  dir.create(log_dir,recursive = T)
}
if(!file.exists(out_dir)){
  dir.create(out_dir,recursive = T)
}
log_downloaded <- paste(log_dir,"WhippletQuant.sample.log",sep = "")

fastq_dir <- paste(data_dir,manifest$file_id[opt$index],sep = "")
### reads
fastq1 <- list.files(path=fastq_dir, pattern='_1.fastq$',full.names = T)
fastq2 <- list.files(path=fastq_dir, pattern='_2.fastq$',full.names = T)
cat(fastq_dir,"\n")
outprefix <- paste(out_dir,manifest$associated_entities[opt$index],sep = "")

#每一个sample的whippet index
Whipplet_index_jls <- paste(root_dir,"anz/result/",cancertype,"/06.buildWhippletIndex/",manifest$associated_entities[opt$index],".jls",sep = "")
if(length(fastq1) == 1 & length(fastq2) == 1){ ##paired-end
  runWhipplet(root_dir,Whipplet_index_jls,Whipplet_quant_jls,fastq1,fastq2,outprefix)
}else{
  fastq <- list.files(path=fastq_dir, pattern='.fastq$',full.names = T)
  runWhipplet(root_dir,Whipplet_index_jls,Whipplet_quant_jls,fastq,fastq2_dir = "",outprefix)
}

# lizx: do not modify manifest
# cat(manifest$file_id[opt$index],"\n",file = log_downloaded,append =TRUE)

