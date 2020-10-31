

suppressMessages(library(tidyverse))
suppressMessages(library(foreach))
suppressMessages(library(doParallel))
suppressMessages(library(optparse))
suppressMessages(library(R.utils))

option_list <- list(
  make_option(c("-w", "--workDir"), type = "character",
              help="Root location of all data and codes (workDir/anz)"),
  make_option(c("-c", "--cancer"), type = "character",
              help="cancer name"),
  make_option(c("-i", "--index"),type="integer"),
  make_option(c("-t", "--thread"),type="integer",default=30)
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
print("options:")
print(opt)
stopifnot(!any(sapply(list(opt$workDir,opt$cancer,opt$index),is.null)),opt$thread>0)
root_dir <- opt$workDir
if(root_dir[length(root_dir)] != "/")
  root_dir <- paste(root_dir,"/",sep = "")
cancertype <- opt$cancer
manifest_all <- read.table(paste(root_dir,"anz/data/manifest/TCGA_all_sample_info.txt",sep = ""),
                           sep = "\t",header = T,check.names = F,quote = "",stringsAsFactors = F)
manifest <- manifest_all[manifest_all$cases == cancertype,]

sample_number <- nrow(manifest)
cat(sprintf("number of samples for this cancer type: %d\n",sample_number))
cat(sprintf("current file_id: %s\n",manifest$file_id[opt$index]))
stopifnot(1<=opt$index,opt$index<=sample_number)



##STAR First Align
STAR_FirstAlign <- function(root_dir,cancertype,reads1,reads2 = "",prefix){
  STAR_Firstindex <- paste(root_dir,"anz/STAR_Firstindex",sep = "")
  stopifnot(system(paste("STAR --runThreadN ",opt$thread,
                         " --genomeDir ",STAR_Firstindex,
                         " --readFilesIn ",reads1, " ", reads2,
                         " --sjdbScore 2",
                         " --alignIntronMax 1000000",
                         " --alignMatesGapMax 1000000",
                         " --alignSJDBoverhangMin 1",
                         " --genomeLoad NoSharedMemory",
                         " --sjdbOverhang 100",
                         " --outSAMstrandField intronMotif",
                         " --outSAMtype None",
                         " --outSAMmode None",
                         " --outSJfilterOverhangMin 8 8 8 8",
                         " --outSJfilterCountUniqueMin 1 1 1 1",
                         " --outSJfilterCountTotalMin 1 1 1 1",
                         " --outSJfilterDistToOtherSJmin 0 0 0 0",
                         " --outSJfilterIntronMaxVsReadN 10000000 20000000",
                         " --outFilterMultimapScoreRange 1",
                         " --outFilterMultimapNmax 20",
                         " --outFilterMismatchNmax 10",
                         " --outFilterMatchNminOverLread 0.33",
                         " --outFilterScoreMinOverLread 0.33",
                         " --outFileNamePrefix ",prefix,
                         sep = ""),wait = T) == 0) 
}

data_dir <- paste(root_dir,"anz/result/",cancertype,"/02.downTCGA/",sep = "")
log_dir <- paste(root_dir,"anz/log/",cancertype,"/03.STAR_1stAlign/",sep = "")
out_dir <- paste(root_dir,"anz/result/",cancertype,"/03.STAR_1stAlign/",sep = "")
if(!file.exists(log_dir)){
  dir.create(log_dir,recursive = T)
}
if(!file.exists(out_dir)){
  dir.create(out_dir,recursive = T)
}
log_downloaded <- paste(log_dir,"STAR.1st.sample.log",sep = "")

fastq_dir <- paste(data_dir,manifest$file_id[opt$index],sep = "")
### reads
fastq1 <- list.files(path=fastq_dir, pattern='_1.fastq$',full.names = T)
fastq2 <- list.files(path=fastq_dir, pattern='_2.fastq$',full.names = T)
outprefix <- paste(out_dir,manifest$associated_entities[opt$index],".",sep = "")
if(length(fastq1) == 1 & length(fastq2) == 1){ ##paired-end
  STAR_FirstAlign(root_dir,cancertype,fastq1,fastq2,outprefix)
}else{
  fastq <- list.files(path=fastq_dir, pattern='.fastq$',full.names = T)
  STAR_FirstAlign(root_dir,cancertype,fastq,reads2 ="",outprefix)
}
cat(manifest$file_id[opt$index],"\n",file = log_downloaded,append =TRUE)
  



