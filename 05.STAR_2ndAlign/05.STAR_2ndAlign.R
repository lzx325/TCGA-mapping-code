suppressMessages(library(tidyverse))
suppressMessages(library(foreach))
suppressMessages(library(doParallel))
suppressMessages(library(optparse))
suppressMessages(library(R.utils))

option_list <- list(
  make_option(c("-w", "--workDir"), type = "character", default=FALSE,
              help="Root location of all data and codes (workDir/anz)"),
  make_option(c("-c", "--cancer"), type = "character", default=FALSE,
              help="cancer name"),
  make_option(c("-i", "--index"),type="integer"),
  make_option(c("-t", "--thread"),type="integer",default=30)
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
root_dir <- opt$workDir

print("options:")
print(opt)
stopifnot(!any(sapply(list(opt$workDir,opt$cancer,opt$index),is.null)),opt$thread>0)

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


STAR_SecondAlign <- function(root_dir,cancertype,reads1,reads2 = "",prefix){
  STAR_Secondindex <- paste(root_dir,"anz/result/",cancertype,"/04.STAR_2ndIndex",sep = "")
  stopifnot(system(paste("STAR --runThreadN 10",
                         " --genomeDir ", STAR_Secondindex,
                         " --readFilesIn ",reads1," ", reads2,
                         " --outFilterMultimapScoreRange 1",
                         " --outFilterMultimapNmax 20",
                         " --outFilterMismatchNmax 10", 
                         " --alignIntronMax 1000000",
                         " --alignMatesGapMax 1000000",
                         " --sjdbScore 2", 
                         " --alignSJDBoverhangMin 1", 
                         " --genomeLoad NoSharedMemory", 
                         " --limitBAMsortRAM 70000000000",
                         " --limitOutSJcollapsed 5000000",
                         " --outFilterMatchNminOverLread 0.33",
                         " --outFilterScoreMinOverLread 0.33",
                         " --sjdbOverhang 100",
                         " --outSJfilterOverhangMin 8 8 8 8",
                         " --outSJfilterCountUniqueMin 1 1 1 1",
                         " --outSJfilterCountTotalMin 1 1 1 1",
                         " --outSJfilterDistToOtherSJmin 0 0 0 0",
                         " --outSJfilterIntronMaxVsReadN 10000000 20000000",
                         " --outSAMstrandField intronMotif", 
                         " --outSAMattributes NH HI NM MD AS XS",
                         " --outSAMunmapped Within", 
                         " --outSAMtype BAM SortedByCoordinate",
                         " --outSAMheaderHD @HD VN:1.4",
                         " --outFileNamePrefix ",prefix,
                         sep = ""),wait = T) == 0)
}


Samtools_rmDup <- function(prefix,paired = TRUE){
  STAR_SecondBam <- paste(prefix,"Aligned.sortedByCoord.out.bam",sep = "")
  samtools_rmdup_file <- paste(prefix,"rmdup.bam",sep = "")
  if(paired == TRUE){
    stopifnot(system(paste("samtools rmdup -S ",STAR_SecondBam," ",samtools_rmdup_file,sep = ""),wait = T) == 0)
  }else{
    stopifnot(system(paste("samtools rmdup -s ",STAR_SecondBam," ",samtools_rmdup_file,sep = ""),wait = T) == 0)
  }
  file.remove(STAR_SecondBam)
}

data_dir <- paste(root_dir,"anz/result/",cancertype,"/02.downTCGA/",sep = "")
log_dir <- paste(root_dir,"anz/log/",cancertype,"/05.STAR_2ndAlign/",sep = "")
out_dir <- paste(root_dir,"anz/result/",cancertype,"/05.STAR_2ndAlign/",sep = "")
if(!file.exists(log_dir)){
  dir.create(log_dir,recursive = T)
}
if(!file.exists(out_dir)){
  dir.create(out_dir,recursive = T)
}
log_downloaded <- paste(log_dir,"STAR.2nd.sample.log",sep = "")

fastq_dir <- paste(data_dir,manifest$file_id[opt$index],sep = "")
### reads
fastq1 <- list.files(path=fastq_dir, pattern='_1.fastq$',full.names = T)
fastq2 <- list.files(path=fastq_dir, pattern='_2.fastq$',full.names = T)
outprefix <- paste(out_dir,manifest$associated_entities[opt$index],".",sep = "")
if(length(fastq1) == 1 & length(fastq2) == 1){ ##paired-end
  STAR_SecondAlign(root_dir,cancertype,fastq1,fastq2,outprefix)
  Samtools_rmDup(outprefix)
}else{
  fastq <- list.files(path=fastq_dir, pattern='.fastq$',full.names = T)
  STAR_SecondAlign(root_dir,cancertype,fastq,reads2 ="",outprefix)
  Samtools_rmDup(outprefix,FALSE)
}
cat(manifest$file_id[opt$index],"\n",file = log_downloaded,append =TRUE)

# STAR_SecondAlign <- function(root_dir,cancertype,reads1,reads2 = "",prefix){
#   STAR_Secondindex <- paste(root_dir,"anz/result/",cancertype,"/04.STAR_2ndIndex",sep = "")
#   stopifnot(system(paste("STAR --runThreadN 10",
#                          " --genomeDir ", STAR_Secondindex,
#                          " --readFilesIn ",reads1," ", reads2,
#                          " --outFilterMultimapScoreRange 1",
#                          " --outFilterMultimapNmax 20",
#                          " --outFilterMismatchNmax 10", 
#                          " --alignIntronMax 1000000",
#                          " --alignMatesGapMax 1000000",
#                          " --sjdbScore 2", 
#                          " --alignSJDBoverhangMin 1", 
#                          " --genomeLoad NoSharedMemory", 
#                          " --limitBAMsortRAM 70000000000",
#                          " --outFilterMatchNminOverLread 0.33",
#                          " --outFilterScoreMinOverLread 0.33",
#                          " --sjdbOverhang 100",
#                          " --outSJfilterOverhangMin 8 8 8 8",
#                          " --outSJfilterCountUniqueMin 1 1 1 1",
#                          " --outSJfilterCountTotalMin 1 1 1 1",
#                          " --outSJfilterDistToOtherSJmin 0 0 0 0",
#                          " --outSJfilterIntronMaxVsReadN 10000000 20000000",
#                          " --outSAMstrandField intronMotif", 
#                          " --outSAMattributes NH HI NM MD AS XS",
#                          " --outSAMunmapped Within", 
#                          " --outSAMtype BAM SortedByCoordinate",
#                          " --outSAMheaderHD @HD VN:1.4",
#                          " --outFileNamePrefix ",prefix,
#                          sep = ""),wait = T) == 0)
# }

# data_dir <- paste(root_dir,"anz/result/",cancertype,"/02.downTCGA/",sep = "")
# log_dir <- paste(root_dir,"anz/log/",cancertype,"/05.STAR_2ndAlign/",sep = "")
# out_dir <- paste(root_dir,"anz/result/",cancertype,"/05.STAR_2ndAlign/",sep = "")
# if(!file.exists(log_dir)){
#   dir.create(log_dir,recursive = T)
# }
# if(!file.exists(out_dir)){
#   dir.create(out_dir,recursive = T)
# }
# log_downloaded <- paste(log_dir,"STAR.2nd.sample.log",sep = "")

# fastq_dir <- paste(data_dir,manifest$file_id[opt$index],sep = "")
# ### reads
# fastq1 <- list.files(path=fastq_dir, pattern='_1.fastq$',full.names = T)
# fastq2 <- list.files(path=fastq_dir, pattern='_2.fastq$',full.names = T)
# outprefix <- paste(out_dir,manifest$associated_entities[opt$index],".",sep = "")
# if(length(fastq1) == 1 & length(fastq2) == 1){ ##paired-end
#   STAR_SecondAlign(root_dir,cancertype,fastq1,fastq2,outprefix)
# }else{
#   fastq <- list.files(path=fastq_dir, pattern='.fastq$',full.names = T)
#   STAR_SecondAlign(root_dir,cancertype,fastq,reads2 ="",outprefix)
# }
# cat(manifest$file_id[opt$index],"\n",file = log_downloaded,append =TRUE)
        



