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
  make_option(c("-t", "--downThread"), type = "integer", default=15,
              help="download thread")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

root_dir <- opt$workDir
downThread <- opt$downThread
if(root_dir[length(root_dir)] != "/")
  root_dir <- paste(root_dir,"/",sep = "")
cancertype <- opt$cancer


manifest_all <- read.table(paste(root_dir,"anz/data/manifest/TCGA_all_sample_info.txt",sep = ""),
                           sep = "\t",header = T,check.names = F,quote = "",stringsAsFactors = F)
manifest <- manifest_all[manifest_all$cases == cancertype,]

sample_number <- nrow(manifest)
if(downThread > sample_number){
  downThread <- sample_number
}

fixTCGA <- function(root_dir,cancertype,i,manifest){
  
  gdc_client <- paste(root_dir,"anz/code/gdc-client",sep = "")

  token_dir <- paste(root_dir,"anz/data/gdc-user-token.txt",sep = "")
  raw_file_dir <- paste(root_dir,"anz/result/",cancertype,"/02.downTCGA/",sep = "")
  removeDot_dir <- paste(root_dir,"anz/code/removeDot_fastq.py",sep = "")
  log_dir <- paste(root_dir,"anz/log/",cancertype,"/02.downTCGA/",sep = "")

  log_downloaded <- paste(log_dir,"downloaded.sample.log",sep = "")
  log_paired_end <- paste(log_dir,"downloaded.paired.info",sep = "")
  log_single_end <- paste(log_dir,"downloaded.single.info",sep = "")
  log_multi_lane <- paste(log_dir,"downloaded.multilane.info",sep = "")
  
  ###4.转换工作目录到下载完的文件
  now_workDir <- paste(raw_file_dir,manifest$file_id[i],sep = "")
  setwd(now_workDir)
  ###6.判断源文件是什么类型的文件格式
  if(manifest$data_format[i] == "TAR"){ ### STAD,OV,LAML,ESCA
    #### tar -> .gz -> fastq
    fastq1_on <- list.files(path=now_workDir, pattern='_1.fastq.gz$',full.names = T)
    fastq2_on <- list.files(path=now_workDir, pattern='_2.fastq.gz$',full.names = T)
    
    if(length(fastq1_on) == 1 & length(fastq2_on) == 1){ ## paired-end each one fastq
      fastq1_on <- list.files(path=now_workDir, pattern='_1.fastq$',full.names = T)
      fastq2_on <- list.files(path=now_workDir, pattern='_2.fastq$',full.names = T)
      fastq1_rn <- paste(now_workDir,"/",manifest$associated_entities[i],"_1.fastq",sep = "")
      fastq2_rn <- paste(now_workDir,"/",manifest$associated_entities[i],"_2.fastq",sep = "")
      stopifnot(system(paste("cat ",fastq1_on," | python ",removeDot_dir," > ",fastq1_rn,sep = ""),wait = TRUE) == 0)
      stopifnot(system(paste("cat ",fastq2_on," | python ",removeDot_dir," > ",fastq2_rn,sep = ""),wait = TRUE) == 0)
      
      file.remove(fastq1_on)
      file.remove(fastq2_on)
      
      cat(manifest$file_id[i],"\n",file = log_downloaded,append =TRUE)
      cat(manifest$associated_entities[i],"\n",file = log_paired_end,append =TRUE)
    }else if(length(fastq1_on) > 1 & length(fastq2_on) > 1){  ##paired-end each multiply fastq with different lane
      stopifnot(system(paste("cat temp_1.fastq | python ",removeDot_dir," > ",paste(manifest$associated_entities[i],"_1.fastq",sep = ""),sep = ""),wait = TRUE) == 0)
      stopifnot(system(paste("cat temp_2.fastq | python ",removeDot_dir," > ",paste(manifest$associated_entities[i],"_2.fastq",sep = ""),sep = ""),wait = TRUE) == 0)
      ##删除所有的fastq.gz文件以及temp_1.fastq和temp_2.fastq
      stopifnot(system(paste("find ",now_workDir," -name '*.fastq.gz' | xargs rm",sep = ""),wait = TRUE) == 0)
      stopifnot(system(paste("rm -f ",now_workDir,"/temp_1.fastq ",now_workDir,"/temp_2.fastq",sep = ""),wait = TRUE) == 0)
      
      cat(manifest$file_id[i],"\n",file = log_downloaded,append =TRUE)
      cat(manifest$associated_entities[i],"\n",file = log_paired_end,append =TRUE)
      cat(manifest$associated_entities[i],"\n",file = log_multi_lane,append =TRUE)
      
    }else if(length(fastq1_on) == 0 & length(fastq2_on) == 0){ ##single-end
      fastq <- list.files(path=now_workDir, pattern='.fastq.gz$',full.names = T)
      if(length(fastq) == 1){##single-end, one fastq
        fastq_on <- list.files(path=now_workDir, pattern='.fastq$',full.names = T)
        fastq_rn <- paste(now_workDir,"/",manifest$associated_entities[i],".fastq",sep = "")
        stopifnot(system(paste("cat ",fastq_on," | python ",removeDot_dir," > ",fastq_rn,sep = ""),wait = TRUE) == 0)
        ##删除原来的reads
        file.remove(fastq_on)
        cat(manifest$file_id[i],"\n",file = log_downloaded,append =TRUE)
        cat(manifest$associated_entities[i],"\n",file = log_single_end,append =TRUE)

      }else{ ### single-end multiple fastq with different lane
        stopifnot(system(paste("cat temp.fastq | python ",removeDot_dir," > ",paste(manifest$associated_entities[i],".fastq",sep = ""),sep = ""),wait = TRUE) == 0)
        #删除temp.fastq和所有的fastq.gz
        stopifnot(system(paste("find ",now_workDir," -name '*.fastq.gz' | xargs rm",sep = ""),wait = TRUE) == 0)
        stopifnot(system(paste("rm -f ",now_workDir,"/temp.fastq ",sep = ""),wait = TRUE) == 0)
        cat(manifest$file_id[i],"\n",file = log_downloaded,append =TRUE)
        cat(manifest$associated_entities[i],"\n",file = log_single_end,append =TRUE)
        cat(manifest$associated_entities[i],"\n",file = log_multi_lane,append =TRUE)
      }
    }
  }
} 

for(j in 1:sample_number){
  cat("processing index: ",j,"\n")
  fixTCGA(root_dir,cancertype,j,manifest)
}