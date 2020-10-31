


suppressMessages(library(tidyverse))
suppressMessages(library(foreach))
suppressMessages(library(doParallel))
suppressMessages(library(optparse))
suppressMessages(library(R.utils))

option_list <- list(
  make_option(c("-w", "--workDir"), type = "character", default=FALSE,
              help="Root location of all data and codes (workDir/anz)"),
  make_option(c("-c", "--cancer"), type = "character", default=FALSE,
              help="cancer type"),
  make_option(c("-t", "--thread"),type="integer",default=30)
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
print("options:")
print(opt)
stopifnot(!any(sapply(list(opt$workDir,opt$cancer),is.null)),opt$thread>0)
root_dir <- opt$workDir
if(root_dir[length(root_dir)] != "/")
  root_dir <- paste(root_dir,"/",sep = "")
cancertype <- opt$cancer




##STAR 2nd pass re-indexing
STAR_Reindexing <- function(root_dir,cancertype,reference_fasta,gtf){
  STAR_2ndindex <- paste(root_dir,"anz/result/",cancertype,"/04.STAR_2ndIndex",sep = "")
  gtf <- paste(root_dir,"anz/ensembl_v97/Homo_sapiens.GRCh38.97.chr.gtf",sep = "")
  reference_fasta <- paste(root_dir,"anz/ensembl_v97/Homo_sapiens.GRCh38.dna.primary_assembly.fa",sep = "")
  bamdir <- paste(root_dir,"anz/result/",cancertype,"/03.STAR_1stAlign/",sep = "")
  ##re-indexing
    print("Starting STAR reindexing.................\n")
    ##set workdir to where all *SJ.out.tab files locate.
    setwd(bamdir)
    
    ##processing all *SJ.out.tab in this cancertype
    ## referred from https://groups.google.com/forum/#!searchin/rna-star/star$202-pass|sort:date/rna-star/Cpsf-_rLK9I/bGWei7QFBAAJ but a little bit different
    
    ## We only extract unannotated splice sites, in the meantime we consider non-canonical motifs.
    # stopifnot(system(paste("cat ",bamdir,"/*.tab | awk '($6==0)' | cut -f1-6 | sort | uniq > ",
    #                        bamdir,cancertype,".all.SJ.filtered.tab",sep = ""),wait = T) == 0)
    if(!file.exists(STAR_2ndindex)){
      dir.create(STAR_2ndindex,recursive = T)
    }
    command=paste("STAR --runThreadN ",opt$thread,
                           " --runMode genomeGenerate",
                           " --genomeDir ",STAR_2ndindex,
                           " --genomeFastaFiles ",reference_fasta,
                           " --sjdbGTFfile ",gtf,
                           " --sjdbFileChrStartEnd ",bamdir,cancertype,".all.SJ.filtered.tab",
                           " --limitSjdbInsertNsj 50000000",
                           sep = "")
    cat(command,"\n")
    # stop("stop")
    stopifnot(system(command,wait = T) == 0)
  
  print("Finish STAR reindexing!!!!!!\n")
}

STAR_Reindexing(root_dir,cancertype,reference_fasta,gtf)
