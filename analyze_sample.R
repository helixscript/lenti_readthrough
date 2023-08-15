#!/usr/bin/Rscript
library(dplyr)
library(optparse)
library(Biostrings)
library(ShortRead)
library(parallel)

tmpFile <- function(){ paste0(paste0(stringi::stri_rand_strings(30, 1, '[A-Za-z0-9]'), collapse = '')) }

option_list = list(
  make_option(c("--R1"), type="character", default=NULL, help="R1 fasta file", metavar="character"),
  make_option(c("--R2"), type="character", default=NULL, help="R2 fasta file", metavar="character"),
  make_option(c("--fastqStreamChunkSize"), type="integer", default=5000000, help="fastqStreamChunkSize", metavar="integer"),
  make_option(c("--CPUsPerChunk"), type="integer", default=2, help="CPUs per parallel chunk", metavar="integer"),
  make_option(c("--parallelChunks"), type="integer", default=15, help="parallel chunks", metavar="integer"),
  make_option(c("--outputDir"), type="character", default=NULL, help="output directory", metavar="character"),
  make_option(c("--hmmFile"), type="character", default=NULL, help="hmm file path", metavar="character"),
  make_option(c("--minHMMscore"), type="numeric", default=20, help="min score to accept an hmm hit", metavar="numeric")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

opt$fastqStreamChunkSize <- 10000000
opt$CPUsPerChunk <- 2 
opt$parallelChunks <- 15 
opt$outputDir <- 'output/16CT023-04_D10'
opt$hmmFile <- '/home/ubuntu/projects/Fraietta/ELPS_CD19_BBzeta.hmm'
opt$minHMMscore <- 5
opt$R1 <- 'data/R1.fastq.gz' 
opt$R2 <- 'data/R2.fastq.gz'


if(! dir.exists(opt$outputDir )) dir.create(opt$outputDir)

R1_strm <- FastqStreamer(opt$R1, n = opt$fastqStreamChunkSize)
R2_strm <- FastqStreamer(opt$R2, n = opt$fastqStreamChunkSize)

nReads <- 0
logFile <- file.path(opt$outputDir, 'log')

write(date(), 'log', append = FALSE)

repeat {
  write(paste0('  ', date(), ' Starting new sequence chunk.'), logFile, append = TRUE)
  
  write(paste0('  ', date(), '   Pulling R1 chunk'), logFile, append = TRUE)
  R1_fq <- yield(R1_strm)
  write(paste0('  ', date(), '   Pulling R2 chunk'), logFile, append = TRUE)
  R2_fq <- yield(R2_strm)
  
  nReads <- nReads + length(R1_fq)
  
  if(length(R1_fq) == 0 | length(R2_fq) == 0) break
  
  i <- dplyr::ntile(1:length(R1_fq), opt$parallelChunks)
  R1 <- split(R1_fq, i)
  R2 <- split(R2_fq, i)
  
  rm(i, R1_fq, R2_fq)
  
  write('    Running HMM.', logFile, append = TRUE)
  o <- mapply(function(R1, R2){
    
    tmp <- tmpFile()
    dir.create(file.path(opt$outputDir, tmp))
    
    writeFasta(R1, file.path(opt$outputDir, tmp, 'R1.fasta'))
    writeFasta(R2, file.path(opt$outputDir, tmp, 'R2.fasta'))
    
    write(c('#!/usr/bin/sh', 
            paste0('./runLTRhmm.R --R1 ', file.path(opt$outputDir, tmp, 'R1.fasta'), 
                   ' --R2 ', file.path(opt$outputDir, tmp, 'R2.fasta'), ' --outputDir ',
                   file.path(opt$outputDir, tmp), ' --CPUs ', opt$CPUsPerChunk, ' --hmmFile ', opt$hmmFile, ' --minHMMscore ', opt$minHMMscore)), 
          file = paste0(tmp, '.sh'))
    
    system(paste0('chmod 755 ./', tmp, '.sh'))
    system(paste0('./', tmp, '.sh'), wait = FALSE)
    file.path(opt$outputDir, tmp, 'done')
  }, R1, R2, SIMPLIFY = TRUE)
  
  write('    HMMs jobs submitted.', logFile, append = TRUE)
  
  write('    Waiting for jobs to complete.', logFile, append = TRUE)
  for(d in o){
    while(! file.exists(d)){
      Sys.sleep(1)
    }
    
    invisible(file.remove(paste0(rev(unlist(stringr::str_split(d, '/')))[2], '.sh')))
  }
  
  r <- bind_rows(lapply(list.files(opt$outputDir, pattern = '^result$', recursive = TRUE, full.names = TRUE), readr::read_tsv, col_names = TRUE))
  
  if(nrow(r) > 0){
    tmp <- tmpFile()
    readr::write_tsv(r, file = file.path(opt$outputDir, paste0('result_', tmp, '.tsv')))

    R1 <- Reduce('append', R1)
    R2 <- Reduce('append', R2)
     
    R1_ids <- sub('\\s.+', '', R1@id)
    R2_ids <- sub('\\s.+', '', R2@id)
     
    writeFastq(R1[R1_ids %in% r$targetName], file.path(opt$outputDir, paste0('R1_', tmp, 'fastq.gz')))
    writeFastq(R2[R2_ids %in% r$targetName], file.path(opt$outputDir, paste0('R2_', tmp, 'fastq.gz')))
  }
  
  invisible(lapply(list.dirs(opt$outputDir, recursive = FALSE), unlink, recursive = TRUE))
}

write(paste0(date(), ' - Done. ', nReads, ' reads processed'), logFile, append = TRUE)
