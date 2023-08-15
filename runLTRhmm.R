#!/usr/bin/Rscript
library(dplyr)
library(optparse)
library(Biostrings)
library(parallel)

option_list = list(
  make_option(c("--R1"), type="character", default=NULL, help="R1 fasta file", metavar="character"),
  make_option(c("--R2"), type="character", default=NULL, help="R2 fasta file", metavar="character"),
  make_option(c("--CPUs"), type="integer", default=10, help="R2 fasta file", metavar="integer"),
  make_option(c("--outputDir"), type="character", default=NULL, help="output directory", metavar="character"),
  make_option(c("--hmmFile"), type="character", default=NULL, help="hmm file path", metavar="character"),
  make_option(c("--minHMMscore"), type="numeric", default=20, help="min score to accept an hmm hit", metavar="numeric")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

tmpFile <- function(){ paste0(paste0(stringi::stri_rand_strings(30, 1, '[A-Za-z0-9]'), collapse = ''), '.tmp') }

# opt$R1 <- '/home/ubuntu/projects/Fraietta/Output-19CT023-04-D10/R2.fasta'
# opt$R2 <- '/home/ubuntu/projects/Fraietta/Output-19CT023-04-D10/R1.fasta'
# opt$outputDir <- '/home/ubuntu/projects/Fraietta/Output-16CT023-04_D10/gxfTaDTpfTsFj8941ChPtDudoCn7zM'
# opt$CPUs <- 3
# opt$label <- '_R1_'
# opt$hmmFile <- '/home/ubuntu/projects/Fraietta/ELPS_CD19_BBzeta.hmm'
# opt$minHMMscore <- 1

cluster <- makeCluster(opt$CPUs)
clusterExport(cluster, c('opt', 'tmpFile'))

R1 <- readDNAStringSet(opt$R1)
R2 <- readDNAStringSet(opt$R2)

i <- ntile(1:length(R1), opt$CPUs)

R1s <- split(R1, i)
R2s <- split(R2, i)

rm(R1, R2)

o <- lapply(1:length(R1s), function(x){
       list(R1 = R1s[[x]], R2 = R2s[[x]])
      })

rm(R1s, R2s)

worker <- function(x){
  library(dplyr)
  library(Biostrings)
  tmp <- tmpFile()
  dir.create(file.path(opt$outputDir, tmp))
  writeXStringSet(x$R1, file.path(opt$outputDir, tmp, 'R1.fasta'))
  writeXStringSet(x$R2, file.path(opt$outputDir, tmp, 'R2.fasta'))

  system(paste0('hmmsearch --tblout ', file.path(opt$outputDir, tmp, 'tbl'), 
                ' --domtblout ', file.path(opt$outputDir, tmp, 'domTbl'), ' ', 
                opt$hmm, ' ', file.path(opt$outputDir, tmp, 'R1.fasta'), 
                ' > ', file.path(opt$outputDir, tmp, 'hmmSearch')))


  r <- readLines(file.path(opt$outputDir, tmp, 'domTbl'))
  r <- r[!grepl('^\\s*#', r)]
  r <- strsplit(r, '\\s+')

  o <- bind_rows(lapply(r, function(x) data.frame(t(x))))

  if(nrow(o) > 0){
    names(o) <- c('targetName', 'targetAcc', 'tlen', 'queryName', 'queryAcc', 'queryLength', 'fullEval', 
                  'fullScore', 'fullBias', 'domNum', 'totalDoms', 'dom_c-Eval', 'dom_i-Eval', 'domScore', 
                  'domBias', 'hmmStart', 'hmmEnd', 'targetStart', 'targetEnd', 'envStart', 'envEnd', 
                  'meanPostProb',  'desc') 
  
    write.table(o, sep = '\t', file = file.path(opt$outputDir, tmp, 'domTbl2'), col.names = TRUE, row.names = FALSE, quote = FALSE)
  
    o <- readr::read_delim(file.path(opt$outputDir, tmp, 'domTbl2'), '\t', col_types = readr::cols())
    o$fullScore <- as.numeric(o$fullScore)
    o$fullEval  <- as.numeric(o$fullEval)
  
    o <- subset(o, fullScore >= as.numeric(opt$minHMMscore))
    
    if(nrow(o) > 0){
      ids <- sub('\\s.+$', '', names(x$R1))
      o$R1_seq <- as.character(x$R1[match(o$targetName, ids)])
      o$R2_seq <- as.character(x$R2[match(o$targetName, ids)])
      
      readr::write_tsv(o, file.path(opt$outputDir, tmp, 'result'))
    }
  }
}

invisible(parLapply(cluster, o, worker))
#invisible(lapply(o, worker))

stopCluster(cluster)

r <- bind_rows(lapply(list.files(file.path(opt$outputDir), pattern = 'result', full.names = TRUE, recursive = TRUE), readr::read_tsv))
readr::write_tsv(r, file.path(opt$outputDir, 'result'))

Sys.sleep(1)
write(date(), file.path(opt$outputDir, 'done'))
q(save = 'no', status = 0)
