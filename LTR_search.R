library(dplyr)
library(ShortRead)
library(parallel)

dataDir   <- '/home/ubuntu/projects/Fraietta/data'
hmmerBin  <- '/home/ubuntu/software/hmmer-3.3.2/bin'
outputDir <- '/home/ubuntu/projects/Fraietta/output'

CPUs <- 38
readStrmLength <- 1000000
minHMMscore <- 15
minReadLengthPostTrim <- 75

f <- list.files(dataDir, recursive = TRUE, full.names = TRUE, pattern = 'fastq', no.. = TRUE)
f <- f[grepl('_R1_', f)]

write(date(), file.path(outputDir, 'log'), append = FALSE)

result <- tibble()

runHMM <- function(fq){
  library(ShortRead)
  library(dplyr)
  
  tmp <- tempfile(pattern = label, tmpdir = file.path(outputDir, 'tmp'))
    
  writeFasta(fq, paste0(tmp, '.fasta'))
    
  comm <- paste0(file.path(hmmerBin, 'hmmsearch'), ' --tblout ',  paste0(tmp, '.tbl'), 
                 ' --domtblout ',  paste0(tmp, '.domTbl'), ' HIV_U5_1-100.hmm ', ' ', paste0(tmp, '.fasta'), ' > ', paste0(tmp, '.hmmSearch'))
  system(comm)
    
  r <- readLines(paste0(tmp, '.domTbl'))
  r <- r[!grepl('^\\s*#', r)]
  r <- strsplit(r, '\\s+')
    
  o <- bind_rows(lapply(r, function(x) data.frame(t(x))))
    
  if(nrow(o) > 0){
      names(o) <- c('targetName', 'targetAcc', 'tlen', 'queryName', 'queryAcc', 'queryLength', 'fullEval', 
                    'fullScore', 'fullBias', 'domNum', 'totalDoms', 'dom_c-Eval', 'dom_i-Eval', 'domScore', 
                    'domBias', 'hmmStart', 'hmmEnd', 'targetStart', 'targetEnd', 'envStart', 'envEnd', 
                    'meanPostProb',  'desc') 
      
      write.table(o, sep = '\t', file = paste0(tmp, '.domTbl2'), col.names = TRUE, row.names = FALSE, quote = FALSE)
      
      o <- readr::read_delim(paste0(tmp, '.domTbl2'), '\t', col_types = readr::cols())
      o$fullScore <- as.numeric(o$fullScore)
      o$fullEval  <- as.numeric(o$fullEval)
      
      o <- subset(o, fullScore >= minHMMscore)
      
      if(nrow(o) > 0) readr::write_tsv(o, paste0(tmp, '.result'))
    }
}


processChunk <- function(fq0, CPUs, label){
  cluster <- makeCluster(CPUs)
  clusterExport(cluster, envir = environment(), c('dataDir', 'hmmerBin', 'outputDir', 'readStrmLength', 'minHMMscore', 'label'))
  
  if(label == '_R2_') fq0 <- reverseComplement(fq0)

  o <- split(fq0, dplyr::ntile(1:length(fq0), CPUs))
  rm(fq0)
  
  write(paste0('    Chunk reads split into ', length(o), ' pieces.'), file.path(outputDir, 'log'), append = TRUE)
  invisible(parLapply(cluster, o, runHMM))
  
  stopCluster(cluster)
}


invisible(lapply(f, function(x){
  write(paste('Starting', x), file.path(outputDir, 'log'), append = TRUE)
  
  R1_strm <- FastqStreamer(x, n = readStrmLength)
  R2_strm <- FastqStreamer(sub('_R1_', '_R2_', x), n = readStrmLength)
  
  repeat {
    R1_fq <- yield(R1_strm)
    R2_fq <- yield(R2_strm)
    write(paste0('  ', date(), ' Starting new sequence chunk.'), file.path(outputDir, 'log'), append = TRUE)
    write('    Triming and preparing reads.', file.path(outputDir, 'log'), append = TRUE)
    R1_fq <- trimTailw(R1_fq, 2, "+", 5)
    R2_fq <- trimTailw(R2_fq, 2, "+", 5)
    
    R1_fq <- R1_fq[width(R1_fq) >= minReadLengthPostTrim]
    R2_fq <- R2_fq[width(R2_fq) >= minReadLengthPostTrim]
    
    R1.ids <- sub('\\s.+$', '', as.character(R1_fq@id))
    R2.ids <- sub('\\s.+$', '', as.character(R2_fq@id))
    
    write('    Syncing reads.', file.path(outputDir, 'log'), append = TRUE)
    n <- base::intersect(R1.ids, R2.ids)
    
    R1_fq <- R1_fq[R1.ids %in% n]
    R2_fq <- R2_fq[R2.ids %in% n]
    
    i <- duplicated(xscat(R1_fq@sread, R2_fq@sread))
    
    R1_fq <- R1_fq[!i]
    R2_fq <- R2_fq[!i]
    
    if(length(R1_fq) == 0 | length(R2_fq) == 0) break
    
    write('    Starting R1 analysis.', file.path(outputDir, 'log'), append = TRUE)
    processChunk(R1_fq, CPUs, '_R1_')
    
    write('    Starting R2 analysis.', file.path(outputDir, 'log'), append = TRUE)
    processChunk(R2_fq, CPUs, '_R2_')
    
    gc()
    
    write('    Gathering results.', file.path(outputDir, 'log'), append = TRUE)
    rf <- list.files(file.path(outputDir, 'tmp'), pattern = 'result', full.names = TRUE)
    R1_tab <- bind_rows(lapply(rf[grepl('_R1_', rf)], readr::read_tsv))
    R2_tab <- bind_rows(lapply(rf[grepl('_R2_', rf)], readr::read_tsv))
    
    if(nrow(R1_tab) > 0 | nrow(R2_tab)){
      R1_ids <- sub('\\s.+$', '', as.character(R1_fq@id)) # R2 ids will be the same.
      indices <- stringr::str_extract(as.character(R1_fq@id), '[ATCGN]+$')
      
      if(nrow(R1_tab) > 0 ){
        R1_tab$LTRhit <- 'R1'
        i <- match(R1_tab$targetName, R1_ids)
        R1_tab$R1_seq <- as.character(R1_fq[i]@sread)
        R1_tab$R2_seq <- as.character(R2_fq[i]@sread)
        R1_tab$index_seq <- indices[i]
      }
      
      if(nrow(R2_tab) > 0 ){
        R2_tab$LTRhit <- 'R2'
        i <- match(R2_tab$targetName, R1_ids)
        R2_tab$R1_seq <- as.character(R1_fq[i]@sread)
        R2_tab$R2_seq <- as.character(R2_fq[i]@sread)
        R2_tab$index_seq <- indices[i]
      }
      
      tab <- bind_rows(R1_tab, R2_tab) %>%
        dplyr::select(targetName, LTRhit, fullEval, fullScore, hmmStart, hmmEnd, targetStart, targetEnd, envStart, envEnd, R1_seq, R2_seq, index_seq) %>%
        dplyr::distinct()
      
      write(paste0('    ', nrow(tab), ' significant hits found.'), file.path(outputDir, 'log'), append = TRUE)
      
      tab$file <- x
      
      result <<- bind_rows(result, tab)
      readr::write_tsv(result, file.path(outputDir, 'readsWithLTRsig.tsv'))
    } else {
      write('    No significant hits found.', file.path(outputDir, 'log'), append = TRUE)
    }
    
    invisible(file.remove(list.files(file.path(outputDir, 'tmp'), full.names = TRUE)))
  }
}))

write(paste0('  ', date(), ' Done.'), file.path(outputDir, 'log'), append = TRUE)

o <- readr::read_tsv(file.path(outputDir, 'readsWithLTRsig.tsv'))
z <- subset(o, LTRhit == 'R1' & fullScore >= 30 & hmmEnd >= 100)

invisible(lapply(split(z, 1:nrow(z)), function(x){
  write(paste0(x$targetName,'\n>R1\n', x$R1_seq, '\n>R2\n', x$R2_seq, '\n\n'), 
        file = file.path(outputDir, 'R1_candidates.txt'), append = TRUE)
}))



