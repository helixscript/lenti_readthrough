library(dplyr)

o <- readr::read_tsv('transfer.txt', col_names = TRUE)

if(file.exists('log')) file.remove('log')
if(! dir.exists('output')) dir.create('output')

invisible(lapply(split(o, o$sample), function(x){
  write(paste(date(), x$sample[1]), 'log', append = TRUE)
  dir.create(paste0('output/', x$sample[1]))
  dir.create('data')

  message(x$url[1])
  system(paste0('wget -q -O data/', x$fileOut[1], ' ', x$url[1]))

  message(x$url[2])
  system(paste0('wget -q -O data/', x$fileOut[2], ' ', x$url[2]))
  
  comm <- paste0('Rscript ./analyze_sample.R ',
                 ' --fastqStreamChunkSize 10000000',
                 ' --CPUsPerChunk 2',
                 ' --parallelChunks 32',
                 ' --outputDir output/', x$sample[1], 
                 ' --hmmFile /home/ubuntu/projects/Fraietta/ELPS_CD19_BBzeta.hmm',
                 ' --minHMMscore 5',
                 ' --R1 data/', x$fileOut[1], 
                 ' --R2 data/', x$fileOut[2]) 
  system(comm)
  
  unlink('data', recursive = TRUE)
}))
