library(dplyr)
library(ShortRead)
source('lib.R')

refGenome        <- '/home/ubuntu/AAVengeR/data/blatDBs/hg38.2bit'
genesReference   <- '/home/ubuntu/AAVengeR/data/genomeAnnotations/hg38.TUs.rds'
exonReference    <- '/home/ubuntu/AAVengeR/data/genomeAnnotations/hg38.exons.rds'

r <- bind_rows(lapply(list.files(pattern = 'HMM_result.tsv', recursive = TRUE, full.names = TRUE), function(x){
       o <- readr::read_tsv(x)
       o$sample <- unlist(strsplit(x, '/'))[2]
       select(o, sample, targetName, tlen, fullEval, fullScore, targetStart, targetEnd, hmmStart, hmmEnd, R1_seq, R2_seq) %>% distinct()
      }))

C1 <- Reduce('append', lapply(list.files(pattern = 'C1.fastq.gz', recursive = TRUE, full.names = TRUE), function(x){
        readFastq(x)
      }))

U1 <- Reduce('append', lapply(list.files(pattern = 'U1.fastq.gz', recursive = TRUE, full.names = TRUE), function(x){
  readFastq(x)
}))

r <- left_join(r, data.frame(targetName = sub('\\s.+$', '', as.character(C1@id)),
                             cellBarCode = as.character(C1@sread)),
               by = 'targetName')

r <- left_join(r, data.frame(targetName = sub('\\s.+$', '', as.character(U1@id)),
                             transcriptUMI = as.character(U1@sread)),
               by = 'targetName')

r$sample <- sub('^Output\\-', '', r$sample)
readr::write_tsv(r, 'full_result.tsv')
openxlsx::write.xlsx(r, 'full_result.xlsx')


r <- filter(r, hmmEnd == 100)

r <- r[vcountPattern('TTTCTTATATGGG', subseq(DNAStringSet(r$R1_seq), 25, 55), max.mismatch = 2) == 1,]



# Create blat input.
if(file.exists('seqs.ff')) file.remove('seqs.ff')

s <- unlist(lapply(split(r, 1:nrow(r)), function(x){
  paste0('>',  x$sample, '___', x$targetName, '___R1\n', substr(x$R1_seq, 40, nchar(x$R1_seq)), 
         '\n>', x$sample, '___', x$targetName, '___R2\n', x$R2_seq)
}))

write(s, file = 'seqs.ff')

# Align candidate read pairs.
system(paste('blat -stepSize=5 -repMatch=5000 -minScore=0 -minIdentity=0 -out=psl -noHead', 
             refGenome, 'seqs.ff', 'seqs.psl'))

alignmentMinSeqId <- 97

b <- parseBLAToutput('seqs.psl') %>% 
     tidyr::separate(qName, into = c('sample', 'readID', 'direction'), sep = '___') %>%
     dplyr::filter(alignmentPercentID >= alignmentMinSeqId & 
                   tNumInsert <= 1 & qNumInsert <= 1 & tBaseInsert <= 2 & qBaseInsert <= 2)

f <- bind_rows(lapply(split(b, b$readID), function(x){
  R1 <- subset(x, direction == 'R1')
  R2 <- subset(x, direction == 'R2')
  if(nrow(R1) == 0 | nrow(R2) == 0) return(tibble())
  
  R2 <- subset(R2, qStart <= 5)
  
  R1 <- dplyr::arrange(R1, desc(matches)) %>% dplyr::select(tName, strand, tStart, tEnd, qStart, qEnd, matches)
  R2 <- dplyr::arrange(R2, desc(matches)) %>% dplyr::select(tName, strand, tStart, tEnd, qStart, qEnd, matches)
  
  names(R1) <- paste0(names(R1), '_R1')
  names(R2) <- paste0(names(R2), '_R2')
  
  j <- full_join(R1, R2, by = c('tName_R1' = 'tName_R2'))
  j <- dplyr::filter(j, abs(j$tStart_R1 - j$tStart_R2) <= 1000 & strand_R1 != strand_R2)
  
  if(nrow(j) > 0){
    return(rowwise(j) %>% 
             mutate(m = sum(matches_R1, matches_R2)) %>% 
             ungroup() %>% 
             dplyr::slice_max(m, with_ties = FALSE) %>%
             dplyr::select(-m) %>%
             mutate(readID = x$readID[1]))
  } else {
    return(tibble())
  }
})) %>% left_join(dplyr::select(r, sample, tlen, targetName, R1_seq, R2_seq), by = c('readID' = 'targetName'))

readr::write_tsv(f, 'fragments.tsv')

# Create a BED file for the UCSC browser.
# Different locations can be accessed by providing a CGI argument, eg. &position=chrX%3A133373794-133374794
# Positive strand: blue; negative strand: red.

bedData <- c(paste0('track name="LTR_readthrough2" description="LTR_readthrough2" visibility=2 itemRgb="On"'),
             unlist(lapply(split(f, 1:nrow(f)), function(x){
               posid <- paste0(x$tName_R1, '_', x$tEnd_R1)
               
               if(x$strand_R1 == '+'){
                 R1_startPos      <- x$tStart_R1 - x$qStart_R1  
                 R1_thickStartPos <- x$tStart_R1
                 R1_thickEndPos   <- x$tEnd_R1
                 R1_endPos        <- x$tEnd_R1 + ((x$tlen-39)-(x$qEnd_R1-0))
                 R1_color         <- '0,0,255'
               } else {
                 R1_startPos      <- x$tStart_R1 - ((x$tlen-39) - (x$qEnd_R1-0)) 
                 R1_thickStartPos <- x$tStart_R1
                 R1_thickEndPos   <- x$tEnd_R1
                 R1_endPos        <- x$tEnd_R1 + x$qStart_R1
                 R1_color         <- '255,0,0'
               }
               
               if(x$strand_R2 == '+'){
                 R2_startPos      <- x$tStart_R2 - x$qStart_R2 
                 R2_thickStartPos <- x$tStart_R2
                 R2_thickEndPos   <- x$tEnd_R2
                 R2_endPos        <- x$tEnd_R2 + (x$tlen-x$qEnd_R2)
                 R2_color         <- '0,0,255'
               } else {
                 R2_startPos      <- x$tStart_R2 - (x$tlen - x$qEnd_R2)
                 R2_thickStartPos <- x$tStart_R2
                 R2_thickEndPos   <- x$tEnd_R2
                 R2_endPos        <- x$tEnd_R2 + x$qStart_R2
                 R2_color         <- '255,0,0'
               }
               
               c(paste0(x$tName_R1, '\t', R1_startPos, '\t', R1_endPos, '\t', 'LTR_read\t0\t', x$strand_R1, '\t', R1_thickStartPos, '\t', R1_thickEndPos, '\t', R1_color),
                 paste0(x$tName_R1, '\t', R2_startPos, '\t', R2_endPos, '\t', 'Break_read\t0\t', x$strand_R2, '\t', R2_thickStartPos, '\t', R2_thickEndPos, '\t', R2_color))
             })))

write(bedData, 'LTR_readthrough2.ucsc')

f$posid <- paste0(f$tName_R1, f$strand_R1, f$tStart_R1)

library(GenomicRanges)
n <- nearestGene(f$posid, readRDS(genesReference), readRDS(exonReference))
n$posid <- paste0(n$chromosome, n$strand, n$position)

f <- left_join(f, distinct(n), by = 'posid')

f$UCSC_link <- paste0('HYPERLINK("http://genome.ucsc.edu/cgi-bin/hgTracks?org=human&db=hg38&hgt.customText=https://microb120.med.upenn.edu/data/export/everett/ucsc/LTR_readthrough2.ucsc&',
                      'position=', f$tName_R1, '%3A', f$tStart_R1-500, '-', f$tStart_R1+500,'", "UCSC link")')

class(f$UCSC_link) <- "formula"      

out <- select(f, sample, UCSC_link, posid, nearestGene, nearestGeneDist, readID, R1_seq, R2_seq)
out$sample <- sub('Output\\-', '', out$sample)

out2 <- left_join(out, distinct(select(r, targetName, targetEnd)), by = c('readID' = 'targetName'))
out2 <- data.frame(out2)
class(out2$UCSC_link) <- "formula" 

out2 <- left_join(out2, distinct(select(r, targetName, cellBarCode, transcriptUMI)),
                  by = c('readID' = 'targetName'))

out2 <- relocate(out2, cellBarCode, transcriptUMI, .after = 'nearestGeneDist') %>% data.frame()
class(out2$UCSC_link) <- "formula"

library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "result")
writeData(wb, "result", out2)
saveWorkbook(wb, "result.xlsx", overwrite = TRUE)
