library(rvest)
library(RJSONIO)
library(R.utils)

getOverSamples <- function(multiqc_file,
                           limitPercentPCR = 10) {
  result <- fromJSON(multiqc_file)
  
  samples <- result[["report_plot_data"]][["fastqc_overrepresented_sequencesi_plot"]][["samples"]][[1]]
  percents1 <- result[["report_plot_data"]][["fastqc_overrepresented_sequencesi_plot"]][["datasets"]][[1]][[1]]$data
  percents2 <- result[["report_plot_data"]][["fastqc_overrepresented_sequencesi_plot"]][["datasets"]][[1]][[2]]$data
  percents <- percents1 + percents2
  
  overSamples <- samples[percents > limitPercentPCR]
  
  overSamples
}

getOverrepresentedTable <- function(file) {
  content <- read_html(file)
  if(grepl("No overrepresented sequences",content,fixed=TRUE)) {
    return(NULL)
  }
  table <- (html_table(html_nodes(content,
                                  "table")[2])[[1]])
  
  table
}

getSubsequences <- function(file,
                            overSample,
                            subseqLength = 25) {
  tab <- getOverrepresentedTable(file)
  
  sequences <- tab[startsWith(tab$`Possible Source`,"RNA PCR Primer"),]$Sequence

  subsequences <- c()
  for(sequence in sequences) {
    for(i in 1:(nchar(sequence)-subseqLength+1)) {
      subseq <- substr(sequence,i,i+subseqLength-1)
      subsequences <- c(subsequences,subseq)
    }
  }
  subsequences
}

getBaseName <- function(overSample) {
  overSampleJustName <- strsplit(overSample,"_")[[1]][1]
  overSampleJustName
}

getCounterpart <- function(overSample) {
  overSampleJustName <- getBaseName(overSample)
  lane <- strsplit(overSample,"_")[[1]][2]
  if(lane == "1") {
    otherLane <- "2"
  } else {
    otherLane <- "1"
  }
  overSampleCounterpart <- paste0(overSampleJustName,"_",otherLane)
  overSampleCounterpart
}

processFile <- function(inFileR1,
                         inFileR2,
                         outFileR1,
                         outFileR2,
                        subsequences) {
  connR1 <- gzfile(inFileR1,"r")
  connR2 <- gzfile(inFileR2,"r")
  # open(conn)
  outConnR1 <- gzfile(outFileR1,"w")
  outConnR2 <- gzfile(outFileR2,"w")
  
  linesR1 <- readLines(connR1)
  linesR2 <- readLines(connR2)
  
  n1 <- length(linesR1)/4
  n2 <- length(linesR2)/4
  if(n1 != n2) {
    print("Error: n1 != n2")
  }
  
  flags <- rep(F,rep(n1))
  # filter 1
  for(subseq in subsequences) {
    flags <- flags | grepl(subseq, linesR1[(1:n1)*4-2], fixed = TRUE)
  }
  # filter 2
  flags2 <- rep(T,rep(n1))
  for(i in 1:n1) {
    s <- linesR1[i*4-2]
    aux <- strsplit(s,"")[[1]]
    pA <- sum(aux == "A")/length(aux) * 100
    pC <- sum(aux == "C")/length(aux) * 100
    pT <- sum(aux == "T")/length(aux) * 100
    pG <- sum(aux == "G")/length(aux) * 100
    if(pA > 50 | pC > 50 | pT > 50 | pG > 50) {
      flags2[i] <- F
    }
  }
  
  
  indexRetainedReads <- which(!flags & flags2)
  allFourIndexes <- rep(F,length(linesR1))
  allFourIndexes[indexRetainedReads*4-3] <- T
  allFourIndexes[indexRetainedReads*4-2] <- T
  allFourIndexes[indexRetainedReads*4-1] <- T
  allFourIndexes[indexRetainedReads*4] <- T
  linesR1 <- linesR1[allFourIndexes]
  linesR2 <- linesR2[allFourIndexes]
  
  writeLines(linesR1,outConnR1)
  writeLines(linesR2,outConnR2)
  
  close(connR1)
  close(connR2)
  close(outConnR1)
  close(outConnR2)
}






DIR_MULTIQC <- "/mnt/scratch/sebastian.ciobanu/data-processed/multiqc_data" #"E:/june-july-2021-data/data3/data-processed/multiqc_report/multiqc_data"
DIR_FQC <- "/mnt/scratch/sebastian.ciobanu/data-processed/fqc" #"E:/june-july-2021-data/data3/data-processed/fqc/fqc"
DIR_DATA_ORIGINAL <- "~/data3/data-original/ena_files"
DIR_DATA_PROCESSED_FINAL <- "/mnt/scratch/sebastian.ciobanu/data-processed-mini1/ena_files"

overSamples <- getOverSamples(file.path(DIR_MULTIQC,"multiqc_data.json"))

for(overSample in overSamples) {
  print(overSample)
  subsequences <- getSubsequences(file.path(DIR_FQC,paste0(overSample,"_fastqc.html")),
                                  overSample)
  overSampleJustName <- getBaseName(overSample)
  overSampleCounterpart <- getCounterpart(overSample)
  
  file <- file.path(DIR_DATA_ORIGINAL,overSampleJustName, paste0(overSample,".fastq.gz"))
  fileCounterpart <- file.path(DIR_DATA_ORIGINAL,overSampleJustName, paste0(overSampleCounterpart,".fastq.gz"))
  if(!dir.exists(file.path(DIR_DATA_PROCESSED_FINAL,overSampleJustName))) {
    dir.create(file.path(DIR_DATA_PROCESSED_FINAL,overSampleJustName),
               recursive = T)
  }
  fileOut <- file.path(DIR_DATA_PROCESSED_FINAL,overSampleJustName, paste0(overSample,".fastq.gz"))
  fileCounterpartOut <- file.path(DIR_DATA_PROCESSED_FINAL,overSampleJustName, paste0(overSampleCounterpart,".fastq.gz"))
  # citit + filtrat AMBELE cu AMBELE criterii; de salvat in DIR_DATA_PROCESSED_FINAL
  
  processFile(inFileR1=file,
              inFileR2=fileCounterpart,
              outFileR1=fileOut,
              outFileR2=fileCounterpartOut,
              subsequences=subsequences)
}

