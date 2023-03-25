#BiocManager::install("ShortRead")
#BiocManager::install("FastqCleaner")
require('Biostrings')
require('ShortRead')
require('FastqCleaner')


#f <- readFastq("sample_JAM16653_S.variegatus_Ilomata.fastq")
#indir <- "~/Downloads/minibar-demux-e2"
#oudir <- "~/Downloads/minibar-demux-e2-filterLength"

length.filter.fastq.directory <- function(indir, outdir = paste0(indir,"/_lengthfiltered"), min, max) {
  setwd(indir)
  dir.create(outdir)
  files <- list.files(pattern=".fastq")
  for (i in 1:length(files)) {
    filename <- files[i]
    f <- readFastq(filename)
    filt <- length_filter(input = f, rm.min = min, rm.max = max)
    outfile <- paste(outdir,filename,sep="/")
    if (length(filt) == 0) {
      print(paste(filename,"had no sequences passing filter"))
    } else {
      writeFastq(filt, file = outfile, compress=F)
    }
  }
}

wd <- getwd()

length.filter.fastq.directory(indir= wd, 
                              outdir=paste0(wd,"/_lengthfiltered"), 
                              min=7000,
                              max=10300)
