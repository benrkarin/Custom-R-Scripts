# bedtools get fasta from bed

#usage
#genomepath="PATH/TO/GENOME.fas"
#setwd("/PATH/TO/OUTFOLDER")
#getFastaFromBed(genome = genomepath, bed = "hypoxia-exons-padded-mergeClose.bed")

#Output will be same name as .bed file but extension changed to .fasta

# REQUIRES bedtools installed in command line path

getFastaFromBed <- function(genome,bed, strandedness = T, bedName=T) {
  filenamenoext <- strsplit(bed,".bed")[[1]]
  output = paste0(filenamenoext,".fasta")
  if (strandedness == T) {
    c <- paste0("bedtools getfasta -fi ", genome," -bed ", bed," -fo ", output, " -nameOnly -s")  
  } else if (bedName == T) {
    c <- paste0("bedtools getfasta -fi ", genome," -bed ", bed," -fo ", output, " -nameOnly")
  } else if (bedName == F) {
    c <- paste0("bedtools getfasta -fi ", genome," -bed ", bed," -fo ", output)
  }
  system(c)
}

