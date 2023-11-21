# bedtools get fasta from bed

#usage
#getFastaFromBed(genome = genome, bed = "hypoxia-exons-padded-mergeClose.bed", output = "hypoxia-TaeGut.fasta", strandedness = T)


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

