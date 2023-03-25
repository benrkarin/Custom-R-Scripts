
## This function trims all fasta alignments in a folder to whatever gaps 
## are present at the beginning/end of sequence 1.

#### USAGE ####
# indir <- c("~/Dropbox/bird_phylogenomics_data/newBirdGenomeAdditions/NEW-GENOMES/Tricholaema_leucomelas/aligned-fastas")
# trimFlankToRef(directory = indir)

require(ape)

trimFlankToRef <- function(directory = getwd()) {
  setwd(directory)
  system(paste("mkdir ", directory,"/_alignments_trimmed",sep="")) # make the fasta extraction directory
  out <- paste(directory,"/_alignments_trimmed",sep="")
  files <- list.files(directory, pattern =".fas")
  for (i in 1:length(files)) { #the number of genes
    locus <- files[i]
    dna.matrix <- read.dna(file=paste(directory,"/",locus,sep=""), format="fasta", as.matrix=T)
    aln.i <- as.alignment(dna.matrix)
    aln.taxa <- rownames(dna.matrix)
    current.taxon <- aln.taxa[i]
    refseq.dnabin <- dna.matrix[1,]
    refseq.aln <- as.character(refseq.dnabin)[1,]
    refseq.rle <- rle(refseq.aln)
    if (refseq.rle$values[1]=="-") { frontclip <- refseq.rle$lengths[1] + 1 }   # set the amount to clip at the front
    else { frontclip <- 1 }
    if (refseq.rle$values[length(refseq.rle$values)]=="-") { backclip <- ncol(dna.matrix) - refseq.rle$lengths[length(refseq.rle$values)] }   # set end clip
    else { backclip <- length(refseq.rle$values) }
    trimmed.aln <- dna.matrix[,frontclip:backclip]
    write.dna(trimmed.aln, file=paste0(out,"/",locus), format="fasta", colsep="",nbcol=-1)
  }
  
}

# Run in terminal directly to get rid of reference and "__"
# remove reference sequences and format names (get rid of anything before __taxonname).
# Set working directory first!
# for filename in *.fas; do
# find "$filename" -type f -exec sed -i '' '/__reference/ { N; d; }' {} \;
# find "$filename" -type f -exec sed -i '' "s/>.*__/>/g" {} \;
# done