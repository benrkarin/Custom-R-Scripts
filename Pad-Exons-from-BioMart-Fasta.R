####################################################################################################
# Load in Functions Below
####################################################################################################
# FastaToTabular was downlaoded from "https://raw.githubusercontent.com/lrjoshi/FastaTabular/master/fasta_and_tabular.R"
# THANKS!

## dat should be a data frame with the appropriate column headers. Can create this from FastaToTabular

#FastaToTabular("NuMt-mart_export.txt")

## clean up output
#system("sed -i \'\' 1d dna_table.csv")
#system("sed -i \'\' \'s/|/\",\"/g\' dna_table.csv")
#system("sed -i \'\' \'s/>//g\' dna_table.csv")

## Rename the column names, 
## download from Biomart in this order and then confirm correct order after swap like so
#colnames(dat) <- c(
#  "order",
#  "Gene.stable.ID",
#  "Transcript.stable.ID.version",
#  "Exon.stable.ID",
#  "Exon.rank.in.transcript",
#  "Strand",
#  "Gene.name",
#  "Exon.region.start..bp.",
#  "Exon.region.end..bp.",
#  "Genomic.coding.start",
#  "Genomic.coding.end",
#  "sequence")

#head(dat[-ncol(dat)],n=5L)

## Sort by Gene/Exon then trim out noncoding Exons
#dat <- dat[order(dat$Gene.name, dat$Transcript.stable.ID.version, dat$Exon.rank.in.transcript),]
#dat <- dat %>% filter(!is.na(Genomic.coding.end))
#
## confirm the apropriate columns are numerics
#dat$Exon.region.start..bp. <- as.numeric(dat$Exon.region.start..bp.)
#dat$Exon.region.end..bp. <- as.numeric(dat$Exon.region.end..bp.)
#dat$Genomic.coding.start <- as.numeric(dat$Genomic.coding.start)
#dat$Genomic.coding.end <- as.numeric(dat$Genomic.coding.end)
#
#write.csv(dat,file="dna_table_clean.csv")

## martFlank is the number of bases collected on each side of the sequence. 
## 100 on front and 100 on back: martFlank = 100
## Default 100

## minLength is the amount to extend the short exons to. 
## Default 200

## Output is just the file name and it will go in the working directory. 
# If the file exists, it will confirm if you want to append or delete the file.

# options for name are "GeneNameThenTranscriptID" (default), "GeneNameThenGeneID", "TranscriptID", or "GeneID"
# the first two give the secondary option if no GeneName is available.



########################################################################################################################
require(ape)
require(devtools)
require(tidyverse)

trimFlankingSeqstoProperPadding <- function(dat,martFlank = 100, minLength = 200, OUTPUT, name="GeneNameThenTranscriptID", trimUTR = T) {
  
  # check if output file exists
  if (file.exists(OUTPUT)) {
    deleteQ <- readline("Output fasta file already exists. Would you like to remove before continuing? 
    (If no, output file will be appended) (y/n) ")
    if(deleteQ=="y") file.remove(OUTPUT)
  }
  
  for (i in 1:nrow(dat)) {
    # find sequence lengths, etc.
    exstart <- dat$Exon.region.start..bp.[i]
    exend <- dat$Exon.region.end..bp.[i]
    cdsstart <- dat$Genomic.coding.start[i]
    cdsend <- dat$Genomic.coding.end[i]
    
    exonLength <- exend - exstart + 1
    CDSLength <- cdsend - cdsstart + 1
    frontUTR <- cdsstart - exstart
    backUTR <- exend - cdsend
    
    dnaString <- dat$sequence[i] 
    seqLength <- nchar(dnaString)
    
    #exStartPos <- martFlank + 1
    #exEndPos <- seqLength - martFlank
    #exLength <- exEndPos - exStartPos +1
    #cdsStartPos <- exStartPos + frontUTR
    #cdsEndPos <- exEndPos
    
    
    strand <- dat$Strand[i]
    if (strand == -1) {
      frontUTR2 <- frontUTR
      frontUTR <- backUTR
      backUTR <- frontUTR2
    }
    
	if (trimUTR == F) {
      CDSLength <- exonLength
      frontUTR <- 0
      backUTR <- 0
    }
    
    
    #Figure out padding on either side if needed to get to minLength, otherwise trim to CDS
    if (CDSLength < minLength) {
      difference <- minLength - CDSLength
      if (difference %% 2 == 0) { # If even number
        frontRMposition <- martFlank + frontUTR - difference/2 + 1
        backRMposition <- frontRMposition + minLength - 1
      } else if (difference %% 2 == 1){ # If odd number
        difference1 <- difference+1
        frontRMposition <- martFlank + frontUTR - difference1/2 + 1
        backRMposition <- frontRMposition + minLength - 1
      }
    } else { 
      frontRMposition <- frontUTR + 100
      backRMposition <- frontRMposition + CDSLength
    }
    
    # output naming
    if (name == "GeneNameThenTransciptID" | name == "GeneNameThenGeneID") {
      geneName <- as.character(dat$Gene.name[i])
      if (is.na(geneName)) {
        if (name == "GeneNameThenTransciptID") {
          geneName <- dat$Transcript.stable.ID.version[i]  
        } else if (name == "GeneNameThenGeneID") {
          geneName <- dat$Gene.stable.ID[i]
        }
      }
    } else if (name == "TranscriptID") {
      geneName <- dat$Transcript.stable.ID.version[i]  
    } else if (name == "GeneID") {
      geneName <- dat$Gene.stable.ID[i]  
    } else if (name != "GeneName" | name != "TrancriptID" | name != "GeneID") {
      return(paste("Name argument",name,"needs to be ne of three options: 'GeneName', 'TranscriptID,"))
    }
    exonNumber <- as.character(dat$Exon.rank.in.transcript[i])
    output <- paste0(geneName,"_e",exonNumber)
    #filename <- paste0(wd,"/",output,".fas") # could use this if you want to write individual fasta files
    
    # clip the sequence
    dnaStringSplit <- strsplit(dnaString,split="",)[[1]]
    dnaStringSplitClip <- dnaStringSplit[frontRMposition:backRMposition]
    dna <- as.DNAbin(dnaStringSplitClip)
    dna <- as.list(dna)
    names(dna) <- output
    
    #append to fasta
    write.FASTA(dna, file=OUTPUT,append = T)
  }
}





trimFlankingSeqstoCDS <- function(dat,martFlank,OUTPUT) {
  
  # check if output file exists
  if (file.exists(OUTPUT)) {
    deleteQ <- readline("Output fasta file already exists. Would you like to remove before continuing? 
    (If no, output file will be appended) (y/n) ")
    if(deleteQ=="y") file.remove(OUTPUT)
  }
  
  for (i in 1:nrow(dat)) {
    # find sequence lengths, etc.
    exstart <- dat$Exon.region.start..bp.[i]
    exend <- dat$Exon.region.end..bp.[i]
    cdsstart <- dat$Genomic.coding.start[i]
    cdsend <- dat$Genomic.coding.end[i]
    
    exonLength <- exend - exstart + 1
    CDSLength <- cdsend - cdsstart + 1
    frontUTR <- cdsstart - exstart
    backUTR <- exend - cdsend
    
    dnaString <- dat$sequence[i] 
    seqLength <- nchar(dnaString)
    
    exStartPos <- martFlank + 1
    exEndPos <- seqLength - martFlank
    exLength <- exEndPos - exStartPos +1
    
    
    strand <- dat$Strand[i]
    if (strand == -1) {
      frontUTR2 <- frontUTR
      frontUTR <- backUTR
      backUTR <- frontUTR2
    }
    
    cdsStartPos <- exStartPos + frontUTR
    cdsEndPos <- exEndPos - backUTR
    
    
    # output naming
    if (name == "GeneNameThenTransciptID" | name == "GeneNameThenGeneID") {
      geneName <- as.character(dat$Gene.name[i])
      if (is.na(geneName)) {
        if (name == "GeneNameThenTransciptID") {
          geneName <- dat$Transcript.stable.ID.version[i]  
        } else if (name == "GeneNameThenGeneID") {
          geneName <- dat$Gene.stable.ID[i]
        }
      }
    } else if (name == "TranscriptID") {
      geneName <- dat$Transcript.stable.ID.version[i]  
    } else if (name == "GeneID") {
      geneName <- dat$Gene.stable.ID[i]  
    } else if (name != "GeneName" | name != "TrancriptID" | name != "GeneID") {
      return(paste("Name argument",name,"needs to be ne of three options: 'GeneName', 'TranscriptID,"))
    }
    #filename <- paste0(wd,"/",output,".fas") # could use this if you want to write individual fasta files
    
    # clip the sequence
    dnaStringSplit <- strsplit(dnaString,split="",)[[1]]
    dnaStringSplitClip <- dnaStringSplit[cdsStartPos:cdsEndPos]
    dna <- as.DNAbin(dnaStringSplitClip)
    dna <- as.list(dna)
    names(dna) <- output
    
    #append to fasta
    write.FASTA(dna, file=OUTPUT,append = T)
  }
}









#################################################################################
#################################################################################

#######################Fasta to Tabular format##################################

FastaToTabular <- function (filename){
  
  #read fasta file
  
  file1 <- readLines(filename)
  
  #find the genename location by grepping >
  
  location <- which((str_sub(file1,1,1))==">")
  
  #start an empty vector to collect name and sequence 
  
  name=c()
  sequence =c()
  
  
  
  #number of genes= number of loops
  #extract name first
  for ( i in 1:length(location)){
    name_line = location[i]
    name1 = file1[name_line]
    name=c(name,name1)
    #extract sequence between the names
    #the last sequence will be missed using this strategy 
    #so, we are using if condition to extract last sequence 
    start= location[i]+1
    end = location[i+1]-1
    if ( i < length (location)){
      
      end=end
      
    } else {
      
      end=length(file1)
    }
    
    lines = start:end
    sequence1= as.character(paste(file1[lines],collapse = ""))
    sequence =c(sequence,sequence1)
  }
  
  #now create table using name and sequence vector 
  
  data <- tibble(name,sequence)
  
  
  
  #finally export the file 
  #before that remove preexisting file
  unlink(c("dna_table.csv"),force=TRUE)
  write.csv(data,"dna_table.csv")
  
  #function ends
}

#usage
#FastaToTabular("dna_fasta.fasta")