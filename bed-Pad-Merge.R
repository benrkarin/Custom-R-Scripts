#Pad Bed coordinates

## DEPENDENCIES
#Must have bedtools installed in path (or adjust line 19)

#USAGE
#setwd("~/PATH/TO/WORK/DIR")
#bedToolsMergethenPad(bed="significant_cor_regions.bed", mergeDistance = 400, shortExonLength = 200, longExonBuffer = 0)

## INPUT - a bed file with a name column (4th) and a strandedness column (6th)
# EXAMPLE:
#NC_056522.1	296640	296712	MPV17_exon_1	NA	-
#NC_056522.1	278635	278750	MPV17_exon_2	NA	-
#NC_056522.1	278429	278521	MPV17_exon_3	NA	-
#NC_056522.1	392565	392649	KHK_exon_1	NA	-
#NC_056522.1	387525	387669	KHK_exon_2	NA	-

#OUTPUT
# same name as input bud with "-merged-padded" added to name

# NEXT STEPS
# use getFastaFromBed() to extract the fasta from the new bed file. See getFastaFromBed.R
# getFastaFromBed(genome = "path/to/genome", bed = bedfile.bed)



bedToolsMergethenPad <- function(bed, mergeDistance = 400, shortExonLength = 200, longExonBuffer = 0, reorder=F) {
  filenamenoext <- strsplit(bed,".bed")[[1]]
  out1 <- paste0(filenamenoext,"-sorted.bed")
  out2 <- paste0(filenamenoext,"-merged.bed")

  
  c1 <- paste0("sort -k1,1 -k2,2n ", bed," > ",out1)
  system(c1)
  c2 <- paste0("bedtools merge -i ", out1," -d 400 -c 4,5,6 -o first,first,first -bed >", out2)
  system(c2)

  

  bedfile <- read.table(out2, header=F)
  bedfile <- bedfile[order(bedfile$V1,bedfile$V2),]

  rownames(bedfile) = NULL
  newbedfile <- bedfile
  deleterows <- c()
  for (i in 1:nrow(bedfile)) {
    startpos <- bedfile[i,2]
    endpos <- bedfile[i,3]
    li <- endpos - startpos + 1
   
    # Padding if below minlength
    if (li <= shortExonLength) {
      diff <- shortExonLength - li + 1
      if (diff %% 2 == 0) {
        pad <- diff/2
        newstartpos <-  startpos - pad
        newendpos <- endpos + pad
      } else if (diff %% 2 == 1) {
        pad <- (diff + 1)/ 2
        newstartpos <-  startpos - pad
        newendpos <- endpos + pad - 1
        newendpos-newstartpos
      }
    } else if (li > shortExonLength) {
      newstartpos <- startpos - longExonBuffer
      newendpos <- endpos + longExonBuffer
    }

      newbedfile[i,2] <- newstartpos
      newbedfile[i,3] <- newendpos

      bf <- newbedfile
  }
      if (reorder == T) {
        
        num <- rep(NA,nrow(bf))
        ex <- rep(NA,nrow(bf))
        for( i in 1:nrow(bf)) {
          s <- strsplit(bf[i,4],"_")
          num[i] <- s[[1]][3]
          ex[i] <- s[[1]][1]
        }
        num <- as.numeric(num)
        
        
        bf2 <- data.frame(bf,ex,num)
        bf3 <- bf2[order(bf2$ex,bf2$num),]
        #bf2$V1,
        newexnums <- c()
        exunique <- levels(factor(bf3$ex))
        for (j in 1:length(exunique)) {
          nmatches <- sum(bf3$ex %in% exunique[j])
          newexnums <- append(newexnums,seq(1:nmatches))
        }
        newexonnames <- c()
        for (i in 1:length(bf3$ex)) {
          newexonnames[i] <- paste(bf3$ex[i],"exon",as.character(newexnums)[i],sep="_")
        }
        bf4 <- data.frame(bf3,newexonnames,newexnums)
        bf4$V4 <- newexonnames
        bf5 <- bf4[order(bf4$V1,bf4$ex,newexnums),]
        
        outfile <- paste0(filenamenoext,"-merged-padded-reordered.bed")
        
        write.table(bf5[1:6], file=outfile, sep="\t", quote=F, col.names=F, row.names = F)
      } else if (reorder ==F) {
        outfile <- paste0(filenamenoext,"-merged-padded.bed")
        write.table(bf, file=outfile, sep="\t", quote=F, col.names=F, row.names = F)
      }
  
}



