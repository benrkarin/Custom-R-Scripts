
library(ape)

setwd("~/Dropbox/Berkeley_Research/MTGENOME_ANALYSES/Analyses/demux-lima-withadapter/25may2022-remove-trash-seqs/create-mtgenome-concat/partitioned")


full.284 <- read.tree("./284-bothfrags/mtgenome-bothfrags-parts.txt.treefile")
frag1.284 <- read.tree("./284-frag1/frag1-parts.txt.treefile")
nd2.284 <- read.tree("./284-nd2-only/parts-nd2.txt.treefile")

mtgenomesupports <- full.284$node.label
frag1supports <- frag1.284$node.label
nd2supports <- nd2.284$node.label


mtgenomealrt <- c()
mtgenomeufboots <- c()
for (i in 1:length(mtgenomesupports)) {
  spl <- strsplit(mtgenomesupports[i],split='/')[[1]]
  mtgenomealrt[i] <- spl[1]
  mtgenomeufboots[i] <- spl[2]
}

frag1alrt <- c()
frag1ufboots <- c()
for (i in 1:length(frag1supports)) {
  spl <- strsplit(frag1supports[i],split='/')[[1]]
  frag1alrt[i] <- spl[1]
  frag1ufboots[i] <- spl[2]
}

nd2alrt <- c()
nd2ufboots <- c()
for (i in 1:length(nd2supports)) {
  spl <- strsplit(nd2supports[i],split='/')[[1]]
  nd2alrt[i] <- spl[1]
  nd2ufboots[i] <- spl[2]
}


mtgenomealrt.rmna <- as.numeric(mtgenomealrt)[is.na(as.numeric(mtgenomealrt))==F]
mtgenomeufboots.rmna <- as.numeric(mtgenomeufboots)[is.na(as.numeric(mtgenomeufboots))==F]
frag1alrt.rmna <- as.numeric(frag1alrt)[is.na(as.numeric(frag1alrt))==F]
frag1ufboots.rmna <- as.numeric(frag1ufboots)[is.na(as.numeric(frag1ufboots))==F]
nd2alrt.rmna <- as.numeric(nd2alrt)[is.na(as.numeric(nd2alrt))==F]
nd2ufboots.rmna <- as.numeric(nd2ufboots)[is.na(as.numeric(nd2ufboots))==F]

a <- data.frame(group="ND2 ALRT", support=nd2alrt.rmna)
b <- data.frame(group="ND2 UFBoot", support=nd2ufboots.rmna)
c <- data.frame(group="Fragment 1 ALRT", support=frag1alrt.rmna)
d <- data.frame(group="Fragment 1 UFBoot", support=frag1ufboots.rmna)
e <- data.frame(group="MtGenome 1 ALRT", support=mtgenomealrt.rmna)
f <- data.frame(group="MtGenome 1 UFBoot", support=mtgenomeufboots.rmna)

test.data <- rbind(d,f)
t.test(support~group,data=test.data)

plot.data <- rbind(a,b,c,d,e,f)
plot.data$group <- factor(plot.data$group, levels =  c("ND2 UFBoot","Fragment 1 UFBoot","MtGenome 1 UFBoot","ND2 ALRT","Fragment 1 ALRT","MtGenome 1 ALRT"),ordered = TRUE)

plot.data2 <- rbind(b,d,f)
plot.data2$group <- factor(plot.data2$group, levels =  c("ND2 UFBoot","Fragment 1 UFBoot","MtGenome 1 UFBoot"),ordered = TRUE)



boxplot(support~group,data=plot.data)


pdf(file="supportvalues-partitioned.pdf",height=4,width = 6.5)
ggplot(plot.data, aes(x=group, y=support, fill = group)) + 
  geom_violin() +
  theme_classic() +
  #scale_color_manual(values=c('#999999','#E69F00')) +
  scale_color_brewer(palette = "PuOr") + 
  xlab("Group") +
  ylab("Support value (ALRT or UFBoot)") +
  #scale_y_log10() +
  stat_summary(fun=mean, geom="point", shape=23, size=2) +
  #stat_summary(fun=median, geom="point", size=2, color="red") +
  theme(legend.position="none")

dev.off()
