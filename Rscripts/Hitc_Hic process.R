### R code from vignette source 'HiC_analysis.Rnw'

###################################################
### code chunk number 1: style
###################################################
BiocStyle::latex(use.unsrturl=FALSE)


###################################################
### code chunk number 2: head
###################################################
require(HiTC)
require(HiCDataHumanIMR90)
data(Dixon2012_IMR90)

show(hic_imr90_40)
class(intdata(hic_imr90_40$chr1chr1))
object.size(hic_imr90_40)


###################################################
### code chunk number 3: describe
###################################################
## Show data
show(hic_imr90_40)
## Is my data complete (i.e. composed of intra + the related inter chromosomal maps)
isComplete(hic_imr90_40)
## Note that a complete object is not necessarily pairwise 
## (is both chr1-chr2 and chr2-chr1 stored ?)
isPairwise(hic_imr90_40)
## Which chromosomes ?
seqlevels(hic_imr90_40)
## Details about a given map
detail(hic_imr90_40$chr6chr6)
## Descriptive statistics
head(summary(hic_imr90_40))


###################################################
### code chunk number 4: plot1
###################################################
## Go back to a smaller dataset (chr21, 22, X) at lower resolution
sset <- reduce(hic_imr90_40, chr=c("chr5","chr6","chr7"))
imr90_500 <- HTClist(mclapply(sset, binningC, 
                              binsize=500000, bin.adjust=FALSE, method="sum", step=1))
mapC(imr90_500)


###################################################
### code chunk number 5: plot2
###################################################
mapC(forcePairwise(imr90_500), maxrange=200)


###################################################
### code chunk number 6: resFrag
###################################################
## Example on chromosome X
## GRanges of restriction fragments after HindIII digestion
resFrag <- getRestrictionFragmentsPerChromosome(resSite="AAGCTT", overhangs5=1, 
                                                chromosomes="chr6", 
                                                genomePack="BSgenome.Hsapiens.UCSC.hg18")
resFrag


###################################################
### code chunk number 7: annot
###################################################
## Annotation of genomic features for LGF normalization
## Example on chromosome 6

## Load mappability track
require(rtracklayer)
##map_hg18 <- import("wgEncodeCrgMapabilityAlign100mer_chr6.bw",format="BigWig")
map_hg18 <- NULL
cutSites <- getAnnotatedRestrictionSites(resSite="AAGCTT", overhangs5=1, 
                                         chromosomes="chr6", 
                                         genomePack="BSgenome.Hsapiens.UCSC.hg18", 
                                         wingc=200, mappability=map_hg18, winmap=500)

head(cutSites)
## Annotation of Hi-C object
imr90_500_chr6annot <- setGenomicFeatures(imr90_500$chr6chr6, cutSites)
y_intervals(imr90_500_chr6annot)


###################################################
### code chunk number 8: normLGF (eval = FALSE)
###################################################
## ## LGF normalization
## imr90_500_LGF <- normLGF(imr90_500_chr6annot)


###################################################
### code chunk number 9: normICE
###################################################
imr90_500_ICE <-  normICE(imr90_500, max_iter=10)
mapC(HTClist(imr90_500_ICE$chr6chr6), trim.range=.95, 
     col.pos=c("white", "orange", "red", "black"))


###################################################
### code chunk number 10: tads
###################################################
hox <- extractRegion(hic_imr90_40$chr6chr6, chr="chr6", from=50e6, to=58e6)
plot(hox, maxrange=50, col.pos=c("white", "orange", "red", "black"))


###################################################
### code chunk number 11: di
###################################################
di<-directionalityIndex(hox)
barplot(di, col=ifelse(di>0,"darkred","darkgreen"))


###################################################
### code chunk number 12: sessionInfo
###################################################
toLatex(sessionInfo(), locale=FALSE)








library(HiTC)
require(Matrix)
indir="D:\\Result\\TAD\\no_contig"
outdir="D:\\Result\\TAD\\no_contig"

data_in=importC(paste(indir,"H01.500000.matrix",sep="\\"),paste(indir,"H01_500000.bed",sep="\\"))

##Visualization of Interaction Maps
print("##Visualization of Interaction Maps")
pdf_file=paste(outdir,"mapC_all.pdf",sep = "/")
pdf(pdf_file)
mapC(data_in)
dev.off()

## Extract region of interest and plot the interaction map
print("## Extract region of interest and plot the interaction map")
hic_chr5A <- extractRegion(data_in$Chr5AChr5A,
                           chr="Chr5A",from=0,to=67500000)
pdf_file=paste(outdir,"Hi-C_interaction_map_of_chromosome_5A.raw.pdf",sep="/")
pdf(pdf_file)
mapC(HTClist(hic_chr5A), maxrange=100)
dev.off()

## Data Normalization by Expected number of Counts
print("## Data Normalization by Expected number of Counts")
hic_chr5A_norm <- normPerExpected(hic_chr5A, method="loess")
#hic_chr5A_norm <- getExpectedCounts(hic_chr5A, method="loess")
pdf_file=paste(outdir,"Interaction_map_of_chromosome_5A_raw_normalized_by_the_expected_interaction_counts.pdf",sep="/")
pdf(pdf_file)
mapC(HTClist(hic_chr5A_norm), log.data=TRUE)
dev.off()

## Correlation Map of Chromosome 2L 
print("## Correlation Map of Chromosome 5A")
#intdata(data_in14norm) <- as(cor(as.matrix(intdata(data_in14norm))),"Matrix")
intdata(hic_chr5A_norm) <- HiTC:::sparseCor(intdata(hic_chr5A_norm))
pdf_file=paste(outdir,"Correlation_map_of_chromosome_5A.raw.pdf",sep="/")
pdf(pdf_file)
mapC(HTClist(hic_chr5A_norm), maxrange=1, minrange=-1,
     col.pos=c("black", "red"), col.neg=c("blue","black"))
dev.off()

## Principal Component Analysis
print("## Principal Component Analysis")
pc <- pca.hic(hic_chr5A_norm, normPerExpected=TRUE, method="loess", npc=1)
pdf_file=paste(outdir,"PCA_analysis_chr5A.raw.pdf",sep="/")
pdf(pdf_file)
plot(start(pc$PC1), score(pc$PC1), type="h", 
     xlab="chr5A", ylab="PC1vec", frame=FALSE,col = ifelse(score(pc$PC1)>0,'darkred','darkgreen'),lwd=2)
dev.off()
barplot(score(pc$PC1),col = ifelse(score(pc$PC1)>0,'red','blue'))


###

library(HiTC)
#require(Matrix)
indir="D:\\Result\\TAD\\no_contig"
outdir="D:\\Result\\TAD\\no_contig"
setwd('D:\\Result\\TAD')
data_in=importC(paste(indir,"Nepal_500000.matrix",sep="\\"),paste(indir,"Nepal_500000_abs.bed",sep="\\"))
#data_in=importC("Nepal_500000.matrix","Nepal_500000_abs.bed")
hic_chr5A <- extractRegion(data_in$Chr5AChr5A,chr="Chr5A",from=0,to=68253030)
# hic_chr5A_norm <- normPerExpected(hic_chr5A, method="loess")
pc <- pca.hic(hic_chr5A, normPerExpected=TRUE, method="loess", npc=1)

plot(start(pc$PC1), score(pc$PC1), type="h", 
     xlab="chr5A", ylab="PC1vec", frame=FALSE,col = ifelse(score(pc$PC1)>0,'darkred','darkgreen'),lwd=2)

barplot(score(pc$PC1),col = ifelse(score(pc$PC1)>0,'darkred','darkgreen'))
#使用ggplot2画图
library(ggplot2)
df <- data.frame(start(pc$PC1),score(pc$PC1))
colnames(df) <- c('start','score')
ggplot(df,aes(start,score))+geom_bar(stat="identity",color = ifelse(score(pc$PC1)>0,'darkred','darkgreen'),fill=ifelse(score(pc$PC1)>0,'darkred','darkgreen'))+theme_classic()


##实际运行
hic_chr5A <- extractRegion(data_in$Chr5AChr5A,chr="Chr5A",from=0,to=68253030)
pc <- pca.hic(hic_chr5A, normPerExpected=TRUE, method="loess", npc=1)
pdf_file=paste(outdir,"PCA_analysis_chr5A.raw.pdf",sep="/")
pdf(pdf_file)
plot(start(pc$PC1), score(pc$PC1), type="h",
     xlab="chr5A", ylab="PC1vec", frame=FALSE,col = ifelse(score(pc$PC1)>0,'darkred','darkgreen'),lwd=2)
dev.off()

#barplot(score(pc$PC1),col = ifelse(score(pc$PC1)>0,'red','blue'))
hic_chr5B <- extractRegion(data_in$Chr5BChr5B,chr="Chr5B",from=0,to=65586031)
pc <- pca.hic(hic_chr5B, normPerExpected=TRUE, method="loess", npc=1)
pdf_file=paste(outdir,"PCA_analysis_chr5B.raw.pdf",sep="/")
pdf(pdf_file)
plot(start(pc$PC1), score(pc$PC1), type="h",
     xlab="chr5B", ylab="PC1vec", frame=FALSE,col = ifelse(score(pc$PC1)>0,'darkred','darkgreen'),lwd=2)
dev.off()

hic_chr5C <- extractRegion(data_in$Chr5CChr5C,chr="Chr5C",from=0,to=61805356)
pc <- pca.hic(hic_chr5C, normPerExpected=TRUE, method="loess", npc=1)
pdf_file=paste(outdir,"PCA_analysis_chr5C.raw.pdf",sep="/")
pdf(pdf_file)
plot(start(pc$PC1), score(pc$PC1), type="h",
     xlab="chr5C", ylab="PC1vec", frame=FALSE,col = ifelse(score(pc$PC1)>0,'darkred','darkgreen'),lwd=2)
dev.off()

hic_chr5D <- extractRegion(data_in$Chr5DChr5D,chr="Chr5D",from=0,to=85924935)
pc <- pca.hic(hic_chr5D, normPerExpected=TRUE, method="loess", npc=1)
pdf_file=paste(outdir,"PCA_analysis_chr5D.raw.pdf",sep="/")
pdf(pdf_file)
plot(start(pc$PC1), score(pc$PC1), type="h",
     xlab="chr5D", ylab="PC1vec", frame=FALSE,col = ifelse(score(pc$PC1)>0,'darkred','darkgreen'),lwd=2)
dev.off()

hic_chr8A <- extractRegion(data_in$Chr8AChr8A,chr="Chr8A",from=0,to=42784806)
pc <- pca.hic(hic_chr8A, normPerExpected=TRUE, method="loess", npc=1)
pdf_file=paste(outdir,"PCA_analysis_chr8A.raw.pdf",sep="/")
pdf(pdf_file)
plot(start(pc$PC1), score(pc$PC1), type="h",
     xlab="chr8A", ylab="PC1vec", frame=FALSE,col = ifelse(score(pc$PC1)>0,'darkred','darkgreen'),lwd=2)
dev.off()

hic_chr8B <- extractRegion(data_in$Chr8BChr8B,chr="Chr8B",from=0,to=54634649)
pc <- pca.hic(hic_chr8B, normPerExpected=TRUE, method="loess", npc=1)
pdf_file=paste(outdir,"PCA_analysis_chr8B.raw.pdf",sep="/")
pdf(pdf_file)
plot(start(pc$PC1), score(pc$PC1), type="h",
     xlab="chr8B", ylab="PC1vec", frame=FALSE,col = ifelse(score(pc$PC1)>0,'darkred','darkgreen'),lwd=2)
dev.off()

hic_chr8C <- extractRegion(data_in$Chr8CChr8C,chr="Chr8C",from=0,to=41941943)
pc <- pca.hic(hic_chr8C, normPerExpected=TRUE, method="loess", npc=1)
pdf_file=paste(outdir,"PCA_analysis_chr8C.raw.pdf",sep="/")
pdf(pdf_file)
plot(start(pc$PC1), score(pc$PC1), type="h",
     xlab="chr8C", ylab="PC1vec", frame=FALSE,col = ifelse(score(pc$PC1)>0,'darkred','darkgreen'),lwd=2)
dev.off()

hic_chr8D <- extractRegion(data_in$Chr8DChr8D,chr="Chr8D",from=0,to=61380726)
pc <- pca.hic(hic_chr8D, normPerExpected=TRUE, method="loess", npc=1)
pdf_file=paste(outdir,"PCA_analysis_chr8D.raw.pdf",sep="/")
pdf(pdf_file)
plot(start(pc$PC1), score(pc$PC1), type="h",
     xlab="chr8D", ylab="PC1vec", frame=FALSE,col = ifelse(score(pc$PC1)>0,'darkred','darkgreen'),lwd=2)
dev.off()

# 对compartment的score进行绘图
df <- read.table('Nipor_compartment_counts.txt',sep = '\t',header = TRUE)
p=ggplot(df,aes(x=chr,y=score,fill=compartment))+geom_boxplot()+scale_fill_manual(values=c("darkred", "darkgreen"))+labs(fill="")+theme_bw()+theme(legend.position = 'top' ,legend.text = element_text(size = 10,face="bold"))+xlab('')+ylab('compartment score')

#AP85-boxplot
p=ggplot(df,aes(x=chr,y=score,fill=compartment))+geom_boxplot()+scale_fill_manual(values=c("darkred", "darkgreen"))+labs(fill="")+theme_bw()+theme(legend.position = 'top' ,axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),legend.text = element_text(size = 10,face="bold"))+xlab('')+ylab('compartment score')



#######################################
#2020.7.02  sorghum compartment analysis
########################################
library(HiTC)
library(ggplot2)
indir = 'C:\\Users\\dell\\Desktop\\Nepal_project\\20200702 Sobic 3D genome'
outdir = indir
data_in=importC(paste(indir,"sorghum_500000_iced.matrix",sep="\\"),paste(indir,"sorghum_500000_abs.bed",sep="\\"))
head(summary(data_in))
# 以下为重复内容
hic_chr01 <- extractRegion(data_in$Chr01Chr01,chr="Chr01",from=0,to=80884392)
pc01 <- pca.hic(hic_chr01, normPerExpected=TRUE, method="loess", npc=1)
write.csv(pc01$PC1,'Sorghum_chr01.compartment.csv')
pc01_df <- as.data.frame(pc01)
pc01_df <- mutate(pc01_df,position=ifelse(start<10,start,start/1000000))
pc01_df[1,9] = 0
ggplot(pc01_df,aes(x=position,y=score,fill = ifelse(score>0,'blue','yellow')))+
  geom_bar(stat = 'identity',width=0.5)+
  xlab('Chr01')+ylab('')+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank())+
  theme_classic()+
  theme(legend.position="none")+
  scale_fill_manual(values=c("orange", "#69b3a2"))
ggsave('Sorghum_Chr01.pdf')

hic_chr02 <- extractRegion(data_in$Chr02Chr02,chr="Chr02",from=0,to=77742459)
pc02 <- pca.hic(hic_chr02, normPerExpected=TRUE, method="loess", npc=1)
write.csv(pc02$PC1,'Sorghum_chr02.compartment.csv')
pc02_df <- as.data.frame(pc02)
pc02_df <- mutate(pc02_df,position=ifelse(start<10,start,start/1000000))
pc02_df[1,9] = 0
ggplot(pc02_df,aes(x=position,y=score,fill = ifelse(score>0,'blue','yellow')))+
  geom_bar(stat = 'identity',width=0.5)+
  xlab('Chr02')+ylab('')+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank())+
  theme_classic()+
  theme(legend.position="none")+
  scale_fill_manual(values=c("orange", "#69b3a2"))
ggsave('Sorghum_Chr02.pdf')

hic_chr03 <- extractRegion(data_in$Chr03Chr03,chr="Chr03",from=0,to=74386277)
pc03 <- pca.hic(hic_chr03, normPerExpected=TRUE, method="loess", npc=1)
write.csv(pc03$PC1,'Sorghum_chr03.compartment.csv')
pc03_df <- as.data.frame(pc03)
pc03_df <- mutate(pc03_df,position=ifelse(start<10,start,start/1000000))
pc03_df[1,9] = 0
ggplot(pc03_df,aes(x=position,y=score,fill = ifelse(score>0,'blue','yellow')))+
  geom_bar(stat = 'identity',width=0.5)+
  xlab('Chr03')+ylab('')+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank())+
  theme_classic()+
  theme(legend.position="none")+
  scale_fill_manual(values=c("orange", "#69b3a2"))
ggsave('Sorghum_Chr03.pdf')

hic_chr04 <- extractRegion(data_in$Chr04Chr04,chr="Chr04",from=0,to=68658214)
pc04 <- pca.hic(hic_chr04, normPerExpected=TRUE, method="loess", npc=1)
write.csv(pc04$PC1,'Sorghum_chr04.compartment.csv')
pc04_df <- as.data.frame(pc04)
pc04_df <- mutate(pc04_df,position=ifelse(start<10,start,start/1000000))
pc04_df[1,9] = 0
ggplot(pc04_df,aes(x=position,y=score,fill = ifelse(score>0,'blue','yellow')))+
  geom_bar(stat = 'identity',width=0.5)+
  xlab('Chr04')+ylab('')+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank())+
  theme_classic()+
  theme(legend.position="none")+
  scale_fill_manual(values=c("orange", "#69b3a2"))
ggsave('Sorghum_Chr04.pdf')

hic_chr05 <- extractRegion(data_in$Chr05Chr05,chr="Chr05",from=0,to=71854669)
pc05 <- pca.hic(hic_chr05, normPerExpected=TRUE, method="loess", npc=1)
write.csv(pc05$PC1,'Sorghum_chr05.compartment.csv')
pc05_df <- as.data.frame(pc05)
pc05_df <- mutate(pc05_df,position=ifelse(start<10,start,start/1000000))
pc05_df[1,9] = 0
ggplot(pc05_df,aes(x=position,y=score,fill = ifelse(score>0,'blue','yellow')))+
  geom_bar(stat = 'identity',width=0.5)+
  xlab('Chr05')+ylab('')+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank())+
  theme_classic()+
  theme(legend.position="none")+
  scale_fill_manual(values=c("orange", "#69b3a2"))
ggsave('Sorghum_Chr05.pdf')

hic_chr06 <- extractRegion(data_in$Chr06Chr06,chr="Chr06",from=0,to=61277060)
pc06 <- pca.hic(hic_chr06, normPerExpected=TRUE, method="loess", npc=1)
write.csv(pc06$PC1,'Sorghum_chr06.compartment.csv')
pc06_df <- as.data.frame(pc06)
pc06_df <- mutate(pc06_df,position=ifelse(start<10,start,start/1000000))
pc06_df[1,9] = 0
ggplot(pc06_df,aes(x=position,y=score,fill = ifelse(score>0,'blue','yellow')))+
  geom_bar(stat = 'identity',width=0.5)+
  xlab('Chr06')+ylab('')+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank())+
  theme_classic()+
  theme(legend.position="none")+
  scale_fill_manual(values=c("orange", "#69b3a2"))
ggsave('Sorghum_Chr06.pdf')

hic_chr07 <- extractRegion(data_in$Chr07Chr07,chr="Chr07",from=0,to=65505356)
pc07 <- pca.hic(hic_chr07, normPerExpected=TRUE, method="loess", npc=1)
write.csv(pc07$PC1,'Sorghum_chr07.compartment.csv')
pc07_df <- as.data.frame(pc07)
pc07_df <- mutate(pc07_df,position=ifelse(start<10,start,start/1000000))
pc07_df[1,9] = 0
ggplot(pc07_df,aes(x=position,y=score,fill = ifelse(score>0,'blue','yellow')))+
  geom_bar(stat = 'identity',width=0.5)+
  xlab('Chr07')+ylab('')+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank())+
  theme_classic()+
  theme(legend.position="none")+
  scale_fill_manual(values=c("orange", "#69b3a2"))
ggsave('Sorghum_Chr07.pdf')

hic_chr08 <- extractRegion(data_in$Chr08Chr08,chr="Chr08",from=0,to=62686529)
pc08 <- pca.hic(hic_chr08, normPerExpected=TRUE, method="loess", npc=1)
write.csv(pc08$PC1,'Sorghum_chr08.compartment.csv')
pc08_df <- as.data.frame(pc08)
pc08_df <- mutate(pc08_df,position=ifelse(start<10,start,start/1000000))
pc08_df[1,9] = 0
ggplot(pc08_df,aes(x=position,y=score,fill = ifelse(score>0,'blue','yellow')))+
  geom_bar(stat = 'identity',width=0.5)+
  xlab('Chr08')+ylab('')+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank())+
  theme_classic()+
  theme(legend.position="none")+
  scale_fill_manual(values=c("orange", "#69b3a2"))
ggsave('Sorghum_Chr08.pdf')

hic_chr09 <- extractRegion(data_in$Chr09Chr09,chr="Chr09",from=0,to=59416394)
pc09 <- pca.hic(hic_chr09, normPerExpected=TRUE, method="loess", npc=1)
write.csv(pc09$PC1,'Sorghum_chr09.compartment.csv')
pc09_df <- as.data.frame(pc09)
pc09_df <- mutate(pc09_df,position=ifelse(start<10,start,start/1000000))
pc09_df[1,9] = 0
ggplot(pc09_df,aes(x=position,y=score,fill = ifelse(score>0,'blue','yellow')))+
  geom_bar(stat = 'identity',width=0.5)+
  xlab('Chr09')+ylab('')+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank())+
  theme_classic()+
  theme(legend.position="none")+
  scale_fill_manual(values=c("orange", "#69b3a2"))
ggsave('Sorghum_Chr09.pdf')

hic_chr10 <- extractRegion(data_in$Chr10Chr10,chr="Chr10",from=0,to=61233695)
pc10 <- pca.hic(hic_chr10, normPerExpected=TRUE, method="loess", npc=1)
write.csv(pc10$PC1,'Sorghum_chr10.compartment.csv')
pc10_df <- as.data.frame(pc10)
pc10_df <- mutate(pc10_df,position=ifelse(start<10,start,start/1000000))
pc10_df[1,9] = 0
ggplot(pc10_df,aes(x=position,y=score,fill = ifelse(score>0,'blue','yellow')))+
  geom_bar(stat = 'identity',width=0.5)+
  xlab('Chr10')+ylab('')+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank())+
  theme_classic()+
  theme(legend.position="none")+
  scale_fill_manual(values=c("orange", "#69b3a2"))
ggsave('Sorghum_Chr10.pdf')


#######################################
#2020.7.14  Os compartment analysis
########################################
library(HiTC)
library(ggplot2)
library(dplyr)
setwd('C:\\Users\\dell\\Desktop\\Nepal_project\\20200714 Os_3D_genome\\compartment')
indir = 'C:\\Users\\dell\\Desktop\\Nepal_project\\20200714 Os_3D_genome'
data_in=importC(paste(indir,"Os_500000_iced.matrix",sep="\\"),paste(indir,"Os_500000_abs.bed",sep="\\"))
#以下为重复任务

hic_chr1 <- extractRegion(data_in$Chr1Chr1,chr="Chr1",from=0,to=45064769)
pc1 <- pca.hic(hic_chr1, normPerExpected=TRUE, method="loess", npc=1)
write.csv(pc1$PC1,'Os_chr1.compartment.csv')
pc1_df <- as.data.frame(pc1)
pc1_df <- mutate(pc1_df,position=ifelse(start<10,start,start/1000000))
pc1_df[1,9] = 0
ggplot(pc1_df,aes(x=position,y=score,fill = ifelse(score>0,'blue','yellow')))+
  geom_bar(stat = 'identity',width=0.5)+
  xlab('Chr1')+ylab('')+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank())+
  theme_classic()+
  theme(legend.position="none")+
  scale_fill_manual(values=c("orange", "#69b3a2"))
ggsave('Os_chr1.pdf')	

hic_chr2 <- extractRegion(data_in$Chr2Chr2,chr="Chr2",from=0,to=36823111)
pc2 <- pca.hic(hic_chr2, normPerExpected=TRUE, method="loess", npc=1)
write.csv(pc2$PC1,'Os_chr2.compartment.csv')
pc2_df <- as.data.frame(pc2)
pc2_df <- mutate(pc2_df,position=ifelse(start<10,start,start/1000000))
pc2_df[1,9] = 0
ggplot(pc2_df,aes(x=position,y=score,fill = ifelse(score>0,'blue','yellow')))+
  geom_bar(stat = 'identity',width=0.5)+
  xlab('Chr2')+ylab('')+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank())+
  theme_classic()+
  theme(legend.position="none")+
  scale_fill_manual(values=c("orange", "#69b3a2"))
ggsave('Os_chr2.pdf')	


hic_chr3 <- extractRegion(data_in$Chr3Chr3,chr="Chr3",from=0,to=37257345)
pc3 <- pca.hic(hic_chr3, normPerExpected=TRUE, method="loess", npc=1)
write.csv(pc3$PC1,'Os_chr3.compartment.csv')
pc3_df <- as.data.frame(pc3)
pc3_df <- mutate(pc3_df,position=ifelse(start<10,start,start/1000000))
pc3_df[1,9] = 0
ggplot(pc3_df,aes(x=position,y=score,fill = ifelse(score>0,'blue','yellow')))+
  geom_bar(stat = 'identity',width=0.5)+
  xlab('Chr3')+ylab('')+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank())+
  theme_classic()+
  theme(legend.position="none")+
  scale_fill_manual(values=c("orange", "#69b3a2"))
ggsave('Os_chr3.pdf')



hic_chr4 <- extractRegion(data_in$Chr4Chr4,chr="Chr4",from=0,to=35863200)
pc4 <- pca.hic(hic_chr4, normPerExpected=TRUE, method="loess", npc=1)
write.csv(pc4$PC1,'Os_chr4.compartment.csv')
pc4_df <- as.data.frame(pc4)
pc4_df <- mutate(pc4_df,position=ifelse(start<10,start,start/1000000))
pc4_df[1,9] = 0
ggplot(pc4_df,aes(x=position,y=score,fill = ifelse(score>0,'blue','yellow')))+
  geom_bar(stat = 'identity',width=0.5)+
  xlab('Chr4')+ylab('')+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank())+
  theme_classic()+
  theme(legend.position="none")+
  scale_fill_manual(values=c("orange", "#69b3a2"))
ggsave('Os_chr4.pdf')



hic_chr5 <- extractRegion(data_in$Chr5Chr5,chr="Chr5",from=0,to=30039014)
pc5 <- pca.hic(hic_chr5, normPerExpected=TRUE, method="loess", npc=1)
write.csv(pc5$PC1,'Os_chr5.compartment.csv')
pc5_df <- as.data.frame(pc5)
pc5_df <- mutate(pc5_df,position=ifelse(start<10,start,start/1000000))
pc5_df[1,9] = 0
ggplot(pc5df,aes(x=position,y=score,fill = ifelse(score>0,'blue','yellow')))+
  geom_bar(stat = 'identity',width=0.5)+
  xlab('Chr5')+ylab('')+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank())+
  theme_classic()+
  theme(legend.position="none")+
  scale_fill_manual(values=c("orange", "#69b3a2"))
ggsave('Os_chr5.pdf')	

hic_chr6 <- extractRegion(data_in$Chr6Chr6,chr="Chr6",from=0,to=32124789)
pc6 <- pca.hic(hic_chr6, normPerExpected=TRUE, method="loess", npc=1)
write.csv(pc6$PC1,'Os_chr6.compartment.csv')
pc6_df <- as.data.frame(pc6)
pc6_df <- mutate(pc6_df,position=ifelse(start<10,start,start/1000000))
pc6_df[1,9] = 0
ggplot(pc6_df,aes(x=position,y=score,fill = ifelse(score>0,'blue','yellow')))+
  geom_bar(stat = 'identity',width=0.5)+
  xlab('Chr6')+ylab('')+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank())+
  theme_classic()+
  theme(legend.position="none")+
  scale_fill_manual(values=c("orange", "#69b3a2"))
ggsave('Os_chr6.pdf')


hic_chr7 <- extractRegion(data_in$Chr7Chr7,chr="Chr7",from=0,to=30357780)
pc7 <- pca.hic(hic_chr7, normPerExpected=TRUE, method="loess", npc=1)
write.csv(pc7$PC1,'Os_chr7.compartment.csv')
pc7_df <- as.data.frame(pc7)
pc7_df <- mutate(pc7_df,position=ifelse(start<10,start,start/1000000))
pc7_df[1,9] = 0
ggplot(pc7_df,aes(x=position,y=score,fill = ifelse(score>0,'blue','yellow')))+
  geom_bar(stat = 'identity',width=0.5)+
  xlab('Chr7')+ylab('')+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank())+
  theme_classic()+
  theme(legend.position="none")+
  scale_fill_manual(values=c("orange", "#69b3a2"))
ggsave('Os_chr7.pdf')	


hic_chr8 <- extractRegion(data_in$Chr8Chr8,chr="Chr8",from=0,to=28530027)
pc8 <- pca.hic(hic_chr8, normPerExpected=TRUE, method="loess", npc=1)
write.csv(pc8$PC1,'Os_chr8.compartment.csv')
pc8_df <- as.data.frame(pc8)
pc8_df <- mutate(pc8_df,position=ifelse(start<10,start,start/1000000))
pc8_df[1,9] = 0
ggplot(pc8_df,aes(x=position,y=score,fill = ifelse(score>0,'blue','yellow')))+
  geom_bar(stat = 'identity',width=0.5)+
  xlab('Chr8')+ylab('')+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank())+
  theme_classic()+
  theme(legend.position="none")+
  scale_fill_manual(values=c("orange", "#69b3a2"))
ggsave('Os_chr8.pdf')	

hic_chr9 <- extractRegion(data_in$Chr9Chr9,chr="Chr9",from=0,to=23661561)
pc9 <- pca.hic(hic_chr9, normPerExpected=TRUE, method="loess", npc=1)
write.csv(pc9$PC1,'Os_chr9.compartment.csv')
pc9_df <- as.data.frame(pc9)
pc9_df <- mutate(pc9_df,position=ifelse(start<10,start,start/1000000))
pc9_df[1,9] = 0
ggplot(pc9_df,aes(x=position,y=score,fill = ifelse(score>0,'blue','yellow')))+
  geom_bar(stat = 'identity',width=0.5)+
  xlab('Chr9')+ylab('')+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank())+
  theme_classic()+
  theme(legend.position="none")+
  scale_fill_manual(values=c("orange", "#69b3a2"))
ggsave('Os_chr9.pdf')	

hic_chr10 <- extractRegion(data_in$Chr10Chr10,chr="Chr10",from=0,to=23661561)
pc10 <- pca.hic(hic_chr10, normPerExpected=TRUE, method="loess", npc=1)
write.csv(pc10$PC1,'Os_chr10.compartment.csv')
pc10_df <- as.data.frame(pc10)
pc10_df <- mutate(pc10_df,position=ifelse(start<10,start,start/1000000))
pc10_df[1,9] = 0
ggplot(pc10_df,aes(x=position,y=score,fill = ifelse(score>0,'blue','yellow')))+
  geom_bar(stat = 'identity',width=0.5)+
  xlab('Chr10')+ylab('')+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank())+
  theme_classic()+
  theme(legend.position="none")+
  scale_fill_manual(values=c("orange", "#69b3a2"))
ggsave('Os_chr10.pdf')	

hic_chr11 <- extractRegion(data_in$Chr11Chr11,chr="Chr11",from=0,to=30828668)
pc11 <- pca.hic(hic_chr11, normPerExpected=TRUE, method="loess", npc=1)
write.csv(pc11$PC1,'Os_chr11.compartment.csv')
pc11_df <- as.data.frame(pc11)
pc11_df <- mutate(pc11_df,position=ifelse(start<10,start,start/1000000))
pc11_df[1,9] = 0
ggplot(pc11_df,aes(x=position,y=score,fill = ifelse(score>0,'blue','yellow')))+
  geom_bar(stat = 'identity',width=0.5)+
  xlab('Chr11')+ylab('')+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank())+
  theme_classic()+
  theme(legend.position="none")+
  scale_fill_manual(values=c("orange", "#69b3a2"))
ggsave('Os_chr11.pdf')	

hic_chr12 <- extractRegion(data_in$Chr12Chr12,chr="Chr12",from=0,to=27757321)
pc12 <- pca.hic(hic_chr12, normPerExpected=TRUE, method="loess", npc=1)
write.csv(pc12$PC1,'Os_chr12.compartment.csv')
pc12_df <- as.data.frame(pc12)
pc12_df <- mutate(pc12_df,position=ifelse(start<10,start,start/1000000))
pc12_df[1,9] = 0
ggplot(pc12_df,aes(x=position,y=score,fill = ifelse(score>0,'blue','yellow')))+
  geom_bar(stat = 'identity',width=0.5)+
  xlab('Chr12')+ylab('')+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank())+
  theme_classic()+
  theme(legend.position="none")+
  scale_fill_manual(values=c("orange", "#69b3a2"))
ggsave('Os_chr12.pdf')