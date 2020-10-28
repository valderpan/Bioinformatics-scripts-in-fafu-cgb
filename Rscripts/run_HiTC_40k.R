library(HiTC)
require(Matrix)
indir="/data/Schedule/20180927_Lesson12_Fengyanjie/Data/raw/40000/"
outdir="/data/Schedule/20180927_Lesson12_Fengyanjie/Analysis/40k/result/"

data_in=importC(paste(indir,"sam1_40000.matrix",sep="/"),paste(indir,"sam1_40000_abs.bed",sep="/"))

##Visualization of Interaction Maps
print("##Visualization of Interaction Maps")
pdf_file=paste(outdir,"mapC_all.pdf",sep = "/")
pdf(pdf_file)
mapC(data_in)
dev.off()

## Extract region of interest and plot the interaction map
print("## Extract region of interest and plot the interaction map")
hic_chr2L <- extractRegion(data_in$chr2Lchr2L,
 chr="chr2L",from=0,to=4000000)
pdf_file=paste(outdir,"Hi-C_interaction_map_of_chromosome_2L.pdf",sep="/")
pdf(pdf_file)
mapC(HTClist(hic_chr2L), maxrange=100)
dev.off()

## Data Normalization by Expected number of Counts
print("## Data Normalization by Expected number of Counts")
hic_chr2L_norm <- normPerExpected(hic_chr2L, method="loess")
pdf_file=paste(outdir,"Interaction_map_of_chromosome_2L_normalized_by_the_expected_interaction_counts.pdf",sep="/")
pdf(pdf_file)
mapC(HTClist(hic_chr2L_norm), log.data=TRUE)
dev.off()

## Correlation Map of Chromosome 2L 
print("## Correlation Map of Chromosome 2L")
#intdata(data_in14norm) <- as(cor(as.matrix(intdata(data_in14norm))),"Matrix")
intdata(hic_chr2L_norm) <- HiTC:::sparseCor(intdata(hic_chr2L_norm))
pdf_file=paste(outdir,"Correlation_map_of_chromosome_2L.pdf",sep="/")
pdf(pdf_file)
mapC(HTClist(hic_chr2L_norm), maxrange=1, minrange=-1,
 col.pos=c("black", "red"), col.neg=c("blue","black"))
dev.off()

## Principal Component Analysis
print("## Principal Component Analysis")
pc <- pca.hic(hic_chr2L_norm, normPerExpected=TRUE, method="loess", npc=1)
pdf_file=paste(outdir,"PCA_analysis_chr2L.pdf",sep="/")
pdf(pdf_file)
plot(start(pc$PC1), score(pc$PC1), type="h", 
 xlab="chr2L", ylab="PC1vec", frame=FALSE)
dev.off()



