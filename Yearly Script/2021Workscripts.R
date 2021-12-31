
###--------------------------------------------------------------------
##2021.01.19 对Baidila八倍体甘蔗初步挂载后的结果进行统计  ##
##这里是将Badila分成了80个group，统计挂载后每个group的长度##
##--------------------------------------------------------------------

setwd('E:/Badila')
library(ggplot2)
library(dplyr)

df <- read.table('groups.asm.fasta.sizes')
df <- mutate(df,size=V2/1000000)
dd <- df[!grepl('tig',df$V1),]
ggplot(dd,aes(V1,size,fill=V1))+geom_bar(stat = 'identity',color='black')+
  theme_classic()+theme(legend.position="none")+
  theme(axis.text.x= element_text(angle=90,hjust=0))+
  scale_y_continuous(expand=c(0,0))+
  theme(axis.text=element_text(size=13),axis.title=element_text(size=20,face="bold"))+xlab('')

###--------------------------------------------------------------------
##2021.03.02 NpX to AP85 的compartment statu switch 的组合与更正 ##
##使用物种见的compartment上的2，5，6，7，8染色体的转变情况组合剩余的几条染色体在npx2AP85间的转变情况进行绘图##
##--------------------------------------------------------------------
library(ggplot2)
library(xlsx)
library(dplyr)
library(reshape2)
library(ggprism)
setwd('D:\\BaiduNetdiskDownload\\Nepal_project\\20200717 Nepal-AP85-compartment-variation')
dd2 <- read.xlsx('NpX2AP85.mergemcorrected.xlsx',sheetIndex = 2)
dd_2_melt <- melt(dd2,id.vars = 'Chr')
ggplot(dd_2_melt,aes(x=as.factor(Chr)))+
  geom_bar(aes(y=value,fill=variable),stat="identity",position="fill")+theme_bw()+
  xlab('Np-X vs AP85-441')+ylab('Ratio of ortholog regions')+theme(legend.position = 'top')


###--------------------------------------------------------------------
##2021.03.09 番茄转录组王刚 ##
##使用edgeR对转录组数据做差异表达##
##---------------------------------------------------------------------
library(edgeR)
library(dplyr)
library(openxlsx)
sampleNames <- c('G3','G6','A')
df <- read.table('E:\\王刚师兄数据\\counts_id.txt',header = TRUE,row.names = 1)
df <- df[,c(6,7,8,9)]
head(df)
dd <- mutate(df,a=rownames(df),tmp=(..A1.bam+..A2.bam)/2)
head(dd)
data <- select(dd,a,..G3.bam,..G6.bam,tmp)
rownames(data) <- data$a
data <- data[,-1]
colnames(data) <- sampleNames
#------G3 vs A----------------
data1 <- data[,c(1,3)]
data1 <- as.matrix(data1)
group_list <- factor(c(rep("Contral",1),rep("Treat",1)))
y <- DGEList(counts=data1, group = group_list)
y <- calcNormFactors(y)
et <- exactTest(y, dispersion = bcv ^ 2)
gene1 <- decideTestsDGE(et, p.value = 0.05, lfc = 0)
summary(gene1)
colnames(gene1) <- "Signifi"
results <- cbind(y$counts,et$table,gene1)
DEG <- results[results$Signifi!=0,]
up <- results[results$Signifi==1,]
down <- results[results$Signifi==-1,]
head(results)
result_name = paste(colnames(data1)[1],colnames(data1)[2],sep = '_vs_')
sheets <- list('DEG_overview_result'=DEG,'up_regulated_gene'=up,'down_regulated_gene'=down)
write.xlsx(sheets,paste(result_name,'.xlsx',sep = ''),row.names = T)

#Treat-Contral
#Down            3931
#NotSig         27465
#Up              4372

#------G3 vs G6--------------
data2 <- data[,c(1,2)]
data2 <- as.matrix(data2)
head(data2)
group_list <- factor(c(rep("Contral",1),rep("Treat",1)))
y <- DGEList(counts=data2, group = group_list)
y <- calcNormFactors(y)
et <- exactTest(y, dispersion = bcv ^ 2)
gene1 <- decideTestsDGE(et, p.value = 0.05, lfc = 0)
summary(gene1)
colnames(gene1) <- "Signifi"
results <- cbind(y$counts,et$table,gene1)
DEG <- results[results$Signifi!=0,]
up <- results[results$Signifi==1,]
down <- results[results$Signifi==-1,]
head(results)
result_name = paste(colnames(data2)[1],colnames(data2)[2],sep = '_vs_')
sheets <- list('DEG_overview_result'=DEG,'up_regulated_gene'=up,'down_regulated_gene'=down)
write.xlsx(sheets,paste(result_name,'.xlsx',sep = ''),row.names = T)

#Treat-Contral
#Down            1805
#NotSig         31431
#Up              2532

#GO/KEGG analysis
setwd('E:\\Slycopersicum_DEG')
go <- read.delim("Sl.allgene.anno.tab",sep='\t',header = TRUE)
head(go)
gene_select <- read.delim('G3vsA_DEG.txt',sep = '\t',header = FALSE)
gene_select <-as.vector(gene_select$V1)
go_enrich <- enricher(gene = gene_select,
  TERM2GENE = go[c('GO', 'gene')], 
  TERM2NAME = go[c('GO', 'description')], 
  pvalueCutoff = 0.05, 
  pAdjustMethod = 'BH')

write.xlsx(as.data.frame(go_enrich),"G3vsA_GO.xlsx",row.names = F)
dotplot(go_enrich,showCategory=20)
barplot(go_enrich,showCategory = 20)

Ko <- read.delim("KOannotation.tsv",sep='\t',header = TRUE)
ko_enrich <- enricher(gene = gene_select,
                      TERM2GENE = Ko[c('KO', 'gene')], 
                      TERM2NAME = Ko[c('KO', 'description')], 
                      pvalueCutoff = 0.1, 
                      pAdjustMethod = 'BH')

write.xlsx(as.data.frame(ko_enrich),"G3vsA_KO.xlsx",row.names = F)
dotplot(ko_enrich,showCategory=20)
barplot(ko_enrich,showCategory = 20)

gene_select <- read.delim('G3vsG6_DEG.txt',sep = '\t',header = F)
gene_select <-as.vector(gene_select$V1)

#03.12 对王刚师兄发送的GO文件结果进行富集
library(openxlsx)
library(ggplot2)
library(dplyr)
df <- read.xlsx('E:\\Slycopersicum_DEG\\DEG_GOKEGG\\GO_analysis_raw.xlsx')
dd <- df[order(-df$Number),]
ds <- dd[1:20,]
ggplot(ds,aes(y=Number,x=Description,fill=Padjust))+geom_bar(stat = 'identity')+coord_flip()+scale_fill_gradient(low="purple",high="red")

ggplot(dd[1:20,],aes(x=Description,y=Number,fill=Term.Type))+geom_bar(stat = 'identity',width = 0.8)+coord_flip()+scale_fill_manual(values = CPCOLS)
CPCOLS <- c("#8DA1CB", "#FD8D62", "#66C3A5")
labess = c('molecular_function','catalytic activity','transferase activity','response to stimulus','oxidoreductase activity','kinase activity',
           'membrane part','intrinsic component of membrane','integral component of membrane','membrane',
           'biological_process','metabolic process','organonitrogen compound metabolic process','biological regulation','protein metabolic process','regulation of cellular process','macromolecule modification',
           'cellular protein modification process','protein modification process','response to stimulus','phosphorus metabolic process')
ds_order <- ds[order(ds$Term.Type),]
ds_order$Description <- factor(ds_order$Description,levels=labess,ordered = TRUE)
ggplot(ds_order,aes(x=Description,y=Number,fill=Term.Type))+geom_bar(stat = 'identity',width = 0.8)+coord_flip()+scale_fill_manual(values = CPCOLS)+xlab('GeneNumber')+ylab('GO term')+theme_bw()


dd_MF <- dd %>% filter(Term.Type %in% 'MF')
dd_MF <- dd_MF[1:10,]
dd_BP <- dd %>% filter(Term.Type %in% 'BP')
dd_BP <- dd_BP[1:10,]
dd_CC <- dd %>% filter(Term.Type %in% 'CC')
dd_CC <- dd_CC[1:10,]
ds <- rbind(dd_MF,dd_BP,dd_CC)
ggplot(ds,aes(x=fct_reorder(Description,Term.Type),y=Number,fill=Term.Type))+geom_bar(stat = 'identity',width = 0.8)+coord_flip()+scale_fill_manual(values = CPCOLS)+xlab('GeneNumber')+ylab('GO term')+theme_bw()






###--------------------------------------------------------------------
##2021.04.02 ATAC-Seq Chipseeker注释##
##---------------------------------------------------------------------
library(GenomicFeatures)
library(ChIPseeker)
setwd('E:\\潘浩然\\调控元件生信任务\\analysis\\callpeak')
spompe_LA <- makeTxDbFromGFF('Soffic.v20191009.corrected.gff3')
spompe_Ss <- makeTxDbFromGFF('Sspon.v20190103.gff3')

setwd('E:\\潘浩然\\调控元件生信任务\\analysis\\Sugarcane\\peaks')
peak_LA <- readPeakFile('LA-pruple_summits.bed')
peak_Ss <- readPeakFile('SES208_peaks.narrowPeak.idr.txt')

peakAnno_LA <- annotatePeak(peak_LA,tssRegion = c(-3000,3000),TxDb = spompe_LA)
peakAnno_Ss <- annotatePeak(peak_Ss,tssRegion = c(-3000,3000),TxDb = spompe_Ss)
write.table(peakAnno_LA, file = 'LA-pruple_peak.txt',sep = '\t', quote = FALSE, row.names = FALSE)
write.table(peakAnno_Ss, file = 'SES208_peak.txt',sep = '\t', quote = FALSE, row.names = FALSE)

# plotAnnoPie(peakAnno_LA)
# plotAnnoBar(peakAnno_LA)
# plotAnnoPie(peakAnno_Ss)
# plotAnnoBar(peakAnno_Ss)

library(openxlsx)
library(dplyr)
peak_Ss <- readPeakFile('SES208_peaks.narrowPeak.txt')
peakAnno_Ss <- annotatePeak(peak_Ss,tssRegion = c(-2000,2000),TxDb = spompe_Ss)
write.table(peakAnno_Ss, file = 'SES208_v2018.2550peak.anno.txt',sep = '\t', quote = FALSE, row.names = FALSE)

df_anno <- read.xlsx('E:\\潘浩然\\调控元件生信任务\\analysis\\Sugarcane\\peaks\\SES208_peak.mutated.xlsx')
df_number <- df_anno %>% group_by(newanno) %>% summarise(number = n())

df_ratio <- df_number %>% arrange(desc(newanno)) %>% mutate(prop = number/sum(number)*100) %>% mutate(ypos = cumsum(prop-0.5*prop))

mycolors <- c("#3C989E","#C40003","#00C19B","#007ED3","#EAC862","#BA415D","#FF9D1E")
mycolors2 <- c("#2C3E50","#E31A1C","#18BC9C","#FF7F00","#A6CEE3","#B2DF8A","#1F78B4","#FB9A99","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A")

ggplot(df_ratio,aes(x="",y=prop,fill=newanno))+ 
  geom_bar(stat = 'identity',width = 1,color='white')+ 
  coord_polar("y", start=0)+theme_void()+scale_fill_manual(values = mycolors2)

#OsZmSb
setwd('E:\\潘浩然\\调控元件生信任务\\analysis\\Sugarcane\\OsZmSb')
spompe_Os <- makeTxDbFromGFF('Osativa_323_v7.0.gene.gff3')
spompe_Zm <- makeTxDbFromGFF()
spompe_Sb <- makeTxDbFromGFF()

# peak_Os <- readPeakFile('Maize_7days_leaf_ACRs.corrected.bed')
# peakAnno_Os <- annotatePeak(peak_Os,tssRegion = c(-3000,3000),TxDb = spompe_Os)
# write.table(peakAnno_Os, file = 'Os_peak.txt',sep = '\t', quote = FALSE, row.names = FALSE)
# peakanno_Os <- read.xlsx('Os_peak_mutated.xlsx')
# Os_number <- peakanno_Os %>% group_by(newanno) %>% summarise(number=n())
# Os_ratio <- Os_number %>% arrange(desc(newanno)) %>% mutate(prop = number/sum(number)*100) %>% mutate(ypos = cumsum(prop-0.5*prop))
# ggplot(Os_ratio,aes(x="",y=prop,fill=newanno))+ 
#   geom_bar(stat = 'identity',width = 1,color='white')+ 
#   coord_polar("y", start=0)+theme_void()+scale_fill_manual(values = mycolors2)
Os_anno <- read.table('Maize_7days_leaf_ACRs.corrected.bed',sep = '\t')
Os_number <- Os_anno %>% group_by(V5) %>% summarise(number=n())
Os_ratio <- Os_number %>% mutate(prop = number/sum(number)*100)

Zm_anno <- 


ggplot(Os_ratio,aes(x='',y=prop,fill=V5))+geom_bar(stat = 'identity',width = 0.5,color='black')


##---------------------------------------------------------------------
#2021.07.15 对Sb/Zm/Os peaks进行注释
##---------------------------------------------------------------------

library(GenomicFeatures)
library(ChIPseeker)
setwd('E:\\潘浩然\\调控元件生信任务\\analysis\\Sugarcane\\OsZmSb')
spompe_Sb <- makeTxDbFromGFF('Sbicolor_454_v3.1.1.gene.gff3')
spompe_Zm <- makeTxDbFromGFF('Zmays_493_RefGen_V4.gene.gff3')
spompe_Os <- makeTxDbFromGFF('Osativa_323_v7.0.gene.gff3')

peak_Sb <- readPeakFile('Sorghum_7days_leaf_ACRs.bed')
peak_Zm <- readPeakFile('Maize_7days_leaf_ACRs.bed')
peak_Os <- readPeakFile('Rice_7days_leaf_ACRs.bed')

peakAnno_Sb <- annotatePeak(peak_Sb,tssRegion = c(-2000,2000),TxDb = spompe_Sb)
peakAnno_Zm <- annotatePeak(peak_Zm,tssRegion = c(-2000,2000),TxDb = spompe_Zm)
peakAnno_Os <- annotatePeak(peak_Os,tssRegion = c(-2000,2000),TxDb = spompe_Os)
write.table(peakAnno_Sb, file = 'Sb.peak.anno.txt',sep = '\t',
            quote = FALSE, row.names = FALSE)
write.table(peakAnno_Zm, file = 'Zm.peak.anno.txt',sep = '\t',
            quote = FALSE, row.names = FALSE)
write.table(peakAnno_Os, file = 'Os.peak.anno.txt',sep = '\t',
            quote = FALSE, row.names = FALSE)







#2021.05.10 直接使用文章中Os\Zm\Sb的bed数据进行绘图
library(openxlsx)
library(tidyverse)
setwd('E:\\潘浩然\\调控元件生信任务\\analysis\\Sugarcane\\OsZmSb')
Zm_acr <- read.delim('Maize_7days_leaf_ACRs.corrected.bed',header = F,sep = '\t')
Os_acr <- read.delim('Rice_7days_leaf_ACRs.bed',sep = '\t',header = F)
Sb_acr <- read.delim('Sorghum_7days_leaf_ACRs.bed',sep = '\t',header = F)
Ses_acr <- read.xlsx('E:\\潘浩然\\调控元件生信任务\\analysis\\Sugarcane\\peaks\\SES208_peak_category.xlsx')
Ses_mono_acr <- read.xlsx('E:\\潘浩然\\调控元件生信任务\\analysis\\Sugarcane\\peaks\\monoploidy\\SES208_mono.peak.anno.mutate3category.xlsx')

##Os
# Os_acr <- Os_acr %>% filter(V1 %in% c(paste('Chr',1:12,sep = '')))
Os <- Os_acr %>% group_by(V5) %>% summarise(number=n()) %>% mutate(ratio = number/sum(number)) %>% mutate(species='Os')

##Zm
Zm_acr <- Zm_acr %>% filter(V1 %in% c(paste('Chr',1:10,sep = '')))
Zm <- Zm_acr %>% group_by(V5) %>% summarise(number=n()) %>% mutate(ratio = number/sum(number))%>% mutate(species='Zm')

##Sb
Sb_acr1 <- Sb_acr %>% filter(V1 %in% c(paste('Chr0',1:9,sep = '')))
Sb_acr2 <- Sb_acr %>% filter(V1 %in% ('Chr10'))
Sb_acr <- rbind(Sb_acr1,Sb_acr2)
Sb <- Sb_acr %>% group_by(V5) %>% summarise(number=n()) %>% mutate(ratio = number/sum(number))%>% mutate(species='Sb')

##SES
ses <- Ses_acr %>% group_by(category) %>% summarise(number=n()) %>% mutate(ratio = number/sum(number),species='SES')
ses <- ses %>% mutate(V5=category) %>% select(V5,number,ratio,species)

#SES mono
ses_mono <- Ses_mono_acr %>% group_by(category) %>% summarise(number=n()) %>% mutate(ratio=number/sum(number))

df <- rbind(Os,Zm,Sb,ses)
head(df)

mycolors <- c("#FF9D1E","#007ED3","#B2DF8A")
Os_pie <- ggplot(df %>% filter(species == 'Os'),aes(x='',y=ratio,fill=V5))+
  geom_bar(stat = 'identity',width = 1,color='white')+
  coord_polar("y", start=0)+theme_void()+
  scale_fill_manual(values = mycolors)+
  theme(legend.position="none")+
  ggtitle('Os')
Zm_pie <- ggplot(df %>% filter(species == 'Zm'),aes(x='',y=ratio,fill=V5))+
  geom_bar(stat = 'identity',width = 1,color='white')+
  coord_polar("y", start=0)+theme_void()+
  scale_fill_manual(values = mycolors)+
  theme(legend.position="none")+
  ggtitle('Zm')
Sb_pie <- ggplot(df %>% filter(species == 'Sb'),aes(x='',y=ratio,fill=V5))+
  geom_bar(stat = 'identity',width = 1,color='white')+
  coord_polar("y", start=0)+theme_void()+
  scale_fill_manual(values = mycolors)+
  theme(legend.position="none")+
  ggtitle('Sb')
ses_pie <- ggplot(df %>% filter(species == 'SES'),aes(x='',y=ratio,fill=V5))+
  geom_bar(stat = 'identity',width = 1,color='white')+
  coord_polar("y", start=0)+theme_void()+
  scale_fill_manual(values = mycolors)+
  theme(legend.position="none")+
  ggtitle('AP85-441')

ses_mono_pie <- ggplot(ses_mono,aes(x='',y=ratio,fill=category))+
  geom_bar(stat = 'identity',width = 1,color='white')+
  coord_polar("y", start=0)+theme_void()+
  scale_fill_manual(values = mycolors)+
  theme(legend.position="none")+ggtitle('AP85-441 mono')
library(gridExtra)
grid.arrange(Os_pie,Sb_pie,Zm_pie,ses_pie,ses_mono_pie
             ,nrow=2,ncol=3)


###ACR number\length\percentage
options(scipen = 999)
Os_total_number <- df %>% filter(species=='Os') %>% 
  mutate(total = sum(number)) %>% group_by(species) %>% select(species,total)
Os_total_number <- Os_total_number[1,]
Sb_total_number <- df %>% filter(species=='Sb') %>% 
  mutate(total = sum(number)) %>% group_by(species) %>% select(species,total)
Sb_total_number <- Sb_total_number[1,]
Zm_total_number <- df %>% filter(species=='Zm') %>% 
  mutate(total = sum(number)) %>% group_by(species) %>% select(species,total)
Zm_total_number <- Zm_total_number[1,]
SES_total_number <- df %>% filter(species=='SES') %>% 
  mutate(total = sum(number)) %>% group_by(species) %>% select(species,total)
SES_total_number <- SES_total_number[1,]


total_number_df <- rbind(Os_total_number,Sb_total_number,Zm_total_number,SES_total_number)
ggplot(total_number_df,aes(x=fct_relevel(species,"Os","Sb","Zm","SES"),y=total))+
  geom_bar(stat='identity',fill="#7cb9e8",width = 0.5)+
  xlab('')+ylab('ACR number')+
  theme_bw()+theme_prism(border = TRUE, base_rect_size = 1) +
  coord_cartesian(clip = "off")

Os_length <- Os_acr %>% mutate(length=V3-V2,species='Os') %>% mutate(total_length=sum(length))%>% select(species,total_length)
Os_length <- Os_length[1,]
Os_length <- Os_length %>% mutate(percent = total_length/374471240*100,Mb=total_length/1000000)

Sb_length <- Sb_acr %>% mutate(length=V3-V2,species='Sb') %>% mutate(total_length=sum(length))%>% select(species,total_length)
Sb_length <- Sb_length[1,]
Sb_length <- Sb_length %>% mutate(percent = total_length/708863705*100,Mb=total_length/1000000)

Zm_length <- Zm_acr %>% mutate(length=V3-V2,species='Zm') %>% mutate(total_length=sum(length))%>% select(species,total_length)
Zm_length <- Zm_length[1,]
Zm_length <- Zm_length %>% mutate(percent = total_length/2135083061*100,Mb=total_length/1000000)

ses_length <- Ses_acr %>% mutate(length=end-start,species='SES') %>% mutate(total_length=sum(length))%>% select(species,total_length)
ses_length <- ses_length[1,]
ses_length <- ses_length %>% mutate(percent = total_length/3140615268*100,Mb=total_length/1000000)


length_df <- rbind(Os_length,Sb_length,Zm_length,ses_length)
ggplot(length_df,aes(x=fct_relevel(species,"Os","Sb","Zm","SES"),y=Mb))+
  geom_bar(stat='identity',fill="#ff9966",width = 0.5)+theme_bw()+
  theme_prism(border = TRUE, base_rect_size = 1) +
  ylab('Total length of ACRs(Mb)')+xlab('')+
  ylim(0,15)+
  coord_cartesian(clip = "off")

ggplot(length_df,aes(x=fct_relevel(species,"Os","Sb","Zm","SES"),y=percent))+
  geom_bar(stat='identity',fill="#ffbf00",width = 0.5)+theme_bw()+
  geom_hline(yintercept = 1.0,color='black',size=1.2,linetype=5)+
  xlab('')+ylab('Percentage of total genome')+
  ylim(0,2)+
  theme_prism(border = TRUE, base_rect_size = 1) +
  coord_cartesian(clip = "off")


###--------------------------------------------------------------------
##2021.04.20 NpX文章 BiTriAll=conserved enrichment                   ##
##---------------------------------------------------------------------
library(ggplot2)
library(clusterProfiler)

setwd('E:\\Nepal_project\\20210420_Article_Revision\\BioTriAll_conserverd_enrichment')

#GO
go_anno <- read.delim('E:\\Nepal_project\\20200528 GOKEGG\\N-go.no1.annot.txt',sep='\t',header = FALSE)
#KEGG
ko_anno <- read.delim('E:\\Nepal_project\\20200528 GOKEGG\\N-kegg.no1.annot.txt',sep='\t',header = FALSE)


#---Bi-----
gene_select <- read.delim('Bi_conserved.txt',sep = '\t',header = FALSE)
gene_select <-as.vector(gene_select$V1)
go_enrich <- enricher(gene = gene_select,
                      TERM2GENE = go_anno[c('V2', 'V1')], 
                      TERM2NAME = go_anno[c('V2', 'V3')], 
                      pvalueCutoff = 0.05, 
                      pAdjustMethod = 'BH')
write.csv(as.data.frame(go_enrich),"Bi_conserved_GO.csv",row.names = F)
barplot(go_enrich,showCategory = 20)+ylab('Gene Number')
dotplot(go_enrich,showCategory=20)


ko_enrich <- enricher(gene = gene_select,
                      TERM2GENE = ko_anno[c('V2', 'V1')], 
                      TERM2NAME = ko_anno[c('V2', 'V3')], 
                      pvalueCutoff = 0.05, 
                      pAdjustMethod = 'BH')

write.csv(as.data.frame(ko_enrich),"Bi_conserved_KO.csv",row.names = F)
barplot(ko_enrich,showCategory = 20)+ylab('Gene Number')
dotplot(ko_enrich,showCategory=20)

#---Tri-----
gene_select <- read.delim('Tri_conserved.txt',sep = '\t',header = FALSE)
gene_select <-as.vector(gene_select$V1)
go_enrich <- enricher(gene = gene_select,
                      TERM2GENE = go_anno[c('V2', 'V1')], 
                      TERM2NAME = go_anno[c('V2', 'V3')], 
                      pvalueCutoff = 0.05, 
                      pAdjustMethod = 'BH')
write.csv(as.data.frame(go_enrich),"Tri_conserved_GO.csv",row.names = F)
barplot(go_enrich,showCategory = 20)+ylab('Gene Number')
dotplot(go_enrich,showCategory=20)


ko_enrich <- enricher(gene = gene_select,
                      TERM2GENE = ko_anno[c('V2', 'V1')], 
                      TERM2NAME = ko_anno[c('V2', 'V3')], 
                      pvalueCutoff = 0.05, 
                      pAdjustMethod = 'BH')

write.csv(as.data.frame(ko_enrich),"Tri_conserved_KO.csv",row.names = F)
dotplot(ko_enrich,showCategory=20)

#---All-----
gene_select <- read.delim('All_conserved.txt',sep = '\t',header = FALSE)
gene_select <-as.vector(gene_select$V1)
go_enrich <- enricher(gene = gene_select,
                      TERM2GENE = go_anno[c('V2', 'V1')], 
                      TERM2NAME = go_anno[c('V2', 'V3')], 
                      pvalueCutoff = 0.05, 
                      pAdjustMethod = 'BH')
write.csv(as.data.frame(go_enrich),"All_conserved_GO.csv",row.names = F)
dotplot(go_enrich,showCategory=20)

ko_enrich <- enricher(gene = gene_select,
                      TERM2GENE = ko_anno[c('V2', 'V1')], 
                      TERM2NAME = ko_anno[c('V2', 'V3')], 
                      pvalueCutoff = 0.05, 
                      pAdjustMethod = 'BH')

write.csv(as.data.frame(ko_enrich),"All_conserved_KO.csv",row.names = F)
dotplot(ko_enrich,showCategory=20)

#---BiTri-----
gene_select <- read.delim('BiTri-conserved.txt',sep = '\t',header = FALSE)
gene_select <-as.vector(gene_select$V1)
go_enrich <- enricher(gene = gene_select,
                      TERM2GENE = go_anno[c('V2', 'V1')], 
                      TERM2NAME = go_anno[c('V2', 'V3')], 
                      pvalueCutoff = 0.05, 
                      pAdjustMethod = 'BH')
write.csv(as.data.frame(go_enrich),"BiTri_conserved_GO.csv",row.names = F)
barplot(go_enrich,showCategory=20)

ko_enrich <- enricher(gene = gene_select,
                      TERM2GENE = ko_anno[c('V2', 'V1')], 
                      TERM2NAME = ko_anno[c('V2', 'V3')], 
                      pvalueCutoff = 0.05, 
                      pAdjustMethod = 'BH')

write.csv(as.data.frame(ko_enrich),"BiTri_conserved_KO.csv",row.names = F)
dotplot(ko_enrich,showCategory=20)


###--------------------------------------------------------------------
##2021.04.21 NpX文章 BiTriAll=conserved 堆积柱状图做显著性分析       ##
##---------------------------------------------------------------------
library(ggplot2)
library(xlsx)
library(reshape2)
library(dplyr)
library(ggsignif)
df <- read.xlsx('E:\\Nepal_project\\Allele-compartment-conserved-non-conserved-count.xlsx',sheetIndex = 2)
df <- melt(df,id.vars = 'chr')

compaired <- list(c("Chr3", "Chr8"))
                  
ggplot(df,aes(chr,fill=variable))+
  geom_bar(aes(y=value),stat = 'identity',position = 'fill',color='black')+
  xlab('')+ylab('Ratio')+
  scale_x_discrete(labels=c('Chr1','Chr2','Chr3','Chr4','Chr5','Chr6','Chr7','Chr8','Chr9','Chr10'))+
  theme_bw()+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank())+theme(legend.title  =  element_blank())+
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,y_position= 2)


###--------------------------------------------------------------------
##2021.04.21 NpX文章 Chr5&8 B2A orthologues RNA-Seq                  ##
##---------------------------------------------------------------------
library(ggplot2)
library(xlsx)
library(reshape2)
library(dplyr)
setwd('E:\\Nepal_project\\20210420_Article_Revision\\NpChr5_8_B2A')
df_npx_chr5_B <- read.xlsx('NpX_Chr5_B_fpkm.xlsx',sheetIndex = 1)
df_npx_chr8_B <- read.xlsx('NpX_Chr8_B_fpkm.xlsx',sheetIndex = 1)
NpXChr5B2A_AP85df <- read.xlsx('NpXChr5B2A_AP85.fpkm.xlsx',sheetIndex = 1)
NpXChr8B2A_AP85df <- read.xlsx('NpXChr8B2A_AP85.fpkm.xlsx',sheetIndex = 1)
df_npx_chr5_B <- df_npx_chr5_B %>% mutate(Chr='Chr5',terms = 'NpX') %>% select(L,S,Chr,terms)
df_npx_chr8_B <- df_npx_chr8_B %>% mutate(Chr='Chr8',terms = 'NpX') %>% select(L,S,Chr,terms)
NpXChr5B2A_AP85df <- NpXChr5B2A_AP85df %>% mutate(Chr='Chr5',terms = 'AP85') %>% select(L,S,Chr,terms)
NpXChr8B2A_AP85df <- NpXChr8B2A_AP85df %>% mutate(Chr='Chr8',terms = 'AP85') %>% select(L,S,Chr,terms)
dd <- rbind(df_npx_chr5_B,df_npx_chr8_B,NpXChr5B2A_AP85df,NpXChr8B2A_AP85df)
#加了errorbar
dd %>% filter(L > 1) %>% ggplot(aes(x=Chr,y=L,fill=terms))+
  geom_boxplot(outlier.colour = NA,width=0.55,position=position_dodge(0.7))+
  stat_boxplot(mapping = aes(x=Chr,y=L,fill=terms),geom = 'errorbar',width=0.2,position=position_dodge(0.7))+
  ylim(0,15)+
  scale_fill_manual(values = c("#B2DF8A","#7FD2FF"))+
  theme_bw()+xlab('')+ylab('FPKM')+theme(legend.position = 'top',legend.title  =  element_blank())


dd %>% filter(S > 1) %>% ggplot(aes(x=Chr,y=S,fill=terms))+
  geom_boxplot(outlier.colour = NA,width=0.55,position=position_dodge(0.7))+
  stat_boxplot(mapping = aes(x=Chr,y=S,fill=terms),geom = 'errorbar',width=0.2,position=position_dodge(0.7))+
  ylim(0,15)+
  scale_fill_manual(values = c("#B2DF8A","#7FD2FF"))+
  theme_bw()+xlab('')+ylab('FPKM')+theme(legend.position = 'top',legend.title  =  element_blank())

#----------------------10chr比较---------------------------
setwd('E:\\Nepal_project\\20210420_Article_Revision\\NpChr5_8_B2A\\10chr比较')
df_npx_chr1 <- read.xlsx('NpX_Chr1_B.fpkm.xlsx',sheetIndex = 1)
df_npx_chr2 <- read.xlsx('NpX_Chr2_B.fpkm.xlsx',sheetIndex = 1)
df_npx_chr3 <- read.xlsx('NpX_Chr3_B.fpkm.xlsx',sheetIndex = 1)
df_npx_chr4 <- read.xlsx('NpX_Chr4_B.fpkm.xlsx',sheetIndex = 1)
df_npx_chr5 <- read.xlsx('NpX_Chr5_B.fpkm.xlsx',sheetIndex = 1)
df_npx_chr6 <- read.xlsx('NpX_Chr6_B.fpkm.xlsx',sheetIndex = 1)
df_npx_chr7 <- read.xlsx('NpX_Chr7_B.fpkm.xlsx',sheetIndex = 1)
df_npx_chr8 <- read.xlsx('NpX_Chr8_B.fpkm.xlsx',sheetIndex = 1)
df_npx_chr9 <- read.xlsx('NpX_Chr9_B.fpkm.xlsx',sheetIndex = 1)
df_npx_chr10 <- read.xlsx('NpX_Chr10_B.fpkm.xlsx',sheetIndex = 1)

df_AP85_chr1 <- read.xlsx('NpX_Chr1_AP85ortho_A.fpkm.xlsx',sheetIndex = 1)
colnames(df_AP85_chr1) <- c('gene_id','L','S1','S','S3')
df_AP85_chr2 <- read.xlsx('NpX_Chr2_AP85ortho_A.fpkm.xlsx',sheetIndex = 1)
colnames(df_AP85_chr2) <- c('gene_id','L','S1','S','S3')
df_AP85_chr3 <- read.xlsx('NpX_Chr3_AP85ortho_A.fpkm.xlsx',sheetIndex = 1)
colnames(df_AP85_chr3) <- c('gene_id','L','S1','S','S3')
df_AP85_chr4 <- read.xlsx('NpX_Chr4_AP85ortho_A.fpkm.xlsx',sheetIndex = 1)
colnames(df_AP85_chr4) <- c('gene_id','L','S1','S','S3')
df_AP85_chr5 <- read.xlsx('NpX_Chr5_AP85ortho_A.fpkm.xlsx',sheetIndex = 1)
colnames(df_AP85_chr5) <- c('gene_id','L','S1','S','S3')
df_AP85_chr6 <- read.xlsx('NpX_Chr6_AP85ortho_A.fpkm.xlsx',sheetIndex = 1)
colnames(df_AP85_chr6) <- c('gene_id','L','S1','S','S3')
df_AP85_chr7 <- read.xlsx('NpX_Chr7_AP85ortho_A.fpkm.xlsx',sheetIndex = 1)
colnames(df_AP85_chr7) <- c('gene_id','L','S1','S','S3')
df_AP85_chr8 <- read.xlsx('NpX_Chr8_AP85ortho_A.fpkm.xlsx',sheetIndex = 1)
colnames(df_AP85_chr8) <- c('gene_id','L','S1','S','S3')
df_AP85_chr9 <- read.xlsx('NpX_Chr9_AP85ortho_A.fpkm.xlsx',sheetIndex = 1)
colnames(df_AP85_chr9) <- c('gene_id','L','S1','S','S3')
df_AP85_chr10 <- read.xlsx('NpX_Chr10_AP85ortho_A.fpkm.xlsx',sheetIndex = 1)
colnames(df_AP85_chr10) <- c('gene_id','L','S1','S','S3')

df_npx_chr1 <- df_npx_chr1 %>% mutate(Chr='Chr1',terms = 'NpX') %>% select(L,S,Chr,terms)
df_npx_chr2 <- df_npx_chr2 %>% mutate(Chr='Chr2',terms = 'NpX') %>% select(L,S,Chr,terms)
df_npx_chr3 <- df_npx_chr3 %>% mutate(Chr='Chr3',terms = 'NpX') %>% select(L,S,Chr,terms)
df_npx_chr4 <- df_npx_chr4 %>% mutate(Chr='Chr4',terms = 'NpX') %>% select(L,S,Chr,terms)
df_npx_chr5 <- df_npx_chr5 %>% mutate(Chr='Chr5',terms = 'NpX') %>% select(L,S,Chr,terms)
df_npx_chr6 <- df_npx_chr6 %>% mutate(Chr='Chr6',terms = 'NpX') %>% select(L,S,Chr,terms)
df_npx_chr7 <- df_npx_chr7 %>% mutate(Chr='Chr7',terms = 'NpX') %>% select(L,S,Chr,terms)
df_npx_chr8 <- df_npx_chr8 %>% mutate(Chr='Chr8',terms = 'NpX') %>% select(L,S,Chr,terms)
df_npx_chr9 <- df_npx_chr9 %>% mutate(Chr='Chr9',terms = 'NpX') %>% select(L,S,Chr,terms)
df_npx_chr10 <- df_npx_chr10 %>% mutate(Chr='Chr10',terms = 'NpX') %>% select(L,S,Chr,terms)


df_AP85_chr1 <- df_AP85_chr1 %>% mutate(Chr='Chr1',terms = 'AP85') %>% select(L,S,Chr,terms)
df_AP85_chr2 <- df_AP85_chr2 %>% mutate(Chr='Chr2',terms = 'AP85') %>% select(L,S,Chr,terms)
df_AP85_chr3 <- df_AP85_chr3 %>% mutate(Chr='Chr3',terms = 'AP85') %>% select(L,S,Chr,terms)
df_AP85_chr4 <- df_AP85_chr4 %>% mutate(Chr='Chr4',terms = 'AP85') %>% select(L,S,Chr,terms)
df_AP85_chr5 <- df_AP85_chr5 %>% mutate(Chr='Chr5',terms = 'AP85') %>% select(L,S,Chr,terms)
df_AP85_chr6 <- df_AP85_chr6 %>% mutate(Chr='Chr6',terms = 'AP85') %>% select(L,S,Chr,terms)
df_AP85_chr7 <- df_AP85_chr7 %>% mutate(Chr='Chr7',terms = 'AP85') %>% select(L,S,Chr,terms)
df_AP85_chr8 <- df_AP85_chr8 %>% mutate(Chr='Chr8',terms = 'AP85') %>% select(L,S,Chr,terms)
df_AP85_chr9 <- df_AP85_chr9 %>% mutate(Chr='Chr9',terms = 'AP85') %>% select(L,S,Chr,terms)
df_AP85_chr10 <- df_AP85_chr10 %>% mutate(Chr='Chr10',terms = 'AP85') %>% select(L,S,Chr,terms)

dd <- rbind(df_npx_chr1,df_npx_chr2,df_npx_chr3,df_npx_chr4,df_npx_chr5,df_npx_chr6,df_npx_chr7,
            df_npx_chr8,df_npx_chr9,df_npx_chr10,df_AP85_chr1,df_AP85_chr2,df_AP85_chr3,df_AP85_chr4,
            df_AP85_chr5,df_AP85_chr6,df_AP85_chr7,df_AP85_chr8,df_AP85_chr9,df_AP85_chr10)

dd %>% filter(L > 0) %>% ggplot(aes(x=Chr,y=L,fill=terms))+
  geom_boxplot(outlier.colour = NA,width=0.55,position=position_dodge(0.7))+
  stat_boxplot(mapping = aes(x=Chr,y=L),geom = 'errorbar',width=0.4,position=position_dodge(0.7))+
  ylim(0,15)+
  scale_x_discrete(labels=c('Chr1','Chr2','Chr3','Chr4','Chr5','Chr6','Chr7','Chr8','Chr9','Chr10'))+
  scale_fill_manual(values = c("#B2DF8A","#7FD2FF"))+
  theme_bw()+xlab('')+ylab('FPKM')+theme(legend.position = 'top',legend.title  =  element_blank())

dd %>% filter(S > 0) %>% ggplot(aes(x=Chr,y=S,fill=terms))+
  geom_boxplot(outlier.colour = NA,width=0.55,position=position_dodge(0.7))+
  stat_boxplot(mapping = aes(x=Chr,y=S),geom = 'errorbar',width=0.4,position=position_dodge(0.7))+
  ylim(0,25)+
  scale_x_discrete(labels=c('Chr1','Chr2','Chr3','Chr4','Chr5','Chr6','Chr7','Chr8','Chr9','Chr10'))+
  scale_fill_manual(values = c("#B2DF8A","#7FD2FF"))+
  theme_bw()+xlab('')+ylab('FPKM')+theme(legend.position = 'top',legend.title  =  element_blank())


###--------------------------------------------------------------------
##2021.05.06 对甘蔗ACR在每条染色体上每个等位进行统计                ##
##---------------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(reshape2)
library(ggprism)

setwd('E:\\潘浩然\\调控元件生信任务\\analysis\\Sugarcane\\peaks')
df <- read.delim('SES208_peak_3category.txt',sep = '\t',header = F)
head(df)

ggplot(df,aes(x=V2,y=V4,fill=V3))+geom_bar(stat='identity')+
  facet_grid(.~V1)+theme_bw()+
  scale_fill_manual(values = c("#FF7F00","#1F78B4","#B2DF8A"))+
  theme_prism(border = TRUE, base_rect_size = 1) +
  coord_cartesian(clip = "off")+xlab('')+ylab('ACR numbers')

#ratio bar
df_chr1A <- df %>% filter(V1=='Chr1',V2=='A') %>% mutate(V5=V4/sum(V4))
df_chr1B <- df %>% filter(V1=='Chr1',V2=='B') %>% mutate(V5=V4/sum(V4))
df_chr1C <- df %>% filter(V1=='Chr1',V2=='C') %>% mutate(V5=V4/sum(V4))
df_chr1D <- df %>% filter(V1=='Chr1',V2=='D') %>% mutate(V5=V4/sum(V4))
df_chr2A <- df %>% filter(V1=='Chr2',V2=='A') %>% mutate(V5=V4/sum(V4))
df_chr2B <- df %>% filter(V1=='Chr2',V2=='B') %>% mutate(V5=V4/sum(V4))
df_chr2C <- df %>% filter(V1=='Chr2',V2=='C') %>% mutate(V5=V4/sum(V4))
df_chr2D <- df %>% filter(V1=='Chr2',V2=='D') %>% mutate(V5=V4/sum(V4))
df_chr3A <- df %>% filter(V1=='Chr3',V2=='A') %>% mutate(V5=V4/sum(V4))
df_chr3B <- df %>% filter(V1=='Chr3',V2=='B') %>% mutate(V5=V4/sum(V4))
df_chr3C <- df %>% filter(V1=='Chr3',V2=='C') %>% mutate(V5=V4/sum(V4))
df_chr3D <- df %>% filter(V1=='Chr3',V2=='D') %>% mutate(V5=V4/sum(V4))
df_chr4A <- df %>% filter(V1=='Chr4',V2=='A') %>% mutate(V5=V4/sum(V4))
df_chr4B <- df %>% filter(V1=='Chr4',V2=='B') %>% mutate(V5=V4/sum(V4))
df_chr4C <- df %>% filter(V1=='Chr4',V2=='C') %>% mutate(V5=V4/sum(V4))
df_chr4D <- df %>% filter(V1=='Chr4',V2=='D') %>% mutate(V5=V4/sum(V4))
df_chr5A <- df %>% filter(V1=='Chr5',V2=='A') %>% mutate(V5=V4/sum(V4))
df_chr5B <- df %>% filter(V1=='Chr5',V2=='B') %>% mutate(V5=V4/sum(V4))
df_chr5C <- df %>% filter(V1=='Chr5',V2=='C') %>% mutate(V5=V4/sum(V4))
df_chr5D <- df %>% filter(V1=='Chr5',V2=='D') %>% mutate(V5=V4/sum(V4))
df_chr6A <- df %>% filter(V1=='Chr6',V2=='A') %>% mutate(V5=V4/sum(V4))
df_chr6B <- df %>% filter(V1=='Chr6',V2=='B') %>% mutate(V5=V4/sum(V4))
df_chr6C <- df %>% filter(V1=='Chr6',V2=='C') %>% mutate(V5=V4/sum(V4))
df_chr6D <- df %>% filter(V1=='Chr6',V2=='D') %>% mutate(V5=V4/sum(V4))
df_chr7A <- df %>% filter(V1=='Chr7',V2=='A') %>% mutate(V5=V4/sum(V4))
df_chr7B <- df %>% filter(V1=='Chr7',V2=='B') %>% mutate(V5=V4/sum(V4))
df_chr7C <- df %>% filter(V1=='Chr7',V2=='C') %>% mutate(V5=V4/sum(V4))
df_chr7D <- df %>% filter(V1=='Chr7',V2=='D') %>% mutate(V5=V4/sum(V4))
df_chr8A <- df %>% filter(V1=='Chr8',V2=='A') %>% mutate(V5=V4/sum(V4))
df_chr8B <- df %>% filter(V1=='Chr8',V2=='B') %>% mutate(V5=V4/sum(V4))
df_chr8C <- df %>% filter(V1=='Chr8',V2=='C') %>% mutate(V5=V4/sum(V4))
df_chr8D <- df %>% filter(V1=='Chr8',V2=='D') %>% mutate(V5=V4/sum(V4))

df <- rbind(df_chr1A,df_chr1B,df_chr1C,df_chr1D,
            df_chr2A,df_chr2B,df_chr2C,df_chr2D,
            df_chr3A,df_chr3B,df_chr3C,df_chr3D,
            df_chr4A,df_chr4B,df_chr4C,df_chr4D,
            df_chr5A,df_chr5B,df_chr5C,df_chr5D,
            df_chr6A,df_chr6B,df_chr6C,df_chr6D,
            df_chr7A,df_chr7B,df_chr7C,df_chr7D,
            df_chr8A,df_chr8B,df_chr8C,df_chr8D)

ggplot(df,aes(x=V2,y=V4,fill=V3))+geom_bar(stat='identity')+
  facet_grid(.~V1)+theme_bw()+
  scale_fill_manual(values = c("#FF7F00","#1F78B4","#B2DF8A"))+
  theme_prism(border = TRUE, base_rect_size = 1) +
  coord_cartesian(clip = "off")+xlab('')+ylab('ACR numbers')

mycolors <- c("#FF9D1E","#007ED3","#B2DF8A")
ggplot(df,aes(x=V2,y=V5,fill=V3))+geom_bar(stat='identity')+
  facet_grid(.~V1)+theme_bw()+scale_fill_manual(values = mycolors)+theme_bw()+
  theme_prism(border = TRUE, base_rect_size = 1) +
  coord_cartesian(clip = "off")+xlab('')+ylab('Proportion of different ACR')+
  theme(legend.key.size=unit(0.8,"cm"),legend.key.width=unit(0.8,"cm"),
        legend.spacing.x=unit(0.3,'cm'),legend.text = element_text(size = 15))

ggsave('E:\\潘浩然\\调控元件生信任务\\analysis\\Sugarcane\\peaks\\SES208_3ACR_ratiobarplot.pdf',
       width=14,height=10)


###------------------------------------------------------------------------
##2021.05.07 对甘蔗ACR在全基因组上7个种类表达量的boxplot并添加显著性标记 ##
##-------------------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(ggprism)
library(openxlsx)
library(ggsignif)
library(forcats)

setwd('E:\\潘浩然\\调控元件生信任务\\analysis\\Sugarcane\\peaks')

#df <- read.delim('SES208_peak_3category.txt',sep = '\t',header = F)
df_onlydACR <- read.xlsx('ACR7category.xlsx',sheet = 1)
df_onlygACR <- read.xlsx('ACR7category.xlsx',sheet = 5)
df_onlypACR <- read.xlsx('ACR7category.xlsx',sheet = 7)
df_gdACR <- read.xlsx('ACR7category.xlsx',sheet = 2)
df_dpACR <- read.xlsx('ACR7category.xlsx',sheet = 3)
df_gdpACR <- read.xlsx('ACR7category.xlsx',sheet = 4)
df_gpACR <- read.xlsx('ACR7category.xlsx',sheet = 6)
allgene <- read.xlsx('ACR7category.xlsx',sheet = 8)
non_ACR <- read.xlsx('ACR7category.xlsx',sheet = 9)

df_onlydACR <- df_onlydACR %>% mutate(type='onlydACR') %>% select(type,FPKM)
df_onlygACR <- df_onlygACR %>% mutate(type='onlygACR') %>% select(type,FPKM)
df_onlypACR <- df_onlypACR %>% mutate(type='onlypACR') %>% select(type,FPKM)
df_gdACR <- df_gdACR %>% mutate(type='gdACR') %>% select(type,FPKM)
df_gpACR <- df_gpACR %>% mutate(type='gpACR') %>% select(type,FPKM)
df_dpACR <- df_dpACR %>% mutate(type='dpACR') %>% select(type,FPKM)
df_gdpACR <- df_gdpACR %>% mutate(type='gdpACR') %>% select(type,FPKM)
df_allgene <- allgene %>% mutate(type='allgene') %>% select(type,FPKM)
df_nonACR <- non_ACR %>% mutate(type='nonACR') %>% select(type,FPKM)



#df <- rbind(df_onlydACR,df_onlygACR,df_onlypACR,df_gdACR,df_dpACR,df_gdpACR,df_gpACR,df_allgene,df_nonACR)
df <- rbind(df_onlydACR,df_onlygACR,df_onlypACR,df_gdpACR,df_allgene,df_nonACR)
head(df_onlydACR)
head(df_onlygACR)


mycolors <- c("#7FD2FF","#B2DF8A","#FFACAA","#FF9D1E","#007ED3","#894FC6")
compaired <- list(c("onlygACR", "onlydACR"), 
                  c("onlydACR","onlypACR"), 
                  c("onlygACR","onlypACR"),
                  c("gdpACR","nonACR"),
                  c("nonACR","allgene"),
                  c("gdpACR","allgene"),
                  c("onlypACR","allgene"),
                  c("onlydACR","allgene"),
                  c("onlygACR","allgene"))


ggplot(df,aes(x=fct_relevel(type,"onlygACR","onlydACR","onlypACR","gdpACR","nonACR","allgene"),y=FPKM,fill=type))+
  geom_boxplot(width=0.5,alpha=0.8,outlier.fill = "black")+scale_fill_manual(values = mycolors)+
  geom_signif(comparisons = compaired,map_signif_level = TRUE,step_increase = 0.1)+
  theme_classic()+xlab("")+ylab('FPKM')+theme(axis.text.x = element_text(hjust = 1,angle=45))


###------------------------------------------------------------------------
##2021.05.07 对含有不同数目ACRs的基因的表达量做boxplot并添加显著性统计   ##
##-------------------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(ggprism)
library(openxlsx)
library(ggsignif)
library(forcats)

setwd('E:\\潘浩然\\调控元件生信任务\\analysis\\Sugarcane\\peaks')



df_1 <- read.xlsx('pergeneACRs_onlygpACR.xlsx',sheet = 1)
df_2 <- read.xlsx('pergeneACRs_onlygpACR.xlsx',sheet = 2)
df_3 <- read.xlsx('pergeneACRs_onlygpACR.xlsx',sheet = 3)
# df_4 <- read.xlsx('pergeneACRs_onlygpACR.xlsx',sheet = 4)
# df_5 <- read.xlsx('pergeneACRs.xlsx',sheet = 5)
# df_6 <- read.xlsx('pergeneACRs.xlsx',sheet = 6)
# df_7 <- read.xlsx('pergeneACRs.xlsx',sheet = 9)
# df_8 <- read.xlsx('pergeneACRs.xlsx',sheet = 8)
# df_9 <- read.xlsx('pergeneACRs.xlsx',sheet = 7)
# df_10<- read.xlsx('pergeneACRs.xlsx',sheet = 10)


df_1 <- df_1 %>% mutate(type = '1ACR')%>% select(type,FPKM)
df_2 <- df_2 %>% mutate(type = '2ACR')%>% select(type,FPKM)
df_3 <- df_3 %>% mutate(type = '3ACR')%>% select(type,FPKM)
# df_4 <- df_4 %>% mutate(type = '4ACR')%>% select(type,FPKM)
# df_5 <- df_5 %>% mutate(type = '5ACR')%>% select(type,FPKM)
# df_6 <- df_6 %>% mutate(type = '6ACR')%>% select(type,FPKM)
# df_7 <- df_7 %>% mutate(type = '7ACR')%>% select(type,FPKM)
# df_8 <- df_8 %>% mutate(type = '8ACR')%>% select(type,FPKM)
# df_9 <- df_9 %>% mutate(type = '9ACR')%>% select(type,FPKM)
# df_10 <- df_10 %>% mutate(type = '>10ACR')%>% select(type,FPKM)

compaired <- list(c("1ACR", "2ACR"), 
                  c("2ACR","3ACR"), 
                  c("1ACR","3ACR"))

df <- rbind(df_1,df_2,df_3)
ggplot(df,aes(x=type,y=FPKM,fill=type))+geom_boxplot(width=0.5,position = 'dodge2')+
  xlab('')+ylab('log10[FPKM]')+
  geom_signif(comparisons = compaired,map_signif_level = TRUE,step_increase = 0.1)+
  theme_classic()+theme_prism(base_rect_size = 1)+theme(legend.position="right")


###------------------------------------------------------------------------
##2021.05.08 对SES208全基因组范围内ACR与最近基因之间距离做密度分布图     ##
##-------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(openxlsx)

setwd('E:\\潘浩然\\调控元件生信任务\\analysis\\Sugarcane\\peaks')
df <- read.xlsx('SES208_peak_category.xlsx')
head(df)
df <- df %>% mutate(distance=abs(distanceToTSS)/1000) %>% select(distance,category)
ggplot(df,aes(distance))+geom_density()



#---------------canu1.9 badila phasing后group长度的统计-------------------
setwd('E:\\Badila')
df <- read.delim('groups.asm.fasta.sizes',header = FALSE)
df <- head(df,80)
df <- df %>% mutate(V3=V2/1000000)
df <- df %>% filter(V1 %in% c('group01','group10','group20','group30','group40','group50','group60','group70','group80'))
ggplot(df,aes(factor(V1),V3,fill=V1))+
  geom_bar(stat = 'identity',width = 0.5)+
  theme_classic()+theme(legend.position='none')+coord_flip()

ggplot(df,aes(x=fct_reorder(V1,V3),V3,fill=V1))+
  geom_bar(stat = 'identity',width = 0.5)+coord_flip()+
  scale_y_continuous(position = 'right',expand = c(0,0))+theme_classic()+
  xlab('')+ylab('Superscaffold length')+
  theme_prism(base_rect_size = 1)+theme(legend.position='none')





###--------------------------------------------------------------------
##2021.06.03 Badila基因组survery##
##---------------------------------------------------------------------
options(scipen = 999)
library(tidyverse)
df <- read.table('E:\\Badila\\survey\\reads.jf.histo',sep = ' ',header = FALSE)
head(df)

ggplot(df,aes(x=V1,y=V2))+geom_line()+xlim(0,2000)+ylim(0,12000000)


#gce 17mer 0603
a <- read.table('E:\\中转站\\freq.stat.2colum',sep = '\t')
ggplot(a,aes(x=V1,y=V2))+geom_line(color="#f46f20",size=0.75)+xlim(0,1000)+ylim(0,6500000)+theme_bw()

#jellyfish 17mer 0602
b <- read.table('E:\\中转站\\17mer_out.histo',sep = ' ')
ggplot(b,aes(x=V1,y=V2))+geom_line(color="#f46f20",size=0.75)+xlim(0,1500)+ylim(0,6500000)+theme_bw()

#jellyfish 21mer 0603
c <- read.table('E:\\中转站\\21mer_out.histo',sep = ' ')
head(c)
c <- c %>% mutate(V3=V2/1000000)
ggplot(c,aes(x=V1,y=V2))+geom_line(color="#007ED3",size=0.75)+
  xlim(0,1500)+ylim(0,13000000)+theme_bw()+
  xlab('Coverage(X)')+ylab('Counts')


cc <- read.table('E:\\中转站\\21zcat_outs2.histo',sep = ' ')
ggplot(cc,aes(x=V1,y=V2))+geom_line(color="#007ED3",size=0.75)+
  xlim(0,500)+ylim(0,20000000)+theme_bw()+
  xlab('Coverage(X)')+ylab('Counts')


#jellyfish 21mer 1004
c <- read.table('E:\\Badila\\survey\\jellyfish\\ZG.21mer.histo',sep = ' ')
head(c)
ggplot(c,aes(x=V1,y=V2))+geom_line(color="#007ED3",size=1)+
  xlim(0,1000)+ylim(0,15000000)+theme_bw()+
  xlab('Coverage(X)')+ylab('Counts')

ggsave('ZGsurvey.21mer.pdf',dpi=300,width = 14,height = 12)


#gce 17mer 1004
setwd('E:\\Badila\\survey\\gce')
df <- read.table('17mer.kmer.freq.stat',comment.char = '#')
head(df)
colnames(df) <- c('Kmer_Frequency','Kmer_Species_Number','Kmer_Species_Ratio',
                  'Kmer_Species_accumulate_Ratio','Kmer_Individual_Number',
                  'Kmer_Individual_Ratio','Kmer_Individual_accumulate_ratio')

head(df)
ggplot(df,aes(x=Kmer_Frequency,y=Kmer_Species_Ratio))+geom_line()+
  xlim(0,250)

###--------------------------------------------------------------------
##2021.06.05 协助秀婷师姐画象限图##
##---------------------------------------------------------------------
library(tidyverse)
library(openxlsx)
library(ggprism)

df <- read.xlsx('E:\\中转站\\Xiuting_allG_vs_allP_v0605.matchType.xlsx',sheet=1)
mycol <- c("#006b7b","#f88421","grey","green","#ef1828")
ggplot(df,aes(x=Plog2FC,y=Glog2FC))+geom_point(aes(color=Type),size=1.5)+
  geom_hline(yintercept = 1,size=1,linetype='dashed')+
  geom_hline(yintercept = -1,size=1,linetype='dashed') + 
  geom_vline(xintercept = 0.584,size=1,linetype='dashed')+
  geom_vline(xintercept = -0.584,size=1,linetype='dashed')+theme_bw()+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank())+
  theme(legend.title  =  element_blank())+
  scale_color_manual(values = mycol)+theme_prism(border = TRUE, base_rect_size = 1) +
  coord_cartesian(clip = "off")+scale_x_continuous( breaks=seq(-7.0,7.0,1.0))+
  scale_y_continuous(breaks = seq(-15,15,5))+
  xlab('log2 ratio of protein')+ylab('log2 ratio of transcript')


df2 <- read.xlsx('E:\\中转站\\Xiuting_allG_vs_allP_v0605.matchType.xlsx',sheet=2)
ggplot(df2,aes(x=Plog2FC,y=Glog2FC))+geom_point(aes(color=Type),size=1.5)+
  geom_hline(yintercept = 1,size=1,linetype='dashed')+
  geom_hline(yintercept = -1,size=1,linetype='dashed') + 
  geom_vline(xintercept = 0.584,size=1,linetype='dashed')+
  geom_vline(xintercept = -0.584,size=1,linetype='dashed')+theme_bw()+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank())+
  theme(legend.title  =  element_blank())+
  scale_color_manual(values = mycol)+theme_prism(border = TRUE, base_rect_size = 1) +
  coord_cartesian(clip = "off")+scale_x_continuous( breaks=seq(-7.0,7.0,1.0))+
  scale_y_continuous(breaks = seq(-15,10,5))+
  xlab('log2 ratio of protein')+ylab('log2 ratio of transcript')


###--------------------------------------------------------------------
##2021.07.11 Subject-binFPKM with ACRs##
##---------------------------------------------------------------------
library(openxlsx)
library(tidyverse)
library(magrittr)
library(ggsci)
library(patchwork)
#pal_mk <- wes_palette("Darjeeling1")
setwd('E:\\潘浩然\\调控元件生信任务\\analysis\\Sugarcane\\peaks\\binFPKMwithACR')
### 基因中含有ACR的比例
dfall <- read.xlsx('binFPKM2ACRratio.xlsx',sheet = 1,colNames = FALSE)
dfg <- read.xlsx('binFPKM2ACRratio.xlsx',sheet = 2,colNames = FALSE)
dfd <- read.xlsx('binFPKM2ACRratio.xlsx',sheet = 3,colNames = FALSE)
dfp <- read.xlsx('binFPKM2ACRratio.xlsx',sheet = 4,colNames = FALSE)

dfall %<>% mutate(type='allACR',binID=paste(X1,'all',sep = ''))
dfg %<>% mutate(type='gACR',binID=paste(X1,'g',sep = ''))
dfd %<>% mutate(type='dACR',binID=paste(X1,'d',sep = ''))
dfp %<>% mutate(type='pACR',binID=paste(X1,'p',sep = ''))

df <- rbind(dfall,dfg,dfd,dfp)
head(df)
#a <- ggplot(df,aes(x=type,y=X2,fill=X1))+geom_bar(stat='identity',position = 'dodge')
a <- ggplot(df,aes(x=fct_relevel(binID,'bin1all','bin2all','bin3all','bin4all','bin5all',
                            'bin1d','bin2d','bin3d','bin4d','bin5d','bin1p','bin2p',
                            'bin3p','bin4p','bin5p','bin1g','bin2g','bin3g','bin4g','bin5g'),
              y=X2,fill=type))+geom_bar(stat='identity',position = 'dodge')+ylim(0,40)+
  theme_bw()+scale_x_discrete(expand = c(0,.9))+scale_fill_npg()+
  xlab('')+ylab('Percentage of genes associated with ACRs(%)')+
  theme(legend.title  =  element_blank(),
        axis.title = element_text(size = 15, color = "black",face = 'bold'),
        axis.text = element_text(size = 15,color = "black"),
        plot.margin = unit(c(1,1,1,1), "cm"),legend.key.size=unit(0.8,"cm"),
        legend.key.width=unit(0.8,"cm"),legend.text = element_text(size = 12,face = 'bold'),
        legend.position = c(0.8, 0.85),panel.background = element_blank(),
        legend.background = element_rect(color = "black",size=0.5),
        legend.margin = margin(t = 5, l = 5, r = 5, b = 5))

### 基因中含有ACR的平均个数
DFall <- read.xlsx('binFPKM2ACRnums.xlsx',sheet = 1,colNames = FALSE)
DFg <- read.xlsx('binFPKM2ACRnums.xlsx',sheet = 2,colNames = FALSE)
DFd <- read.xlsx('binFPKM2ACRnums.xlsx',sheet = 3,colNames = FALSE)
DFp <- read.xlsx('binFPKM2ACRnums.xlsx',sheet = 4,colNames = FALSE)

DFall %<>% mutate(type='allACR',binID=paste(X1,'all',sep = ''))
DFg %<>% mutate(type='gACR',binID=paste(X1,'g',sep = ''))
DFd %<>% mutate(type='dACR',binID=paste(X1,'d',sep = ''))
DFp %<>% mutate(type='pACR',binID=paste(X1,'p',sep = ''))

DF <- rbind(DFall,DFg,DFd,DFp)
head(DF)
# b <- ggplot(DF,aes(x=type,y=X2,group=X1))+geom_line()
# b <- ggplot(DF,aes(x=X1,y=X2,group=1))+geom_line()+facet_wrap(~type)
# b <- ggplot(DF,aes(x=fct_relevel(binID,'bin1all','bin2all','bin3all','bin4all','bin5all',
#                             'bin1d','bin2d','bin3d','bin4d','bin5d','bin1p','bin2p',
#                             'bin3p','bin4p','bin5p','bin1g','bin2g','bin3g','bin4g','bin5g'),
#               y=X2,group=1))+geom_line()
b <- ggplot(DF,aes(x=fct_relevel(binID,'bin1all','bin2all','bin3all','bin4all','bin5all',
                                 'bin1d','bin2d','bin3d','bin4d','bin5d','bin1p','bin2p',
                                 'bin3p','bin4p','bin5p','bin1g','bin2g','bin3g','bin4g','bin5g'),
                   y=X2,group=type,))+geom_line(size=0.7)+geom_point(aes(shape=type))+ylim(0.5,2.5)+
  theme_bw()+ylab('Average ACRs number per gene')+xlab('')+
  theme(legend.title  =  element_blank(),
                   axis.title = element_text(size = 15, color = "black",face = 'bold'),
                   axis.text = element_text(size = 15,color = "black"),
                   plot.margin = unit(c(1,1,1,1), "cm"),legend.key.size=unit(0.8,"cm"),
                   legend.key.width=unit(0.8,"cm"),legend.text = element_text(size = 12,face = 'bold'),
                   legend.position = c(0.8, 0.85),panel.background = element_blank(),
                   legend.background = element_rect(color = "black",size=0.5),
                   legend.margin = margin(t = 5, l = 5, r = 5, b = 5))




library(gtable)
library(grid)

y2_plot <- function(p1, p2) {
  p1 <- ggplotGrob(p1)
  p2 <- ggplotGrob(p2)
  
  # Get the location of the plot panel in p1.
  # These are used later when transformed elements of p2 are put back into p1
  pp <- c(subset(p1$layout, name == 'panel', se = t:r))
  
  # Overlap panel for second plot on that of the first plot
  p1 <- gtable_add_grob(p1, p2$grobs[[which(p2$layout$name == 'panel')]], pp$t, pp$l, pp$b, pp$l)
  
  # Then proceed as before:
  
  # ggplot contains many labels that are themselves complex grob; 
  # usually a text grob surrounded by margins.
  # When moving the grobs from, say, the left to the right of a plot,
  # Make sure the margins and the justifications are swapped around.
  # The function below does the swapping.
  # Taken from the cowplot package:
  # https://github.com/wilkelab/cowplot/blob/master/R/switch_axis.R 
  
  hinvert_title_grob <- function(grob){
    
    # Swap the widths
    widths <- grob$widths
    grob$widths[1] <- widths[3]
    grob$widths[3] <- widths[1]
    grob$vp[[1]]$layout$widths[1] <- widths[3]
    grob$vp[[1]]$layout$widths[3] <- widths[1]
    
    # Fix the justification
    grob$children[[1]]$hjust <- 1 - grob$children[[1]]$hjust 
    grob$children[[1]]$vjust <- 1 - grob$children[[1]]$vjust 
    grob$children[[1]]$x <- unit(1, 'npc') - grob$children[[1]]$x
    grob
  }
  
  # Get the y axis title from p2
  index <- which(p2$layout$name == 'ylab-l') # Which grob contains the y axis title?
  ylab <- p2$grobs[[index]]                # Extract that grob
  ylab <- hinvert_title_grob(ylab)         # Swap margins and fix justifications
  
  # Put the transformed label on the right side of p1
  p1 <- gtable_add_cols(p1, p2$widths[p2$layout[index, ]$l], pp$r)
  p1 <- gtable_add_grob(p1, ylab, pp$t, pp$r + 1, pp$b, pp$r + 1, clip = 'off', name = 'ylab-r')
  
  # Get the y axis from p2 (axis line, tick marks, and tick mark labels)
  index <- which(p2$layout$name == 'axis-l')  # Which grob
  yaxis <- p2$grobs[[index]]                  # Extract the grob
  
  # yaxis is a complex of grobs containing the axis line, the tick marks, and the tick mark labels.
  # The relevant grobs are contained in axis$children:
  #   axis$children[[1]] contains the axis line;
  #   axis$children[[2]] contains the tick marks and tick mark labels.
  
  # First, move the axis line to the left
  yaxis$children[[1]]$x <- unit.c(unit(0, 'npc'), unit(0, 'npc'))
  
  # Second, swap tick marks and tick mark labels
  ticks <- yaxis$children[[2]]
  ticks$widths <- rev(ticks$widths)
  ticks$grobs <- rev(ticks$grobs)
  
  # Third, move the tick marks
  ticks$grobs[[1]]$x <- ticks$grobs[[1]]$x - unit(1, 'npc') + unit(3, 'pt')
  
  # Fourth, swap margins and fix justifications for the tick mark labels
  ticks$grobs[[2]] <- hinvert_title_grob(ticks$grobs[[2]])
  
  # Fifth, put ticks back into yaxis
  yaxis$children[[2]] <- ticks
  
  # Put the transformed yaxis on the right side of p1
  p1 <- gtable_add_cols(p1, p2$widths[p2$layout[index, ]$l], pp$r)
  p1 <- gtable_add_grob(p1, yaxis, pp$t, pp$r + 1, pp$b, pp$r + 1, clip = 'off', name = 'axis-r')
  grid.newpage()
  grid.draw(p1)
}

c <- y2_plot(b,a)
ggsave('sample_merge.pdf',width = 12,height = 6)
a+b
ggsave('binFPKMwithACRs.pdf',width = 18,height = 8)



###--------------------------------------------------------------------
##2021.07.16 Subject-Allele diff fpkm heatmap##
##---------------------------------------------------------------------
library(tidyverse)
library(ggtree)
library(aplot)
library(openxlsx)
library(edgeR)


### DEG两两等位之间
setwd('E:\\潘浩然\\调控元件生信任务\\analysis\\Sugarcane\\peaks\\allele')
allele4 = read.xlsx('4allele.readscount.xlsx')
colnames(allele4) <- c('geneID','fpkm1','fpkm2','fpkm3','fpkm4')
rownames(allele4) <- allele4[,1]
allele4 <- allele4[,-1]
head(allele4)
allele4 <- allele4[rowSums(allele4)>1,]

allele_num <- 'allele4' #TODO change name
query_df <- allele4 %>% select(fpkm3,fpkm4) #TODO change col
group <- 1:2

y <- DGEList(counts=query_df, group = group)
keep <- rowSums(cpm(y)>1) >= 1
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y_bcv <- y
bcv <- 0.1
et <- exactTest(y_bcv, dispersion = bcv ^ 2)
gene1 <- decideTestsDGE(et, p.value = 0.05, lfc = 1)
summary(gene1)

colnames(gene1) <- "Signifi"
results <- cbind(y$counts,et$table,gene1)
DEG <- results[results$Signifi!=0,]
up <- results[results$Signifi==1,]
down <- results[results$Signifi==-1,]
result_name = paste(allele_num,colnames(query_df)[1],colnames(query_df)[2],sep = '_')
sheets <- list('DEG_overview_result'=DEG,'up_regulated_gene'=up,'down_regulated_gene'=down)
write.xlsx(sheets,paste(result_name,'.xlsx',sep = ''),row.names = T)

##画热图
allele4df2heat <- read.xlsx('E:\\潘浩然\\调控元件生信任务\\analysis\\Sugarcane\\peaks\\allele\\DEG\\allele4jvenn_60DEG.xlsx')
colnames(allele4df2heat) <- c('geneID','fpkm1','fpkm2','fpkm3','fpkm4')
rownames(allele4df2heat) <- allele4df2heat[,1]
allele4df2heat <- allele4df2heat[,-1]
# scale(allele4df2heat) %>% head()
a <- allele4df2heat %>%  mutate_if(is.numeric,function(x) x+ 1) %>% log10()
rownames(a) <- rownames(allele4df2heat)
a %>% mutate(ID=rownames(a))
aa <- a %>% mutate(ID=rownames(a)) %>% gather(.,1:4,key="condition", value='expr')

phr <- hclust(dist(a)) %>% ggtree(layout="rectangular", branch.length="none")
phc <- hclust(dist(t(a))) %>% ggtree() + layout_dendrogram()
P <- ggplot(aa,aes(x=fct_relevel(condition,'fpkm1','fpkm2','fpkm3','fpkm4'),y=ID,fill=expr)) + geom_tile()+
  theme_minimal()+
  scale_fill_viridis_c() +
  scale_y_discrete(position="right")+
  xlab(NULL) + ylab(NULL)
PP <- P %>% insert_left(phr,width = .2) 
ggsave('allele4_60DEG.fpkm.heatmap.pdf',PP,width = 14,height = 16)
# P %>% insert_left(phr, width=.1) %>%
#   insert_top(phc, height=.1) 


set.seed(56489)                                    
df <- data.frame(x = rep(LETTERS[1:20], each =20),
                 y = letters[1:20],
                 value= runif(50,0,20)) %>% 
  pivot_wider(names_from=x,values_from=value) %>% 
  column_to_rownames("y")

heatmap <- df %>% mutate_if(is.numeric,function(x) x+ 1) %>%
  log10() %>% 
  rownames_to_column("id") %>% 
  pivot_longer(-id) %>% 
  ggplot(.,aes(x=name,y=id,color=value))+
  geom_tile(color="grey70",fill="white",size=0.6)+
  geom_point(shape=19)+
  geom_point(aes(size=abs(value)),show.legend = F)+
  labs(x = NULL,y = NULL,color=NULL) + 
  scale_color_viridis_c()+
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0),position="right")+
  theme(text=element_text(family="Roboto"),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks.x = element_blank(),
        axis.ticks.y=element_blank(),
        panel.border = element_rect(fill=NA,color="grey80",
                                    size=1, linetype="solid"))+
  scale_size(range=c(1,5),guide=NULL)+
  guides(color=guide_colorbar(direction = "vertical",
                              reverse = F,barwidth = unit(.6, "cm"),
                              barheight = unit(10,"cm")))


###--------------------------------------------
###2021.08.05 NpX-v0630 call compartment 500kb
###--------------------------------------------
library(HiTC)
library(tidyverse)
library(wesanderson)

setwd('D:\\FAFU-CGB\\2021summervacation\\NpX\\20210802_matrix_bed\\AB\\500kb')
wrkdir='.'
outdir='compartment'
pal_mk <- wes_palette("Darjeeling1")
data_in=importC(paste(wrkdir,"NpX_500000_iced.modified.removeCtg.matirx",sep="\\"),
                paste(wrkdir,"NpX_500000_abs.changedID.bed",sep="\\"))
fai <- read.table('NpX.v0630.changedID.size.txt',sep = '\t')

chrID <- 'Chr10A'
hic_chr <- extractRegion(data_in$Chr10AChr10A,chr=chrID,from=0,
                         to=fai %>% filter(V1==chrID) %>% .[1,2])
pc1 <- pca.hic(hic_chr, normPerExpected=TRUE, method="loess", npc=1)
write.csv(pc1$PC1,paste(outdir,paste('NpX_v0630',chrID,'compartment.csv',sep = '.'),
                        sep = '\\'),row.names = FALSE)

pc1_df <- as.data.frame(pc1)
pc1_df <- mutate(pc1_df,position=ifelse(start<10,start,start/1000000))
pc1_df[1,9] = 0
p1 <- ggplot(pc1_df,aes(x=position,y=score,fill = ifelse(score>0,'blue','yellow')))+
  geom_bar(stat = 'identity',width=0.5)+
  xlab(chrID)+ylab('')+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank())+
  theme_classic()+
  theme(legend.position="none")+
  #scale_fill_manual(values=c("#6cc396","#9896ca"))
  #scale_fill_manual(values = c("#f46f20","#156077"))
  scale_fill_manual(values = c(pal_mk[4],pal_mk[2]))
p1
ggsave(p1,paste(outdir,paste('NpX_v0630',chrID,'compartment.pdf',sep = '.'),
             sep = '\\'))
ggsave(p1,paste(outdir,paste('NpX_v0630',chrID,'compartment.png',sep = '.'),
             sep = '\\'),dpi=300)

###------------------200kb-------------------------
library(HiTC)
library(tidyverse)
library(wesanderson)

setwd('D:\\FAFU-CGB\\2021summervacation\\NpX\\20210802_matrix_bed\\AB\\200kb')
wrkdir='.'
outdir='compartment'
pal_mk <- wes_palette("Darjeeling1")
data_in=importC(paste(wrkdir,"NpX_200000_iced.modified.noctg.matrix",sep="\\"),
                paste(wrkdir,"NpX_200000_abs.changeID.bed",sep="\\"))
fai <- read.table('D:\\FAFU-CGB\\2021summervacation\\NpX\\20210802_matrix_bed\\NpX.v0630.changedID.size.txt',sep = '\t')

chrID <- 'Chr1A'
hic_chr <- extractRegion(data_in$Chr1AChr1A,chr=chrID,from=0,
                         to=fai %>% filter(V1==chrID) %>% .[1,2])
pc1 <- pca.hic(hic_chr, normPerExpected=TRUE, method="loess", npc=1)


pc1_df <- as.data.frame(pc1)
pc1_df <- mutate(pc1_df,position=ifelse(start<10,start,start/1000000))
pc1_df[1,9] = 0
p1 <- ggplot(pc1_df,aes(x=position,y=score,fill = ifelse(score>0,'blue','yellow')))+
  geom_bar(stat = 'identity',width=0.2)+
  xlab(chrID)+ylab('')+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank())+
  theme_classic()+
  theme(legend.position="none")+
  #scale_fill_manual(values=c("#6cc396","#9896ca"))
  #scale_fill_manual(values = c("#f46f20","#156077"))
  scale_fill_manual(values = c(pal_mk[4],pal_mk[2]))
p1
ggsave(p1,paste(outdir,paste('NpX_v0630',chrID,'compartment.pdf',sep = '.'),
             sep = '\\'))
ggsave(p1,paste(outdir,paste('NpX_v0630',chrID,'compartment.png',sep = '.'),
             sep = '\\'),dpi=300)
write.csv(pc1$PC1,paste(outdir,paste('NpX_v0630',chrID,'compartment.csv',sep = '.'),
                        sep = '\\'),row.names = FALSE)



###--------------------------------------------
###2021.08.09 NpX-v0630 chr-ks boxplot(Fig2f)
###--------------------------------------------
library(tidyverse)
library(tidyquant)
library(readxl)
library(ggsignif)

x1 <- palette_light()
x2 <- palette_dark()

setwd('D:\\FAFU-CGB\\Nepal_project\\20200508 Nepal_chr_Ks\\NpX.v0804.allele_among_KaKs')
df <- read.table('NpX.v0804.allele_among.Ks.tab',header=TRUE)
colnames(df) <- c('Chr','type','Ks')
head(df)
rorder <- paste('chr',1:10,sep = '')
dd <- df %>% filter(Ks>0)
compaired <-list(c("chr5",'chr8'),c("chr5",'chr1'))

#original
p1 <- ggplot(dd,aes(x=fct_relevel(Chr,rorder),y=Ks,color=Chr))+
  geom_boxplot(outlier.colour = NA,width=0.7,size=1)+ylim(0,0.15)+
  scale_color_manual(values = matrix(x2)[,1])+
  stat_boxplot(geom = 'errorbar',width=0.3)+theme_classic()+
  theme(legend.position="none",panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title = element_text(size = 15, color = "black"),
        axis.text.y = element_text(size = 14,face = 'bold'),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
        panel.grid.major=element_blank(),
        axis.text.x = element_text(size = 14,face = 'bold'))+
  xlab('')+ylab('Synonymous substitution rate (Ks)')
  
p1
ggsave('NpX.v0804.allele_among.Ks.boxplot.pdf',p1,width = 16,height = 10,dpi=300)

#coord_flip
p2 <- ggplot(dd,aes(x=fct_relevel(Chr,rorder),y=Ks,color=Chr))+
  geom_boxplot(outlier.colour = NA,width=0.7,size=1)+ylim(0,0.2)+
  scale_color_manual(values = matrix(x2)[,1])+
  stat_boxplot(geom = 'errorbar',width=0.3)+theme_classic()+
  theme(legend.position="none",panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title = element_text(size = 15, color = "black"),
        axis.text.y = element_text(size = 12,face = 'bold'),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
        panel.grid.major=element_blank(),
        axis.text.x = element_text(size = 12,face = 'bold'))+
  xlab('')+ylab('Synonymous substitution rate (Ks)')+
  geom_hline(yintercept=0.05,linetype=2)+coord_flip()
p2
ggsave('NpX.v0804.allele_among.Ks.boxplot.coordflip.pdf',p2,width = 12,height = 15,dpi=300)


###------------------------------------------------
###2021.08.17
###1.NpX2AP85 ABSwitch barplot(Fig2d)
###2.NpX allele BiTriAll-conserved barplot (FigS15)
###3.NpX in four speces ABswitch barplot (FgS18)
###------------------------------------------------
library(tidyverse)
library(readxl)
library(tidyquant)
library(RColorBrewer)
colors<-brewer.pal(name="Set3",9)
colors

setwd('E:\\Nepal_project\\2021summervacation_NpX\\20210812_Nepal-AP85-compartment-variation')
x1 <- palette_light()
x2 <- palette_dark()

##1
df <- read_excel('NPX2AP85.AB_switch.BiTriAll.xlsx',sheet = 1)
qdf <- df %>% select(ChrID,A2B_ratio,B2A_ratio,Con_ratio)
qdc <- qdf %>% pivot_longer(cols = c(A2B_ratio,B2A_ratio,Con_ratio),names_to = 'type',values_to = 'value')

new_o <- paste('Chr',1:10,sep = '')

#mycolor <- c(matrix(x1)[,1][2],matrix(x2)[,1][9],matrix(x2)[,1][6])
#mycolor <- c('#3cb346','#eeb401',"#2e409a")
mycolor <- c(colors[7],colors[6],colors[5])
ggplot(qdc,aes(x=fct_relevel(ChrID,new_o),y=value,fill=type))+
  geom_col(width = 0.9)+theme_bw()+
  theme(axis.line = element_line(color = "black"),
        legend.position = 'top',
        axis.title = element_text(size = 40, color = "black",face = "bold"),
        axis.text.y = element_text(size = 18,color = "black",face = 'bold'),
        axis.text.x = element_text(size = 22,color = "black",face = 'bold'),
        legend.key.size=unit(0.8,"cm"),
        legend.text = element_text(size = 15),
        legend.spacing.x=unit(0.3,'cm'))+
  labs(fill='')+
  scale_fill_manual(values = mycolor)+
  #scale_fill_manual(values = )+
  ylab('Ratio of ortholog regions')+xlab('')
ggsave('NpX2AP85.ABswitch.barplot.pdf',height = 10,width = 12,dpi = 300)


##2.
df2 <- read_excel('NPX2AP85.AB_switch.BiTriAll.xlsx',sheet = 2)
qdf2 <- df2 %>% select(ChrID,Bi,Tri,All)
qdc2 <- pivot_longer(qdf2,cols = c(Bi:All),
                     names_to = 'type',values_to = 'value')
new_o <- paste('Chr',1:10,sep = '')
ggplot(qdc2,aes(x=fct_relevel(ChrID,new_o),y=value,
                fill=fct_relevel(type,'Bi','Tri','All')))+
  geom_col(width = 0.9)+theme_bw()+
  theme(axis.line = element_line(color = "black"),
        legend.position = 'top',
        axis.title = element_text(size = 40, color = "black",face = "bold"),
        axis.text.y = element_text(size = 18,color = "black",face = 'bold'),
        axis.text.x = element_text(size = 22,color = "black",face = 'bold'),
        legend.key.size=unit(1.5,"cm"),
        legend.text = element_text(size = 15),
        legend.spacing.x=unit(0.3,'cm'))+
  labs(fill='')+
  scale_fill_manual(values = mycolor)+
  #scale_fill_manual(values = )+
  ylab('Ratio of compartment')+xlab('')
ggsave('NpXallel.BiTriAllcomserved.flipped.barplot.pdf',height = 10,width = 12,dpi = 300)


##3.
df3 <- read_excel('NPX2AP85.AB_switch.BiTriAll.xlsx',sheet = 3,skip = 1)
head(df3)
qdf3 <- df3 %>% pivot_longer(cols = c(Chr6_A2B:Chr8_Conserved),
                     names_to = 'type',values_to = 'value')
qdf3 %>% head()
qdc3 <- qdf3 %>% mutate(chr=str_sub(type,start = 1,end=4),
                        status=str_sub(type,start = 6)) %>% 
  select(species,chr,status,value)
write.csv(qdc3,'tmp.csv')
df3 <- read_excel('NPX2AP85.AB_switch.BiTriAll.xlsx',sheet = 4)
ggplot(df3,aes(x=fct_relevel(species,'Os vs Sb','Sb vs NpX','NpX vs AP85-441'),
               y=ratio,fill=status))+
  geom_col()+facet_grid(~chr)+
  theme_bw()+
  theme(panel.border = element_rect(color = 'transparent'))+
  theme(axis.text.x=element_text(hjust=1,angle=45))+
  theme(legend.position = 'top',
        axis.title = element_text(size = 40, color = "black",face = "bold"),
        axis.text.y = element_text(size = 18,color = "black",face = 'bold'),
        legend.key.size=unit(1.5,"cm"),
        legend.text = element_text(size = 15),
        legend.spacing.x=unit(0.3,'cm'))+
  labs(fill='')+
  scale_fill_manual(values = mycolor)+xlab('')+ylab('Ratio')
  
ggsave('Fourspecies_ABswitch.barplot.pdf',height = 10,width = 12,dpi = 300)


###------------------------------------------------
###2021.08.20
###NpX2AP85 AB length ratio barplot(FigS14)
###------------------------------------------------
setwd('E:\\Nepal_project\\2021summervacation_NpX\\20210802_matrix_bed\\AB\\500kb\\compartment_ex')
library(tidyverse)
library(readxl)
library(wesanderson)

pal_mk <- wes_palette("Darjeeling1")
df <- read.table('NpX_v0630.allchr.flipped.compartment.bed',sep = '\t',col.names = T)
dd <- df %>% mutate(AB=ifelse(score>0,'A','B')) %>% 
  group_by(seqnames,AB) %>% 
  summarise(n = n(),t_w = sum(width)) %>% 
  group_by(seqnames) %>% summarise(ratio=t_w/sum(t_w))
write.table(dd,'NpX_v0630.allchr.compartment.LengthRatio.bed',sep = '\t')
write.csv(dd,'NpX_v0630.allchr.compartment.LengthRatio.csv',row.names = FALSE)

#s手动在excel中处理后
df <- read.csv('NpX_v0630.allchr.compartment.LengthRatio.csv')
new_o <- paste('Chr',1:10,sep = '')
ggplot(df,aes(allele,ratio,fill=type))+
  geom_col()+theme_bw()+
  facet_grid(~as.factor(seqnames))+
  scale_fill_manual(values = c(pal_mk[4],pal_mk[2]))+
  theme(panel.border = element_rect(color = 'transparent'))+
  theme(panel.spacing = unit(0,'lines'),
        strip.background = element_blank(),strip.placement = 'outside')+
  xlab('')+ylab('Ratio')+
  theme(legend.position="top")+labs(fill = '')
  # theme(axis.text = element_text(size = 15,color = "black",face='bold'),
  #       axis.title = element_text(size = 35, color = "black",face = "bold"),
  #       legend.key.size=unit(1.0,"cm"),legend.spacing.x=unit(0.3,'cm'),
  #       legend.text = element_text(size = 15))

ggsave('NpX.v0804.ChrABratio.baplot1.pdf',dpi=300,width = 16,height = 10)




###------------------------------------------------
###2021.08.24
###NpX 根据王茂军老师意见修改
###1.review1
###2.review2
###------------------------------------------------
library(tidyverse)
library(wesanderson)
library(tidyquant)

#1.review1
x1 <- palette_dark()
pal_mk <- wes_palette("Darjeeling1")
# colors<-brewer.pal(name="Set3",9)
setwd('D:\\FAFU-CGB\\2021summervacation\\NpX\\20210824_WMJreview\\review1')
ABdf <- read.table('NpX_v0630.allchr.compartment.bed',header = TRUE)
head(ABdf)
ABdf <- mutate(ABdf,position=ifelse(start<10,start,start/1000000))
# p1 <- ABdf %>% filter(seqnames=='Chr10A') %>% ggplot(aes(x=position,y=score,fill = ifelse(score>0,'blue','yellow')))+
#   geom_bar(stat = 'identity',width=0.3)+theme(legend.position="none")+
#   scale_fill_manual(values = c(pal_mk[4],pal_mk[2]))+
#   xlab('')+ylab('PC1 score')


numdf <- read.table('NpX.v0630.genome-wide_num.bed',
                    col.names = c('seqnames','start','end','score'))
numdf <- numdf %>% mutate(position=ifelse(start<10,start,start/1000000))
  
# p2 <- ggplot(numdf %>% filter(seqnames=='Chr10A'),aes(x=position,y=score,fill=seqnames))+
#   geom_col()+scale_fill_manual(values =matrix(x1)[,1][2])+
#   theme(legend.position="none")+xlab('')+ylab('Gene num')


fpkmdf <- read.table('NpX.v0630.genome-wide_fpkm.bed',
                     col.names = c('seqnames','start','end','L','S'))

fpkmdf <- fpkmdf %>% mutate(position=ifelse(start<10,start,start/1000000)) 
  
head(fpkmdf)
Ldf <- fpkmdf %>% select(seqnames:L,position) %>% mutate(L=log2(L+1))
# p3 <- ggplot(Ldf %>% filter(seqnames=='Chr10A'),aes(position,L,fill=seqnames))+geom_col()+
#   scale_fill_manual(values = matrix(x1)[,1][6])+
#   theme(legend.position = 'none')+xlab('')+ylab('L fpkm')

Sdf <- fpkmdf %>% select(seqnames:end,S,position) %>% mutate(S=log2(S+1))
# p4 <- ggplot(Sdf %>% filter(seqnames=='Chr10A'),aes(position,S,fill=seqnames))+geom_col()+
#   scale_fill_manual(values = matrix(x1)[,1][1])+
#   theme(legend.position = 'none')+xlab('')+ylab('S fpkm')

plot_grid(p1,p2,p3,p4,labels = c('a','b','c','d'),ncol = 1,nrow = 4)
ggsave('Chr10A_4tracks.pdf',dpi = 300,width = 20,height = 16)

for (i in 1:10) {
  for (j in c('A','B','C','D')) {
    ChrID <- paste('Chr',i,j,sep = '')
    print(ChrID)
    p1 <- ABdf %>% filter(seqnames==ChrID) %>% ggplot(aes(x=position,y=score,fill = ifelse(score>0,'blue','yellow')))+
      geom_bar(stat = 'identity',width=0.3)+theme(legend.position="none")+
      scale_fill_manual(values = c(pal_mk[4],pal_mk[2]))+
      xlab('')+ylab('PC1 score')
    p2 <- ggplot(numdf %>% filter(seqnames==ChrID),aes(x=position,y=score,fill=seqnames))+
      geom_col()+scale_fill_manual(values =matrix(x1)[,1][2])+
      theme(legend.position="none")+xlab('')+ylab('Gene num')
    p3 <- ggplot(Ldf %>% filter(seqnames==ChrID),aes(position,L,fill=seqnames))+geom_col()+
      scale_fill_manual(values = matrix(x1)[,1][6])+
      theme(legend.position = 'none')+xlab('')+ylab('L fpkm')
    p4 <- ggplot(Sdf %>% filter(seqnames==ChrID),aes(position,S,fill=seqnames))+geom_col()+
      scale_fill_manual(values = matrix(x1)[,1][1])+
      theme(legend.position = 'none')+xlab('')+ylab('S fpkm')
    pt <- plot_grid(p1,p2,p3,p4,labels = c('a','b','c','d'),ncol = 1,nrow = 4)
    ggsave(paste(ChrID,'_4tracks.pdf',sep = ''),pt,dpi = 300,width = 20,height = 16)
    ggsave(paste(ChrID,'_4tracks.png',sep = ''),pt,dpi = 300,width = 20,height = 16)
    
  }
  
}


#----根据某个区间画图
seqid <- 'Chr1A'
i_s <- 200000
i_e <- 15500000

iabdf <- ABdf %>% filter(seqnames==seqid,start > i_s,end < i_e)
inumdf <- numdf %>% filter(seqnames==seqid,start > i_s,end < i_e)
iLdf <- Ldf %>% filter(seqnames==seqid,start > i_s,end < i_e)
iSdf <- Sdf %>% filter(seqnames==seqid,start > i_s,end < i_e)
p1 <- ggplot(iabdf,aes(x=position,y=score,fill = ifelse(score>0,'blue','yellow')))+
  geom_col()



#1.2 review1
##-----全局比较AB与fpkm boxplot
library(readxl)
library(ggsignif)
library(tidyverse)
library(wesanderson)
library(tidyquant)


x1 <- palette_dark()
pal_mk <- wes_palette("Darjeeling1")
setwd('D:\\FAFU-CGB\\2021summervacation\\NpX\\20210824_WMJreview\\review2')

adf <- read_excel('NpX.v0804.500kb.Afpkm.xlsx')
adf_L <- adf %>% filter(L>0) %>% mutate(LL=log(L+1)) %>% 
  select(AB,LL) %>% pivot_longer(cols = LL,names_to = 'type',values_to = 'value')
adf_S <- adf %>% filter(S>0)%>% mutate(SS=log(S+1)) %>% 
  select(AB,SS) %>% pivot_longer(cols = SS,names_to = 'type',values_to = 'value')

bdf <- read_excel('NpX.v0804.500kb.Bfpkm.xlsx')
bdf_L <- bdf %>% filter(L>0)%>% mutate(LL=log(L+1)) %>% 
  select(AB,LL) %>% pivot_longer(cols = LL,names_to = 'type',values_to = 'value')
bdf_S <- bdf %>% filter(S>0)%>% mutate(SS=log(S+1)) %>% 
  select(AB,SS) %>% pivot_longer(cols = SS,names_to = 'type',values_to = 'value')
df <- rbind(adf_L,adf_S,bdf_L,bdf_S)


ggplot(df,aes(x=type,y=value,fill=as.factor(AB)))+
  geom_boxplot(width=0.4,position=position_dodge(0.5),
               outlier.colour="black",
               outlier.shape=1,
               outlier.size=0.7,notch=F)+theme_classic()+
  scale_fill_manual(values = c(pal_mk[4],pal_mk[2]))+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = 'top')+labs('')
ggsave('NpX.v0814.500kb.ABfpkm.boxplor.pdf',dpi=300,width = 12,height = 10)
##dodge 箱线图无法计算每组之间的显著性，故先将其分开计算，最后手动添加
compaired <- list(c('A_LL','B_LL'),
                  c('A_SS','B_SS'))
ggplot(df,aes(x=paste(AB,type,sep = '_'),y=value))+
  geom_boxplot()+theme_bw()+
  geom_signif(comparisons = compaired,map_signif_level = TRUE,step_increase = 0.1)+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank())+
  xlab('')
ggsave('NpX.v0814.500kb.ABfpkm.boxplot.sample.pdf',dpi=300,width = 12,height = 10)
##-----单条染色体比对AB与fpkm boxplot
setwd('D:\\FAFU-CGB\\2021summervacation\\NpX\\20210824_WMJreview\\review2')
for (i in 1:10) {
  for (j in c('A','B','C','D')) {
    ChrID <- paste('Chr',i,j,sep = '')
    print(ChrID)
    
  }
  
}

adf <- read_excel(paste('NpX.v0804.500kb.Chr1.Afpkm.xlsx',sep = ''))
ladf <- adf %>% select(chrID,AB,L)
sadf <- adf %>% select(chrID,AB,S)
bdf <- read_excel(paste('NpX.v0804.500kb.Chr1.Bfpkm.xlsx',sep = ''))
lbdf <- bdf %>% select(chrID,AB,S)
lbdf <- bdf %>% select(chrID,AB,S)
for (i in 1:10) {
  ChrID <- paste('Chr',i,sep = '')
  print(ChrID)
  # df = read_excel('NpX.v0804.500kb.Chr1.Afpkm')
  adf <- read_excel(paste('NpX.v0804.500kb.',ChrID,'.Afpkm.xlsx',sep = ''))
  ladf <- adf %>% select(chrID,AB,L)
  sadf <- adf %>% select(chrID,AB,S)
  bdf <- read_excel(paste('NpX.v0804.500kb.',ChrID,'.Bfpkm.xlsx',sep = ''))
  lbdf <- bdf %>% select(chrID,AB,S)
  lbdf <- bdf %>% select(chrID,AB,S)
}



##-----------review3
library(readxl)
library(ggsignif)
library(tidyverse)
library(wesanderson)
library(tidyquant)

x1 <- palette_dark()
pal_mk <- wes_palette("Darjeeling1")

setwd('E:\\2021summervacation\\NpX\\20210824_WMJreview\\review3')

Adf <- read_excel('NpX.v0814.alleleDiff.Afpkm.xlsx')
Bdf <- read_excel('NpX.v0814.alleleDiff.Bfpkm.xlsx')
adf <- Adf %>% filter(L>0,S>0) %>% 
  mutate(type='A',newL = log10(L),newS = log10(S)) %>% 
  select(type,newL,newS) %>% 
  pivot_longer(cols = newL:newS,names_to = 'sample',values_to = 'value')
bdf <- Bdf %>% filter(L>0,S>0) %>%
  mutate(type='B',newL = log10(L),newS = log10(S)) %>% 
  select(type,newL,newS) %>% 
  pivot_longer(cols = newL:newS,names_to = 'sample',values_to = 'value')
df <- rbind(adf,bdf)
head(df)

ggplot(df,aes(x=sample,y=value,fill=type))+
  geom_boxplot(width=0.5,position=position_dodge(0.6))+
  scale_fill_manual(values = c(pal_mk[4],pal_mk[2]))+
  theme_bw()+
  theme(plot.margin = unit(c(1,1,1,1), "cm"),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),legend.position = 'top',
        legend.title  = element_blank(),
        legend.key.size=unit(0.8,"cm"))+
  xlab('')+ylab('log10[FPKM]')
ggsave('NpX.v0814.alleleDiffexp.boxplot.pdf',dpi = 300,width = 8,height = 10)

compaired <- list(c('A_newL','B_newL'),
                  c('A_newS','B_newS'))
ggplot(df,aes(x=paste(type,sample,sep = '_'),y=value,fill=type))+
  geom_boxplot(width=0.5,position=position_dodge(0.6))+
  scale_fill_manual(values = c(pal_mk[4],pal_mk[2]))+
  theme_bw()+
  geom_signif(comparisons = compaired,step_increase = 0.5,map_signif_level = TRUE)+
  theme(plot.margin = unit(c(1,1,1,1), "cm"),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),legend.position = 'none',
        legend.title  = element_blank(),
        legend.key.size=unit(0.8,"cm"))+
  xlab('')+ylab('')
ggsave('NpX.v0814.alleleDiffexp.boxplot.signif.pdf',dpi = 300,width = 8,height = 10)



##--------review4-------
library(readxl)
library(tidyverse)
library(wesanderson)
library(tidyquant)
library(ggsignif)
library(ggprism)


x1 <- palette_dark()
pal_mk <- wes_palette("Darjeeling1")


setwd('E:\\2021summervacation\\NpX\\20210824_WMJreview\\review4')
conAP85 <- read_excel('AP85.matureLS-Conserved.AP85geneID.xlsx')
conAP85 <- conAP85 %>% mutate(FPKM_S = (FPKM_s3+FPKM_s6+FPKM_s9)/3) %>% 
  select(gene_id,FPKM_L,FPKM_S) %>% 
  filter(FPKM_L>0,FPKM_S>0) %>% 
  mutate(L = log10(FPKM_L),S = log10(FPKM_S),type='con',species='AP85') %>% 
  select(species,type,L,S)
a2bAP85 <- read_excel('AP85.matureLS-NpX2AP85.A2B.AP85geneID.xlsx')
a2bAP85 <- a2bAP85 %>% mutate(FPKM_S = (FPKM_s3+FPKM_s6+FPKM_s9)/3) %>% 
  select(gene_id,FPKM_L,FPKM_S) %>% 
  filter(FPKM_L>0,FPKM_S>0) %>% 
  mutate(L = log10(FPKM_L),S = log10(FPKM_S),type='AtoB',species='AP85') %>% 
  select(species,type,L,S)
b2aAP85 <- read_excel('AP85.matureLS-NpX2AP85.B2A.AP85geneID.xlsx')
b2aAP85 <- b2aAP85 %>% mutate(FPKM_S = (FPKM_s3+FPKM_s6+FPKM_s9)/3) %>% 
  select(gene_id,FPKM_L,FPKM_S) %>% 
  filter(FPKM_L>0,FPKM_S>0) %>% 
  mutate(L = log10(FPKM_L),S = log10(FPKM_S),type='BtoA',species='AP85') %>% 
  select(species,type,L,S)

conNpX <- read_excel('NpX.v0804.cleanFPKM-Conserved.NpXgeneID.xlsx')
conNpX <- conNpX %>% filter(L>0,S>0) %>% 
  mutate(L = log10(L),S = log10(S),type='con',species='NpX') %>% 
  select(species,type,L,S)
a2bNpX <- read_excel('NpX.v0804.cleanFPKM-NpX2AP85.A2B.NpXgeneID.xlsx')
a2bNpX <- a2bNpX %>% filter(L>0,S>0) %>%
  mutate(L = log10(L),S = log10(S),type='AtoB',species='NpX') %>% 
  select(species,type,L,S)
b2aNpX <- read_excel('NpX.v0804.cleanFPKM-NpX2AP85.B2A.NpXgeneID.xlsx')
b2aNpX <- b2aNpX %>% filter(L>0,S>0) %>%
  mutate(L = log10(L),S = log10(S),type='BtoA',species='NpX') %>% 
  select(species,type,L,S)

df <- rbind(conAP85,a2bAP85,b2aAP85,conNpX,a2bNpX,b2aNpX)
ggplot(df,aes(x=type,y=L,fill=species))+
  geom_boxplot(width=0.5,position=position_dodge(0.6),size=1.2)+
  scale_fill_manual(values = c(pal_mk[4],pal_mk[2]))+
  theme_prism(border = T)+
  theme(plot.margin = unit(c(1,1,1,1), "cm"),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),legend.position = 'top',
        legend.title  = element_blank(),
        legend.key.size=unit(1.4,"cm"))+
  xlab('')+ylab('log10[FPKM]')
ggsave('NpX.v0804.npx2ap85_ABswitch_withfpkm_L.boxplot.pdf',dpi = 300,width = 8,height = 10)
ggplot(df,aes(x=type,y=S,fill=species))+
  geom_boxplot(width=0.5,position=position_dodge(0.6),size=1.2)+
  scale_fill_manual(values = c(pal_mk[4],pal_mk[2]))+
  theme_prism(border = T)+
  theme(plot.margin = unit(c(1,1,1,1), "cm"),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),legend.position = 'top',
        legend.title  = element_blank(),
        legend.key.size=unit(1.4,"cm"))+
  xlab('')+ylab('log10[FPKM]')
ggsave('NpX.v0804.npx2ap85_ABswitch_withfpkm_S.boxplot.pdf',dpi = 300,width = 8,height = 10)

#做比较
compaired <- list(c('AP85_con','NpX_con'),
                  c('AP85_AtoB','NpX_AtoB'),
                  c('AP85_BtoA','NpX_BtoA'))
# compaired <- list(c('AP85_con','AP85_AtoB'),
#                   c('AP85_AtoB','AP85_BtoA'),
#                   c('AP85_BtoA','AP85_con'),
#                   c('NpX_con','NpX_AtoB'),
#                   c('NpX_con','NpX_BtoA'),
#                   c('NpX_AtoB','NpX_BtoA'))

ggplot(df,aes(x=paste(species,type,sep = '_'),y=L,fill=type))+
  geom_boxplot(width=0.5,position=position_dodge(0.6))+
  # scale_fill_manual(values = c(pal_mk[4],pal_mk[2]))+
  theme_bw()+
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = TRUE)+
  theme(plot.margin = unit(c(1,1,1,1), "cm"),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),legend.position = 'none',
        legend.title  = element_blank(),
        legend.key.size=unit(0.8,"cm"))+
  xlab('')+ylab('')
ggsave('NpX.v0804.npx2ap85_ABswitch_withfpkm_L.signif.boxplot.pdf',dpi = 300,width = 8,height = 10)

ggplot(df,aes(x=paste(species,type,sep = '_'),y=S,fill=type))+
  geom_boxplot(width=0.5,position=position_dodge(0.6))+
  # scale_fill_manual(values = c(pal_mk[4],pal_mk[2]))+
  theme_bw()+
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = TRUE)+
  theme(plot.margin = unit(c(1,1,1,1), "cm"),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),legend.position = 'none',
        legend.title  = element_blank(),
        legend.key.size=unit(0.8,"cm"))+
  xlab('')+ylab('')
ggsave('NpX.v0804.npx2ap85_ABswitch_withfpkm_S.signif.boxplot.pdf',dpi = 300,width = 8,height = 10)



###------------------------------------------------
###2021.08.31
### AP85 HiC reads map to NpX
###------------------------------------------------
library(HiTC)
setwd('E:\\2021summervacation\\NpX\\20210831_AP85maptoNpX')
wrkdir='.'
data_in=importC(paste(wrkdir,"AP85_500000_iced.modified.noctg.matrix",sep="\\"),
                paste(wrkdir,"AP85_500000_abs.noctg.bed",sep="\\"))
fai <- read.table('D:\\FAFU-CGB\\2021summervacation\\NpX\\20210802_matrix_bed\\NpX.v0630.changedID.size.txt',sep = '\t')

sset <- reduce(data_in,chr=c('Chr5A','Chr5B','Chr5C','Chr5D','Chr6A','Chr6B','Chr6C','Chr6D'))
sset <- reduce(data_in,chr=c('Chr5A','Chr6A','Chr6B'))
data_500 <- HTClist(mclapply(sset,binningC,binsize=500000,bin.adjust=FALSE,method='sum',step=1))
mapC(data_500)
mapC(forcePairwise(data_500))



###------------------------------------------------
###2021.09.02
### DLF_NQ Ks reads mapping plot
###------------------------------------------------
library(tidyverse)
library(readxl)
library(ggalt)
library(ggprism)
library(tidyquant)
library(wesanderson)
library(RColorBrewer)
library(ggsci)

colors<-brewer.pal(name="Set3",9)
colors
x1 <- palette_dark()
x2 <- palette_light()
pal_mk <- wes_palette("Darjeeling1")

setwd('E:\\niqiu\\Comparative_genome\\ks')
df <- read_excel('Ks.reads_mapping.xlsx')


#dd <- df %>% pivot_longer(cols = PdPDKaKs:PdXlKaKs,names_to = 'type',values_to = 'value')
dd <- df %>% pivot_longer(cols = c(PdPDKaKs:PdTtKaKs,PdBpKaKs:PdXlKaKs),
                          names_to = 'type',values_to = 'value')
#线图
ggplot(dd,aes(x=index,y=value,color=type),alpha=0.2)+
  geom_xspline(spline_shape = 0.5,size=0.8)+
  #geom_line()+
  theme_prism(border = T,base_size = 14)+
  xlab('Synonymous substitution rate (Ks)')+ylab('Percentage of gene pairs')+
  #scale_color_manual(values = matrix(x2)[,1])+
  #scale_color_manual(values = colors)+
  #scale_color_d3()+
  theme(legend.position = c(0.8,0.8))+
  scale_y_continuous(guide="prism_minor")

ggsave('NQ_9Ks_line_1.pdf',width = 14,height = 10)
ggsave('NQ_9Ks_line_4.pdf',width = 14,height = 10)

#柱状图
ggplot(dd,aes(x=index,y=value,fill=type),color='black',alpha=0.2)+
  #geom_xspline(spline_shape = 0.5,size=0.8)+
  geom_col()+
  theme_prism(border = T,base_size = 14)+
  xlab('Synonymous substitution rate (Ks)')+ylab('Percentage of gene pairs')+
  scale_fill_manual(values = matrix(x1)[,1])+
  #scale_fill_manual(values = colors)+
  #scale_fill_simpsons()+
  theme(legend.position = c(0.8,0.8))+
  scale_y_continuous(guide="prism_minor")

ggsave('NQ_9Ks_bar.pdf',dpi = 300,width = 14,height = 10)


###------------------------------------------------
###2021.09.09
### FLLuo Ks distribution plot
###------------------------------------------------
library(tidyverse)
library(readxl)
library(ggalt)
library(ggprism)
library(tidyquant)
library(wesanderson)
library(RColorBrewer)
library(ggsci)


colors<-brewer.pal(name="Accent",7)
colors
x1 <- palette_dark()
x2 <- palette_light()
pal_mk <- wes_palette("Darjeeling1")

setwd('E:\\FLL\\Comparative_genome\\Ks')
#df <- read_excel('Ks.reads_mapping.xlsx')
df <- read_excel('Ks.distribution.xlsx')


dd <- df %>% pivot_longer(cols = c(BaBaKaKs,BaLnKaKs:BaAfKaKs),
                          names_to = 'type',values_to = 'value')
dd <- df %>% pivot_longer(cols = c(BaBaKaKs:BaPmKaKs),
                          names_to = 'type',values_to = 'value')
colors<-brewer.pal(name="Accent",7)
#线图
order <- c('BaBaKaKs','BaLnKaKs','BaLvKaKs','BaMcKaKs','BaPcKaKs','BaPmKaKs','BaAfKaKs')

ggplot(dd,aes(x=index,y=value,color=fct_relevel(type,order)))+
  geom_xspline(spline_shape = 1.5,size=1.0)+
  #geom_line()+
  theme_prism(border = T,base_size = 14)+
  xlab('Synonymous substitution rate (Ks)')+ylab('Percentage of gene pairs')+
  #scale_color_manual(values = matrix(x1)[,1])+
  scale_color_manual(values = colors)+
  #scale_color_ucscgb()+
  theme(legend.position = c(0.8,0.8))+
  scale_y_continuous(guide="prism_minor")+
  ylim(0,4)



#带有自定义区间的线图
df <- read_excel('Ba_species_ksv1016_4_0.1.xlsx')
df <- read_excel('Ba_species_ksv1016_4.2_0.2.xlsx')
#df <- read_excel('Ba_species_ksv1016_4_0.25.xlsx')
# df <- read_excel('Ks_interval.0.03.xlsx')
# dd <- df %>% pivot_longer(cols = c(BaBa,BaLn:BaAf),
#                           names_to = 'type',values_to = 'value')
dd <- df %>% pivot_longer(cols = c(BaBa,BaLn:BaAf,LnMc:McPc),
                          names_to = 'type',values_to = 'value')
dd <- df %>% pivot_longer(cols = c(BaBa:AfLv,LnMc:McPc),
                          names_to = 'type',values_to = 'value')
# order <- c('BaBa','BaLn','BaLv','BaMc','BaPc','BaPm','BaAf')
dd <- df %>% pivot_longer(cols = BaBa:McPc,
                          names_to = 'type',values_to = 'value')
dd <- df %>% pivot_longer(cols = c(BaBa,BaDr:BaPc,BaAf:McPc),
                          names_to = 'type',values_to = 'value')

ggplot(dd,aes(x=index,y=value,color=type))+
  geom_xspline(spline_shape = 1.5,size=1.1)+
  scale_color_manual(values = matrix(x1)[,1])+
  theme_bw()+theme(legend.position = c(0.8,0.7))+
  scale_x_continuous(limits=c(0,4), breaks=seq(0,4,0.2))+
  scale_y_continuous(limits=c(0,80), breaks=seq(0,80,10))


ggplot(dd,aes(x=index,y=value,color=fct_relevel(type,order)))+
  geom_xspline(spline_shape = 1.5,size=1.2)+
  #geom_line()+
  theme_prism(border = T,base_size = 14)+
  xlab('Synonymous substitution rate (Ks)')+ylab('Percentage of gene pairs')+
  scale_color_manual(values = matrix(x1)[,1])+
  #scale_color_manual(values = colors)+
  #scale_color_ucscgb()+
  theme(legend.position = c(0.8,0.8))+
  scale_y_continuous(guide="prism_minor")+
  ylim(0,4)

ggsave('Ks_interval.0.03.pdf',dpi=300,width = 14,height = 10)


#密度图
library(ggridges)
library(viridis)
df <- read.table('Ks.density.tab',sep='\t',header = TRUE)
colnames(df) <- c('name','Ks')
ggplot(df,aes(x=Ks,fill=name,color=name))+
  geom_density(alpha=0.3)+xlim(0,5)+
  #scale_color_manual(values = colors)+
  #scale_fill_manual(values = colors)+
  #scale_fill_viridis(discrete = TRUE)+
  #scale_color_viridis(discrete = TRUE)+
  theme_prism(border = T,base_size = 14)+
  theme(legend.position = c(0.8,0.7),
        legend.key.size=unit(0.8,"cm"),
        legend.key.width=unit(0.9,"cm"),
        legend.spacing.x=unit(0.3,'cm'),
        legend.background = element_rect(color = "black",size=0.5),
        legend.text = element_text(size = 15),
        legend.margin = margin(t = 4, l = 5, r = 5, b = 4))+
  scale_y_continuous(guide="prism_minor")+
  xlab('Synonymous substitution rate (Ks)')+ylab('Density')
#峰峦图
ggplot(df,aes(x=Ks,y=name,fill=name,color=name))+
  geom_density_ridges(alpha=0.3,scale = 1.6,rel_min_height=0.02)+xlim(0,5)+
  #scale_color_manual(values = matrix(x1)[,1])+
  #scale_fill_manual(values = matrix(x1)[,1])+
  theme_prism(border = T,base_size = 14)+
  theme(legend.position = c(0.95,0.8))


###------------------------------------------------
###2021.09.14
### FLLuo LTR insert time
###------------------------------------------------

library(dplyr)
library(ggplot2)
library(viridis)
library(ggsci)
library(ggprism)
library(viridis)
library(tidyquant)

x1 <- palette_dark()

setwd('E:\\FLL\\Repeat_Sequence\\LTR')
FLL_df <- read.table('HLL.v0711.asm.fasta.pass.list')
Af_df <- read.table('Af.genome.fasta.mod.pass.list')
Dr_df <- read.table('danio_rerio.genome.fasta.mod.pass.list')
Ln_df <- read.table('Lnyas.genome.fasta.pass.list')
Lv_df <- read.table('Lvent.genome.fasta.pass.list')
Mc_df <- read.table('Mcorn.genome.fasta.pass.list')
Pc_df <- read.table('Pcana.genome.fasta.pass.list')

FLL_dd <- FLL_df %>% mutate(species='Ba',MyA=V12/1000000) %>% select(species,MyA) 
Af_dd <- Af_df %>% mutate(species='Af',MyA=V12/1000000) %>% select(species,MyA) 
Dr_dd <- Dr_df %>% mutate(species='Dr',MyA=V12/1000000) %>% select(species,MyA) 
Ln_dd <- Ln_df %>% mutate(species='Ln',MyA=V12/1000000) %>% select(species,MyA) 
Lv_dd <- Lv_df %>% mutate(species='Lv',MyA=V12/1000000) %>% select(species,MyA) 
Mc_dd <- Mc_df %>% mutate(species='Mc',MyA=V12/1000000) %>% select(species,MyA) 
Pc_dd <- Pc_df %>% mutate(species='Pc',MyA=V12/1000000) %>% select(species,MyA) 

df <- rbind(FLL_dd,Af_dd,Dr_dd,Ln_dd,Lv_dd,Mc_dd,Pc_dd)
dff <- df %>% filter(MyA>0)
write.table(df,'FLL.LTR_insertTime.bed',sep = '\t',row.names = FALSE,quote=F)
ggplot(df,aes(MyA,fill=species,color=species))+
  geom_density(adjust=1.5, alpha=0.4)+xlim(0,0.5)+
  theme_prism(border = TRUE)+
  #scale_fill_viridis(discrete = TRUE)+scale_color_viridis(discrete = TRUE)
  #scale_fill_aaas()+scale_color_aaas()
  scale_fill_manual(values = matrix(x1)[,1])+scale_color_manual(values = matrix(x1)[,1])



##将上面的图用山峦图进行展示
library(ggridges)
library(ggplot2)

mycolor <- c("#0077c1","#00a99e","#6bc72b","#ff5a20","#ff1620","#752995","#FF9D1E")
mytheme <- theme(panel.background = element_rect(fill = NA),
                 plot.margin = margin(t=10,r=10,b=5,l=5,unit = "mm"),
                 axis.ticks.y = element_blank(),
                 axis.ticks.x = element_line(colour = "grey40",size = 0.5),
                 axis.line.x = element_line(colour = "grey40",size = 0.5),
                 axis.text.x = element_text(size = 10),
                 axis.title.x = element_text(size = 12),
                 panel.grid.major.y = element_line(colour = rev(mycolor),size = 0.5),
                 panel.grid.major.x = element_blank())

ggplot(df,aes(MyA,species,fill=species,color=species))+
  geom_density_ridges(scale=1.9,alpha = 0.5,rel_min_height=0.01,
                      quantile_lines = TRUE,quantiles = 0.5,
                      vline_size=0.8,vline_linetype="dashed")+
  coord_cartesian(expand=FALSE, clip = "off")+
  ylab('Density')+xlab('Time(Mya)')+
  scale_fill_manual(values = rev(mycolor), guide = "none") +
  scale_color_manual(values = rev(mycolor), guide = "none")+mytheme+
  theme(panel.grid.major.y = element_line(colour = NA,size = 0.5)) 

ggsave('FLL.LTR_insertTime.pdf',dpi=300,width = 14,height = 12)





###------------------------------------------------
###2021.09.24
### 测试juicer call compartment与HiTC的区别
###------------------------------------------------
library(HiTC)
library(tidyverse)
library(ggprism)
library(wesanderson)

setwd('E:\\潘浩然\\脚本测试\\20210923_juicercallAB\\Juicer_result')
pal_mk <- wes_palette("Darjeeling1")

chrID <- 'PGA_SCAFFOLD_40__295_CONTIGS__LENGTH_56163402'
df <- read.table(paste(chrID,'.AB.result',sep = ''),sep = '\t')
df <- mutate(df,position=ifelse(V2<10,V2,V2/1000000))
p1 <- ggplot(df,aes(x=position,y=V4,fill = ifelse(V4>0,'blue','yellow')))+
  geom_bar(stat = 'identity',width=0.5)+
  xlab(chrID)+ylab('')+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank())+
  theme_classic()+
  theme(legend.position="none")+
  theme(legend.position="none")+
  scale_y_continuous(guide = "prism_offset")+
  theme(plot.margin = unit(c(1,1,1,1), "cm"),
        plot.title = element_text(size = 22, face = "bold", #璁剧疆鏍囬瀛椾綋
                                  hjust = 0.5,
                                  margin = margin(b = 15)),
        axis.line = element_line(color = "black"),
        axis.title = element_text(size = 15, color = "black",  #璁剧疆X,Y杞存爣棰樼殑瀛椾綋
                                  face = "bold"),
        axis.text = element_text(size = 12,color = "black"),
        axis.text.y = element_text(size = 15))+
  #scale_fill_manual(values=c("#6cc396","#9896ca"))
  #scale_fill_manual(values = c("#f46f20","#156077"))
  scale_fill_manual(values = c(pal_mk[4],pal_mk[2]))
p1
ggsave(paste(chrID,'compartment.png',sep = '.'),dpi=300,width = 14,height = 10)




###--------------------------------------------------------------------------------------------------
##2021.10.01 NpX文章 flip后对同源染色体之间ALL-conserverd与Non-All-conserved gene做GO富集分析     ##
##---------------------------------------------------------------------------------------------------
library(ggplot2)
library(clusterProfiler)

setwd('E:\\Nepal_project\\20210420_Article_Revision\\BioTriAll_conserverd_enrichment\\flipped')

#GO
go_anno <- read.delim('E:\\Nepal_project\\20200528 GOKEGG\\N-go.no1.annot.txt',sep='\t',header = FALSE)
#KEGG
ko_anno <- read.delim('E:\\Nepal_project\\20200528 GOKEGG\\N-kegg.no1.annot.txt',sep='\t',header = FALSE)


#---All-----
gene_select <- read.delim('AllChr.Allconserved.geneID.rename.txt',sep = '\t',header = FALSE)
gene_select <-as.vector(gene_select$V1)
go_enrich <- enricher(gene = gene_select,
                      TERM2GENE = go_anno[c('V2', 'V1')], 
                      TERM2NAME = go_anno[c('V2', 'V3')], 
                      pvalueCutoff = 0.05, 
                      pAdjustMethod = 'BH')
write.csv(as.data.frame(go_enrich),"All_conserved_GO.csv",row.names = F)

barplot(go_enrich,showCategory=20)
ggsave('All_conserved_GO.pdf',dpi=300,width = 16,height = 12)
ko_enrich <- enricher(gene = gene_select,
                      TERM2GENE = ko_anno[c('V2', 'V1')], 
                      TERM2NAME = ko_anno[c('V2', 'V3')], 
                      pvalueCutoff = 0.05, 
                      pAdjustMethod = 'BH')

write.csv(as.data.frame(ko_enrich),"All_conserved_KO.csv",row.names = F)
dotplot(ko_enrich,showCategory=20)

#---NonAll-----
gene_select <- read.delim('AllChr.NonAllconserved.geneID.rename.txt',sep = '\t',header = FALSE)
gene_select <-as.vector(gene_select$V1)
go_enrich <- enricher(gene = gene_select,
                      TERM2GENE = go_anno[c('V2', 'V1')], 
                      TERM2NAME = go_anno[c('V2', 'V3')], 
                      pvalueCutoff = 0.05, 
                      pAdjustMethod = 'BH')
write.csv(as.data.frame(go_enrich),"NonAllChr_conserved_GO.csv",row.names = F)
barplot(go_enrich,showCategory=20)
ggsave('NonAll_conserved_GO.pdf',dpi=300,width = 16,height = 12)
ko_enrich <- enricher(gene = gene_select,
                      TERM2GENE = ko_anno[c('V2', 'V1')], 
                      TERM2NAME = ko_anno[c('V2', 'V3')], 
                      pvalueCutoff = 0.05, 
                      pAdjustMethod = 'BH')

write.csv(as.data.frame(ko_enrich),"BiTri_conserved_KO.csv",row.names = F)
dotplot(ko_enrich,showCategory=20)


###--------------------
##2021.10.01 Fig3A 部分数据  ##
##---------------------

library(tidyverse)
library(readxl)
library(ggalt)
library(ggprism)
library(tidyquant)
library(wesanderson)
library(RColorBrewer)
library(ggsci)


colors<-brewer.pal(name="Accent",7)
colors
x1 <- palette_dark()
x2 <- palette_light()
pal_mk <- wes_palette("Darjeeling1")


setwd('E:\\中转站')
df <- read_excel('Npx_species_ks.0.01.xlsx')
dd <- df %>% pivot_longer(cols = SoMs:NpNp,names_to = 'type',values_to='value')
ggplot(dd,aes(x=index,y=value,color=type))+
  geom_xspline(spline_shape = 1.5,size=1.1)+
  scale_color_manual(values = matrix(x1)[,1])+theme_bw()+xlim(0,0.3)
ggsave('px_species_ks.0.02.pdf',dpi=300,width = 16,height = 14)



###---------------------------------
##2021.10.19 Np-x ONT reads sort  ##
##----------------------------------
library(ggplot2)
library(plotly)
setwd('E:\\Nepal_project\\20211019_ONTreads_sort')
alignments<- read.table('test.bed', stringsAsFactors = F, fill = T,header = TRUE)
ggplot(alignments)+
  geom_point(mapping = aes(x=refStart2,y=queryStart2),size=0.009)+
  geom_point(mapping = aes(x = refEnd2, y = queryEnd2),size=0.009)+
  geom_segment(aes(
    x = refStart2,
    xend = refEnd2,
    y = queryStart2,
    yend = queryEnd2,
    text = sprintf(
      'Query ID: %s<br>Query Start Pos: %s<br>Query End Pos: %s<br>Target ID: %s<br>Target Start Pos: %s<br>Target End Pos: %s<br>Length: %s kb',
      queryID,
      queryStart,
      queryEnd,
      refID,
      refStart,
      refEnd,
      round(lenAln / 1000, 1)
    )
  ))+
  scale_x_continuous(
                   labels = levels(alignments$refID)) +
  theme_bw()


###--------------------------------------------
##2021.11.04 比较我和王刚师兄谁的注释结果好  ##
##---------------------------------------------
library(tidyverse)
library(ggsci)
library(ggbreak)

setwd('E:\\Badila\\PHR_anno')
a1 <- read.table('PHR.Anno_compare.bed',header = T)
a2 <- read.table('WG.Anno_compare.bed',header = T)

aa1 <- a1 %>% select(gLenDiff,cLenDiff) %>% 
  pivot_longer(cols = c(gLenDiff,cLenDiff),names_to='type',values_to='value') %>% 
  mutate(Author='PHR')

aa2 <- a2 %>% select(gLenDiff,cLenDiff) %>% 
  pivot_longer(cols = c(gLenDiff,cLenDiff),names_to='type',values_to='value') %>% 
  mutate(Author='WG')

head(aa1)
head(aa2)

aa <- rbind(aa1,aa2)
ggplot(aa,aes(x=Author,y=value,fill=type))+geom_boxplot()+
  theme_bw()+ylim(-1000,1000)

#gene length
PgMean <- mean(a1$gLenDiff)
WgMean <- mean(a2$gLenDiff)
PgMean
WgMean
#cds length
PcMean <- mean(a1$cLenDiff)
PcMean <- mean(a2$cLenDiff)
PcMean
PcMean
Author <- c('PHR','WG')
gDiff <- c(PgMean,WgMean)
cDiff <- c(PcMean,PcMean)
TotalgNum <- c(267678,237211)
dd <- data.frame(Author,gDiff,cDiff,TotalgNum) 
dd <- dd %>% 
  pivot_longer(cols=gDiff:TotalgNum,names_to = 'type',values_to = 'values')
head(dd)
# ggplot(dd,aes(x=fct_relevel(type,c('TotalgNum','gDiff','cDiff')),y=values,fill=Author))+
#   geom_bar(stat = 'identity',position = 'dodge',width = 0.5)+
#   #scale_y_continuous(expand=c(0,0))+
#   scale_fill_npg()+theme_classic()

ggplot(dd,aes(x=Author,y=values,fill=Author))+
  geom_bar(stat = 'identity',position = 'dodge',width = 0.5)+
  facet_wrap(~fct_relevel(type,c('TotalgNum','gDiff','cDiff')))+scale_fill_npg()+
  scale_y_break(breaks = c(200, 170000), scales = 0.6)+#, #ticklabels = c(170000, 500000))+
  ylim(0,280000)+theme_bw()+xlab('')+ylab('Values')+theme(legend.position="none")



###---------------------------------------------------
##2021.12.02 跟着NC做批量箱线图+散点图+显著性分析  ##
##---------------------------------------------------
readdf <- function(x){
  return(read_excel('E:\\潘浩然\\文献\\01附代码文献\\NC新冠患者血液贯穿组学特征\\Figure2C.xlsx',sheet=x,skip=1))
}
df_L <- map(seq(1,6),readdf)

Tdf <- function(x){
  return(df_L[[x]] %>% pivot_longer(cols = 2:dim(.)[2],
                                    names_to = 'Sample',
                                    values_to = 'value') %>% 
           pivot_wider(names_from = 'GeneSymbol',values_from = 'value'))
}
df_LT <- map(seq(1:6),Tdf)

alldf_LT <- bind_rows(df_LT)
alldf_LT <- alldf_LT %>% column_to_rownames('Sample')
exp<-log2(alldf_LT+1)
bar_mat<-t(exp)



anno <- read_excel('E:\\潘浩然\\文献\\01附代码文献\\NC新冠患者血液贯穿组学特征\\Figure2C样品分组.xlsx')
anno$type2<-anno$Type

df_p <- exp %>% mutate(gene=rownames(exp)) %>% 
  pivot_longer(cols = 1:dim(exp)[2],names_to = 'sample',values_to = 'value')


mapdf <- function(x){
  return(filter(anno,Sample==x)[1,'Type'] %>% c() %>% .[[1]])
}

tmp <- map_chr(df_p[,'sample'] %>% c() %>% .[[1]],mapdf)
df_p <- df_p %>% mutate(sample=tmp)

color <-c("#5CB85C","#337AB7","#F0AD4E","#D9534F")
compaired <- list(c("Asymptomatic", "Mild"),
                       c("Asymptomatic", "Severe"),
                       c("Asymptomatic", "Critical"),
                       c("Mild", "Severe"),
                       c("Mild", "Critical"),
                       c("Severe", "Critical"))

sample_order <- c("CXCL8","CXCR1","CXCR2","IL1R2","IL1R1",
                  "TLR4","TLR6","MMP8","MMP9","S100A12",
                  "S100A8","CD28","CD3D","CD8A","LCK",
                  "ZAP70","GATA3","EOMES","IL23A","UBE2E3",
                  "NEDD4L","CCDC34","GABARAPL2","PINK1","FOXO3")

ggplot(df_p,aes(x=fct_relevel(sample,"Asymptomatic","Mild",'Severe','Critical'),y=value,fill=sample))+
  geom_boxplot(outlier.shape = NA,alpha=0.7,width=0.5)+geom_point(position="jitter",size=0.8,aes(color=sample),alpha=0.7)+facet_wrap(~fct_relevel(gene,sample_order),scales = "free_y")+
  scale_fill_manual(values = color)+
  scale_color_manual(values = color)+
  theme_bw()+
  geom_signif(comparisons = compaired,map_signif_level = TRUE,step_increase = 0.12)+
  theme(legend.position = 'none')+
  theme(axis.text.x = element_text(size = 10,angle = 45,vjust = 1,hjust = 1))+
  xlab('')+ylab('')


  theme(axis.line=element_line(colour="black"))+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_blank())+

  theme(axis.text.y = element_text(size = 15))+
  theme(legend.position = "NA")

  
  
###---------------------------------------------------
##2021.12.03 统计蔗茅与高粱4号染色体的AB compartm含量
##以及AB switch的比例
##---------------------------------------------------
  
setwd('E:\\蔗茅\\Chr04')
library(tidyverse)
Sball <- read_delim('Sb.result',delim='\t',
                    col_names = c('seqid','start','end','ID','AB'))
Sballratio <- Sball %>% group_by(AB) %>% 
  summarise(n=n()) %>% mutate(ratio=n/sum(n)*100) %>% head(2)

ZMall <- read_delim('ZM.result',delim='\t',
                    col_names = c('seqid','start','end','ID','AB'))
ZMallratio <- ZMall %>% group_by(AB) %>% 
  summarise(n=n()) %>% mutate(ratio=n/sum(n)*100) %>% head(2)

Sbinvert <- read_delim('Sb.invertion.result',delim='\t',
                       col_names = c('seqid','start','end','ID','AB'))

Sbinvertratio <- Sbinvert %>% group_by(AB) %>% 
  summarise(n=n()) %>% mutate(ratio=n/sum(n)*100)

ZMinvert <- read_delim('ZM.invertion.result',delim='\t',
                       col_names = c('seqid','start','end','ID','AB'))


ZMinvertratio <- ZMinvert %>% group_by(AB) %>% 
  summarise(n=n()) %>% mutate(ratio=n/sum(n)*100)

invert_anchors <- read_delim('Srufi.Chr04.Sbicolor.Chr04.invertion.anchors',delim='\t',col_names = c('ZMid','Sbid','Score'))

mapSbdf <- function(x){
  return(filter(Sball,ID==x)[1,'AB'] %>% c() %>% .[[1]])
}
mapZMdf <- function(x){
  return(filter(ZMall,ID==x)[1,'AB'] %>% c() %>% .[[1]])
}
Sbtmp <- map_chr(invert_anchors[,'Sbid'] %>% c() %>% .[[1]],mapSbdf)
ZMtmp <- map_chr(invert_anchors[,'ZMid'] %>% c() %>% .[[1]],mapZMdf)
invert_df <- invert_anchors %>% mutate(Sbab=Sbtmp,ZMab=ZMtmp) 

#判断AB swith类型
add_comma <- function(x, y) {
  if (x == "A" & y == "A") {
    q <- 'Conserved'
  } else if (x=='A' & y== 'B'){
    q <- 'A2B'
  } else if (x=='B' & y== 'A'){
    q <- 'B2A'
  } else {
    q <- 'Conserved'
  }
  return(q)
 }
arg_list <- list(x=invert_df$ZMab,y=invert_df$Sbab)
invertdf <- invert_df %>% mutate(Swithch=unlist(purrr::pmap(arg_list,add_comma)))
invertdf2p <- invertdf %>% group_by(Swithch) %>% summarise(n=n()) %>% 
  mutate(ratio=n/sum(n)*100,type='invert')
all_anchors <- read_delim('Srufi.Chr04.Sbicolor.Chr04.anchors',
                          delim='\t',col_names = c('ZMid','Sbid','Score'),
                          comment = '#')
Sbtmp <- map_chr(all_anchors[,'Sbid'] %>% c() %>% .[[1]],mapSbdf)
ZMtmp <- map_chr(all_anchors[,'ZMid'] %>% c() %>% .[[1]],mapZMdf)
all_df <- all_anchors %>% mutate(Sbab=Sbtmp,ZMab=ZMtmp) 
arg_list <- list(x=all_df$ZMab,y=all_df$Sbab)
alldf <- all_df %>% mutate(Swithch=unlist(purrr::pmap(arg_list,add_comma)))
alldf2p <- alldf %>% group_by(Swithch) %>% summarise(n=n()) %>%
  mutate(ratio=n/sum(n)*100,type='all')

df2p <- rbind(invertdf2p,alldf2p)

# library(RColorBrewer)
# colors<-brewer.pal(name="Set3",2)
library(ggsci)
library(wesanderson)
pal_mk <- wes_palette("Darjeeling1")
ggplot(df2p,aes(x=type,y=ratio,fill=Swithch))+
  geom_bar(stat = 'identity',position = 'stack',width = 0.5)+
  theme_classic()+
  scale_fill_manual(values = c(pal_mk[2],pal_mk[4],pal_mk[5]))+
  scale_y_continuous(expand=c(0,0))



###---------------------------------------------------
##2021.12.18 环棱螺及其近源物种基因组的LAI比较########
##---------------------------------------------------
library(tidyverse)
library(ggsignif)
library(ggsci)
library(wesanderson)
a <- wes_palette("Darjeeling1")
setwd('E:\\FLL\\Repeat_Sequence\\LAI')


compaired <- list(c("Ba", "Dr"), 
                  c("Ba","Lv"), 
                  c("Ba","Mc"),
                  c("Dr","Lv"),
                  c("Dr","Mc"),
                  c("Lv","Mc"))


Ba <- read_delim('HLL.v0711.asm.fasta.out.LAI',delim ='\t')
Ba <- Ba[-1,] %>% mutate(Species='Ba')
head(Ba)

Dr <- read_delim('danio_rerio.genome.fasta.mod.out.LAI',delim='\t')
Dr <- Dr[-1,] %>% mutate(Species='Dr')
head(Dr)

Lv <- read_delim('Lvent.genome.fasta.out.LAI',delim = '\t')
Lv <- Lv[-1,] %>% mutate(Species='Lv')
head(Lv)

Mc <- read_delim('Mcorn.genome.fasta.out.LAI',delim = '\t')
Mc <- Mc[-1,] %>% mutate(Species='Mc')
head(Mc)

df <- rbind(Ba,Dr,Lv,Mc) %>% select(LAI,Species) %>% filter(LAI>0)
ggplot(df,aes(x=Species,y=LAI,fill=Species,color=Species),alpha=0.8)+
  geom_violin(alpha=0.8,width=0.8)+
  geom_boxplot(width=0.1,alpha=0.8,fill='white',outlier.shape = NA)+
  #scale_fill_npg(alpha = 0.8)+
  scale_fill_manual(values = a)+
  scale_color_manual(values = a)+
  geom_signif(comparisons = compaired,map_signif_level = TRUE,step_increase = 0.1)+
  theme_test()+theme_prism(border = T,base_size = 14)+
  xlab('')+theme(legend.position="top")


# 2021.12.26 方棱螺Ks分布图重新绘制 ----------------------------------------------------------
library(tidyverse)
library(readxl)
library(ggalt)
library(ggprism)
library(tidyquant)

x1 <- palette_dark()
x2 <- palette_light()

setwd('E:\\FLL\\Comparative_genome\\Ks')
df <- read_excel('Ba_species_ksv1226_5.2_0.2.xlsx')
dd <- df %>% pivot_longer(cols = BaPc:LvLv,names_to = 'type',values_to = 'value')

p1 <- ggplot(dd,aes(x=index,y=value,color=type))+
  geom_xspline(spline_shape = 1.0,size=1.5)+
  theme_classic()+
  theme_prism(border = T,base_size = 14)+
  xlab('Synonymous substitution rate (Ks)')+ylab('Percentage of gene pairs')+
  scale_color_manual(values = matrix(x1)[,1])+
  #scale_color_manual(values = colors)+
  #scale_color_ucscgb()+
  theme(legend.position = c(0.8,0.8))+
  scale_x_continuous(limits=c(0,5), breaks=seq(0,5,0.5))+
  scale_y_continuous(limits=c(0,80), breaks=seq(0,80,10))
p1
