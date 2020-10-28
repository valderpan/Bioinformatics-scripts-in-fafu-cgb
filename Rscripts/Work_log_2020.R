setwd('D:\\R\\Rproject')
rm(list = ls())
################################
#2020.4.23 AP85 reads_mapping coverage to Nepal 
################################
setwd('C:\\Users\\dell\\Desktop\\Nepal_project\\20200429 reads_mapping\\AP85_mapping_NIP')
library(ggplot2)
df_chr5 <- read.table('AP85_MAP_NIP_W500kb_St100Kb.Chr5.coverage.txt',sep = '\t')
#ggplot(df_chr5,aes(x=V4,y=V5,color = V1))+geom_line(size=0.5)+facet_grid(.~V1)+xlab('Interval position on Chr5')+ylab('coverage')+theme_bw()+theme(legend.title=element_blank())
ggplot(df_chr5,aes(x=V4,y=V5,color = V1))+geom_line(size=1)+facet_grid(rows = vars(V1))+xlab('Interval position in Chr5')+ylab('coverage')+theme_bw()+theme(legend.title=element_blank())+theme(legend.position="none")
#两个物种8个图，分开来画
df_chr5 <- read.table('../all2species.chr5.txt',sep = '\t')
ggplot(df_chr5,aes(x=V4,y=V5,color = V1))+geom_line(size=1)+facet_grid(rows = vars(V1))+xlab('Interval position on Chr5 in Nepal and AP85')+ylab('coverage')+theme_bw()+theme(legend.title=element_blank())+theme(legend.position="none")
df_chr8 <- read.table('all2species.chr8.txt',sep = '\t')
ggplot(df_chr8,aes(x=V4,y=V5,color = V1))+geom_line(size=1)+facet_grid(rows = vars(V1))+xlab('Interval position on Chr8 in Nepal and AP85')+ylab('coverage')+theme_bw()+theme(legend.title=element_blank())+theme(legend.position="none")
#两个物种相同染色体在同一张图，共四张图
df_chr5<- read.table('2species_in1plot.chr5.txt',sep = '\t',header = TRUE)
ggplot(df_chr5,aes(x=num))+geom_line(aes(y=AP85),color = 'darkred',lwd=1)+geom_line(aes(y=Nepal),color = 'darkgreen',lwd=1)+facet_grid(rows = vars(Chr))+theme_bw()+xlab('Interval position on Chr5 in Nepal and AP85')+ylab('coverage')
df_chr8<- read.table('2species_in1plot.chr8.txt',sep = '\t',header = TRUE)
ggplot(df_chr8,aes(x=num))+geom_line(aes(y=AP85),color = 'darkred',lwd=1)+geom_line(aes(y=Nepal),color = 'darkgreen',lwd=1)+facet_grid(rows = vars(Chr))+theme_bw()+xlab('Interval position on Chr8 in Nepal and AP85')+ylab('coverage')


#########################################
#2020.4.27 collienty compartment A to B 
#########################################
setwd('D:\\Result\\TAD')
library(HiTC)
data_in=importC("Nepal_500000.matrix","Nepal_500000_abs.bed")

hic_chr2C <- extractRegion(data_in$Chr2CChr2C,chr="Chr2C",from=0,to=90787606)
pc <- pca.hic(hic_chr2C, normPerExpected=TRUE, method="loess", npc=1)
write.csv(pc$PC1,'Chr2C.compartment.csv')

hic_chr5D <- extractRegion(data_in$Chr5DChr5D,chr="Chr5D",from=0,to=85924935)
pc <- pca.hic(hic_chr5D, normPerExpected=TRUE, method="loess", npc=1)
write.csv(pc$PC1,'Chr5D.compartment.csv')

hic_chr6C <- extractRegion(data_in$Chr6CChr6C,chr="Chr6C",from=0,to=64567053)
pc <- pca.hic(hic_chr6C, normPerExpected=TRUE, method="loess", npc=1)
write.csv(pc$PC1,'Chr6C.compartment.csv')

hic_chr7B <- extractRegion(data_in$Chr7BChr7B,chr="Chr7B",from=0,to=59761120)
pc <- pca.hic(hic_chr7B, normPerExpected=TRUE, method="loess", npc=1)
write.csv(pc$PC1,'Chr7B.compartment.csv')

hic_chr8D <- extractRegion(data_in$Chr8DChr8D,chr="Chr8D",from=0,to=61380726)
pc <- pca.hic(hic_chr8D, normPerExpected=TRUE, method="loess", npc=1)
write.csv(pc$PC1,'Chr8D.compartment.csv')

hic_chr9C <- extractRegion(data_in$Chr9CChr9C,chr="Chr9C",from=0,to=64561144)
pc <- pca.hic(hic_chr9C, normPerExpected=TRUE, method="loess", npc=1)
write.csv(pc$PC1,'Chr9C.compartment.csv')

data_in=importC("AP85_500000.matrix","AP85_500000_abs.bed")

hic_chr2C <- extractRegion(data_in$Chr2CChr2C,chr="Chr2C",from=0,to=126636275)
pc <- pca.hic(hic_chr2C, normPerExpected=TRUE, method="loess", npc=1)
write.csv(pc$PC1,'AP85.Chr2C.compartment.csv')

hic_chr5B <- extractRegion(data_in$Chr5BChr5B,chr="Chr5B",from=0,to=92524012)
pc <- pca.hic(hic_chr5B, normPerExpected=TRUE, method="loess", npc=1)
write.csv(pc$PC1,'AP85.Chr5B.compartment.csv')

hic_chr6A <- extractRegion(data_in$Chr6AChr6A,chr="Chr6A",from=0,to=105899069)
pc <- pca.hic(hic_chr6A, normPerExpected=TRUE, method="loess", npc=1)
write.csv(pc$PC1,'AP85.Chr6A.compartment.csv')

hic_chr7C <- extractRegion(data_in$Chr7CChr7C,chr="Chr7C",from=0,to=92112115)
pc <- pca.hic(hic_chr7C, normPerExpected=TRUE, method="loess", npc=1)
write.csv(pc$PC1,'AP85.Chr7C.compartment.csv')

########################################
# 2020.4.28 lizhen
########################################
setwd('C:\\Users\\dell\\Desktop')
library(xlsx)
df <- read.xlsx('allele_TPM.xlsx',sheetIndex = 3)

########################################
#2020.4.29 Nepal conserved/no_conserved compartment plot
########################################
library(xlsx)
library(reshape2)
data <- read.xlsx('C:\\Users\\dell\\Desktop\\Nepal_project\\20200611 Nepal_compartment保守非保守\\Nepal_conserved_nonconserved_compartment.xlsx',sheetName = 1)
data_melt <- melt(data,id.vars = 'Chr')
ggplot(data_melt,aes(x=as.factor(Chr)))+geom_bar(aes(y=value,fill=variable),stat="identity",position="fill")+theme_bw()+xlab('Chr')+ylab('Proportion')


#######################################
#2020.5.07  Nepal conserved/no_conserved compartment plot and A/B compartment ratio barplot
########################################
df <- read.xlsx('C:\\Users\\dell\\Desktop\\Nepal_project\\20200717 Nepal-AP85-compartment-variation\\Nepal-AP85-compartment-variation.xlsx',sheetIndex = 3)
#ͼһ
df_1 <- read.xlsx('C:\\Users\\dell\\Desktop\\Nepal_project\\20200717 Nepal-AP85-compartment-variation\\df_1_melt.xlsx',sheetIndex = 1)
ggplot(df_1,aes(Type,fill=Part))+geom_bar(aes(y=Value),stat = 'identity',position = 'fill')+
  facet_wrap(~Chr, nrow = 1,strip.position = 'bottom')+theme_bw()+
  theme(panel.border = element_rect(color = 'transparent'))+
  theme(panel.spacing = unit(0,'lines'),strip.background = element_blank(),strip.placement = 'outside')+
  theme(axis.text.x=element_blank())+
  xlab('Chr')+ylab('Proportion')

#Nepal 10 chr A2B/B2A/conserved ratio barplot
df_2 <- select(df,1,6:8)
df_2_melt <- melt(df_2,id.vars = 'Chr')
ggplot(df_2_melt,aes(x=as.factor(Chr)))+geom_bar(aes(y=value,fill=variable),stat="identity",position="fill")+theme_bw()+xlab('Chr')+ylab('Proportion')




############################
##2020.07.17 A2B/B2A/conserved ratio barplot in the evolution of chromosome 2, 5, 6, 7, 8, 9 of the four species
############################
library(ggplot2)
library(xlsx)
library(reshape2)
library(dplyr)

data <- read.xlsx('C:\\Users\\dell\\Desktop\\Nepal_project\\20200717 Nepal-AP85-compartment-variation\\Four-species_compartment_count.xlsx',sheetIndex = 3)
ggplot(data,aes(Species,fill=type))+
  geom_bar(aes(y=value),stat = 'identity',position = 'fill',color='black')+
  facet_wrap(~Chr,nrow = 1,strip.position = 'top')+
  theme_bw()+theme(panel.border = element_rect(color = 'transparent'))+
  theme(panel.spacing = unit(0,'lines'),strip.background = element_blank(),strip.placement = 'outside')+
  theme(axis.text.x=element_text(angle=45,hjust=1))+
  xlab('')+ylab('')+labs(fill = '')
# 把barplot中每个小bar按照演化顺序进行排序
d <-  rep(c(1,2,3),18)
d <- sort(d)
data_id <-  mutate(data,id=d)
ggplot(data_id,aes(id,fill=type))+
  geom_bar(aes(y=value),stat = 'identity',position = 'fill',color='black')+
  facet_wrap(~Chr,nrow = 1,strip.position = 'top')+
  theme_bw()+theme(panel.border = element_rect(color = 'transparent'))+
  theme(panel.spacing = unit(0,'lines'),strip.background = element_blank(),strip.placement = 'outside')+
  xlab('')+ylab('Ratio')+labs(fill = '')




############################
##2020.07.20 Nepal 10 chr ks_median boxplot
############################
library(ggplot2)
library(dplyr)
library(reshape2)
library(xlsx)
setwd('C:\\Users\\dell\\Desktop\\Nepal_project\\20200508 Nepal_chr_Ks\\Median\\nepal_ks_median')
data_chr <- read.xlsx('Nepal_ks_median.xlsx',sheetIndex = 2)
head(data_chr)
data <- melt(data_chr,id.vars = 'Chr')
data <- data[,-2]
ggplot(data,aes(Chr,value,fill=Chr))+
  geom_boxplot()+
  theme_classic()+theme(legend.position="none")+
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  xlab('')+ylab('Ks median')+
  scale_x_discrete(labels=c('Chr1','Chr2','Chr3','Chr4','Chr5','Chr6','Chr7','Chr8','Chr9','Chr10'))
ggsave('Ks_median_boxplot_v1.pdf')

# y轴截取0-0.05
# data <- read.xlsx('C:\\Users\\dell\\Desktop\\Nepal_project\\Nepal_chr_Ks\\Median\\nepal_ks_chr\\Nepal_allele_chr_ks.xlsx',sheetIndex = 4)
data <- read.table('C:\\Users\\dell\\Desktop\\Nepal_project\\20200508 Nepal_chr_Ks\\Median\\nepal_ks_chr\\0_005.txt',sep = '\t',header = TRUE)

ggplot(data,aes(chr,value,fill=chr),color='black')+theme_bw()+
  geom_boxplot()+theme(legend.position="none")+xlab('')+
  scale_x_discrete(labels=c('Chr1','Chr2','Chr3','Chr4','Chr5',
                            'Chr6','Chr7','Chr8','Chr9','Chr10'))+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank())

# 加上显著性检验pvalue值
library(ggpubr)
p <- ggplot(data,aes(chr,value,fill=chr),color='black')+theme_bw()+
  geom_boxplot()+theme(legend.position="none")+xlab('')+
  scale_x_discrete(labels=c('Chr1','Chr2','Chr3','Chr4','Chr5',
                            'Chr6','Chr7','Chr8','Chr9','Chr10'))+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank())
  # p + stat_compare_means(method = "anova", label.y = 0.055)+
  # stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.")
  # p + geom_signif(comparisons = list(c("chr5", "chr8"), c("chr5","chr7"), c("chr1","chr8")),y_position = c(0.052,0.054,0.056),
  #             tip_length = c(0),map_signif_level = T,test = anova)
my_comparisons <- list(c("chr5", "chr8"), c("chr5","chr7"), c("chr1","chr8"))
p+stat_compare_means(comparisons = my_comparisons)+stat_compare_means(label.y = 0.07)

#######################################
## 2020.07.21 Nepal 10chr compartment 4 allele conserved/non-conserved
########################################
library(ggplot2)
library(xlsx)
library(reshape2)
library(dplyr)
df <- read.xlsx('C:\\Users\\dell\\Desktop\\Nepal_project\\Allele-compartment-conserved-non-conserved-count.xlsx',sheetIndex = 2)
df <- melt(df,id.vars = 'chr')
ggplot(df,aes(chr,fill=variable))+
  geom_bar(aes(y=value),stat = 'identity',position = 'fill',color='black')+
  xlab('')+ylab('Ratio')+
  scale_x_discrete(labels=c('Chr1','Chr2','Chr3','Chr4','Chr5','Chr6','Chr7','Chr8','Chr9','Chr10'))+
  theme_bw()+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank())+theme(legend.title  =  element_blank())
    
ggsave('Allele-compartment-conserved-non-conserved-count_barplot.pdf')



#######################################
## 2020.07.21 Ks reads mapping distribution line plot !!!!
########################################
library(ggalt)
library(ggplot2)
library(xlsx)
library(reshape2)

data <- read.xlsx('C:\\Users\\dell\\Desktop\\Nepal_project\\20200722 KS\\allele_minKS_v2\\LA_0_02_reads_mapping.xlsx',sheetIndex = 2)
head(data)
colnames(data)[1] <- 'reads_mapping'
df <- melt(data,id.vars = 'reads_mapping')

ggplot(df,aes(x=reads_mapping,y=value,color=variable))+
  geom_xspline(spline_shape = 0.5,size=1.1)+
  theme_classic()+
  theme(legend.position=c(0.8,0.8))+
  theme(text=element_text(size=13))+
  theme(legend.title  =  element_blank())+
  xlab('Synonymous substitution rate (Ks)')+
  ylab('Percentage of gene pairs')
# change color
cbPalette <- c("#999999", '#FF0000',"#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",'#990099','#00FF00','#0099CC','#333333','#FF0033')
ggplot(df,aes(x=reads_mapping,y=value,color=variable))+
  geom_xspline(spline_shape = 0.5,size=1)+
  theme_bw()+
  theme(legend.position=c(0.8,0.65))+
  theme(text=element_text(size=13))+
  theme(legend.title  =  element_blank())+
  xlab('Synonymous substitution rate (Ks)')+
  ylab('Percentage of gene pairs')+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank())+
  scale_colour_manual(values=cbPalette)

ggsave('ks_reads_mappingv3.pdf')



#######################################
## 2020.07.21 Nepal 10chr ABCD allele AB compartment ratio barplot
########################################

data1 <- read.table('C:\\Users\\dell\\Desktop\\Nepal_project\\20200729 four_species_AB_ratio_info\\Nepal_AB_compartment_ratio_info.txt',sep = '\t')
colnames(data1) <- c('Chr','allele','type','value')

ggplot(data1,aes(x=allele,fill=type))+geom_bar(aes(y=value),stat = 'identity',position = 'fill',color='black')+
  facet_grid(~factor(Chr))+theme_bw()+
  theme(panel.border = element_rect(color = 'transparent'))+
  theme(panel.spacing = unit(0,'lines'),strip.background = element_blank(),strip.placement = 'outside')+
  xlab('')+ylab('Ratio')+
  theme(legend.position="right")+labs(fill = '')+
  scale_fill_manual(values=c("orange", "#69b3a2"))



#######################################
## 2020.08.03 Nepal reads mapping compartment bins reads count barplot
########################################

data <- read.table('C:\\Users\\dell\\Desktop\\Nepal_project\\20200429 reads_mapping\\bins_reads_mapping\\chr_bin_region\\bins_reads\\df2ggplot.txt',sep = '\t',header = TRUE)
colnames(data) <- c('Chr','L','S','bin')
df_L <- data[,-2]
dd <- melt(data,id.vars = c('Chr','bin'))

ggplot(df,aes(bin,S))+geom_bar(aes(y=S),fill='#339966',stat = 'identity',color='black')+theme_bw()+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank())+theme(legend.position="top")+xlab('')+ylab('')


##################################################################
## 2020.08.17 Nepal compartment gene in A2B/B2A/conserved GO/KEGG enrich
##################################################################

### Himalaya-AP85
setwd('C:\\Users\\dell\\Desktop\\Nepal_project\\20200817 four_species_compartment_conserved_GOenrich')
library(clusterProfiler)

go <- read.table('N-go.annot',sep = '\t',quote = "")
kegg <- read.table('Nepal-KEGG_anno.txt',sep = '\t',quote = "")

gene_select <- read.table('Himalaya.AP85.A2B_compartment.HimalayageneID.txt')
#gene_select <- read.table('Himalaya.AP85.B2A_compartment.HimalayageneID.txt')
#gene_select <- read.table('Himalaya.AP85.conserved_compartment.HimalayageneID.txt')
go_term2gene <- data.frame(go$V2,go$V1)
kegg_term2gene <- data.frame(kegg$V2,kegg$V1)
go_term2name <- data.frame(go$V2,go$V3)
kegg_term2name <- data.frame(kegg$V2,kegg$V3)

names(go_term2gene) <- c("go_term","gene")
names(go_term2name) <- c("go_term","name")
names(kegg_term2gene) <- c("ko_term","gene")
names(kegg_term2name) <- c("ko_term","name")

gene_select <- as.vector(gene_select$V1)
go_enrich <- enricher(gene=gene_select,pvalueCutoff = 0.05,
                      pAdjustMethod = "BH",
                      TERM2GENE = go_term2gene,
                      TERM2NAME = go_term2name)
barplot(go_enrich,showCategory=20)
ggsave('Himalaya.AP85.A2B.GO_enrich_barplot.pdf',height = 10,width = 12)
kegg_enrich <- enricher(gene=gene_select,pvalueCutoff = 0.05,
                        pAdjustMethod = "BH",
                        TERM2GENE = kegg_term2gene,
                        TERM2NAME = kegg_term2name)
barplot(kegg_enrich,showCategory=20)
ggsave('Himalaya.AP85.A2B.KEGG_enrich_barplot.pdf',height = 10,width = 15)

##################################################################
## 2020.09.09 ADONIS permutational MANOVA 微生物分析
##################################################################
library(vegan)
setwd('D:\\BaiduNetdiskDownload\\190322-R包vegan执行群落结构差异检验之置换多元方差分析（PERMANOVA）\\data\\Table')
##读入文件
#OTU 丰度表
otu <- read.delim('Table-1(1).txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu <- data.frame(t(otu))
#样本分组文件
group <- read.delim('Table-1-group.txt',sep = '\t', stringsAsFactors = FALSE)
##PERMANOVA 分析（所有分组间比较，即整体差异）
#根据 group$site 这一列分组进行 PERMANOVA 分析检验组间差异，基于 999 次置换，详情 ?adonis

#直接输入 OTU 丰度表，在参数中指定距离类型
#使用 Bray-Curtis 距离测度
adonis_result <- adonis(otu~site, group, distance = 'bray', permutations = 999)
#查看结果
adonis_result
summary(adonis_result)
adonis_result$aov.tab

#可选输出 这里输出的是整体结果
otuput <- data.frame(adonis_result$aov.tab, check.names = FALSE, stringsAsFactors = FALSE)
otuput <- cbind(rownames(otuput), otuput)
names(otuput) <- c('', 'Df', 'Sums of squares', 'Mean squares', 'F.Model', 'Variation (R2)', 'Pr (>F)')
write.table(otuput, file = 'Table_1_PERMANOVA.result_all.txt', row.names = FALSE, sep = '\t', quote = FALSE, na = '')

##组间两两比较
##PERMANOVA 分析（使用循环处理，进行小分组间比较，如两组间）
#推荐使用 OTU 丰度表作为输入数据，每次筛选分组后重新计算样本距离，避免由于样本数减少可能导致的距离变动而造成误差
group_name <- unique(group$site)
adonis_result_two <- NULL
for (i in 1:(length(group_name) - 1)) {
  for (j in (i + 1):length(group_name)) {
    group_ij <- subset(group, site %in% c(group_name[i], group_name[j]))
    otu_ij <- otu[group_ij$names, ]
    adonis_result_otu_ij <- adonis(otu_ij~site, group_ij, permutations = 999, distance = 'bray')     #Bray-Curtis 距离测度，基于 999 次置换
    adonis_result_two <- rbind(adonis_result_two, c(paste(group_name[i], group_name[j], sep = '/'), 'Bray-Curtis', unlist(data.frame(adonis_result_otu_ij$aov.tab, check.names = FALSE)[1, ])))
  }
}
adonis_result_two <- data.frame(adonis_result_two, stringsAsFactors = FALSE)
names(adonis_result_two) <- c('group', 'distance', 'Df', 'Sums of squares', 'Mean squares', 'F.Model', 'Variation (R2)', 'Pr (>F)')
#可选添加 p 值校正，例如 Benjamini 校正
adonis_result_two$'Pr (>F)' <- as.numeric(adonis_result_two$'Pr (>F)')
adonis_result_two$P_adj_BH <- p.adjust(adonis_result_two$'Pr (>F)', method = 'BH')
#输出
write.table(adonis_result_two, 'PERMANOVA.result_two.txt', row.names = FALSE, sep = '\t', quote = FALSE, na = '')


##############################################################################
## 2020.10.22 对含有3/4个基因的等位基因组文件中每一行的3/4个基因的CRE进行绘图
##############################################################################
library(ggplot2)
library(gggenes)
setwd('D:\\调控元件生信任务\\数据\\row_CRE')
df <- read.delim('row_0.bed')
ggplot(df, aes(xmin = Start, xmax = End, y = GeneID, fill = CRE,forward = direction)) +
  geom_gene_arrow() +
  facet_wrap(~ GeneID,scales = "free", ncol = 1) +
  theme_genes()
ggsave('test.pdf',width = 30, height = 6)