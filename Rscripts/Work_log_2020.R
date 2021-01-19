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



##-----------------------------------------------------
##2020.11.04 对3/4个等位基因上的CRE分布进行柱状图/饼图可视化
##-----------------------------------------------------
library(ggplot2)
library(reshape2)
library(tidyquant)
setwd('E:\\潘浩然\\调控元件生信任务\\数据\\allele_CRE_compare')
df <- read.delim('row_332.allele4_4.info')
df_melt <- melt(df,id.vars='GeneID')
# df_melt <- df_melt[-c(21:24),]
ggplot(df_melt,aes(x=GeneID,fill=variable))+
  geom_bar(aes(y=value),stat = 'identity',color='black',width = 0.55)+
  theme_classic()+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank())+
  scale_x_discrete(labels=unique(df_melt$GeneID))+labs(fill='cis-element')+
  xlab('')+ylab('Number of elements')+scale_y_continuous(expand=c(0,0))+
  # scale_fill_manual(values = c('yellow','darkolivegreen1','lightskyblue','darkgreen',
  #                                                         'deeppink','khaki2','firebrick','brown1','darkorange1',
  #                                                         'cyan1','royalblue4','darksalmon','darkgoldenrod1',
  #                                                         'darkseagreen','darkorchid','#3C5488B2',"#00A087B2",
  #                                                         "#F39B7FB2","#91D1C2B2"))
  scale_fill_manual(values = c("#0055AA","#C40003","#00C19B","#EAC862","#7FD2FF","#007ED3","#B2DF8A","#FFACAA","#FF9D1E",
                               "#C3EF00","#CAB2D6","#894FC6","#2C3E50","#E31A1C","#18BC9C","#CCBE93","#A6CEE3","#1F78B4",
                               "#B2DF8A","#FB9A99","#FDBF6F","#FF7F00"))


rm(list = ls())
setwd('E:\\潘浩然\\调控元件生信任务\\数据\\allele_CRE_compare')
file_path <- paste0(paste0('allele',1:4),'CRE_infotop20.tab')
for (j in seq(1,4)){
  allele <- read.delim(file_path[j])
  allele <- allele %>% mutate(Ratio=Count/sum(allele$Count)*100)
  allele_new <- allele %>% 
    arrange(desc(Site.Name)) %>%
    mutate(ypos = cumsum(Ratio)- 0.5*Ratio )
  ggplot(allele_new, aes(x="", y=Ratio, fill=Site.Name)) +
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    theme_void()  +
    geom_text(aes(y = ypos, label = Site.Name), color = "white", size=5)+
    geom_text(aes(y = ypos+1, label = round(Ratio,2)), color = "black", size=3)+
  scale_fill_manual(values = c("#0055AA","#FF9900","#00C19B","#EAC862","#7FD2FF","#99FFCC","#B2DF8A","#FFACAA","#FFFF66",
                               "#C3EF00","#CAB2D6","#894FC6","#2C3E50","#E31A1C","#bdbdbd","#CC9966","#A6CEE3","#1F78B4",
                               "#33FF66","#FF3366","#CC3399","#8491B4FF"))
  ggsave(paste0('allele',j,'_pie.pdf'),width = 12,height = 9)
  }


# allele <- read.delim('allele1CRE_infotop20.tab')
# allele <- allele %>% mutate(Ratio=Count/sum(allele$Count)*100)
# allele_new <- allele %>% 
#   arrange(desc(Site.Name)) %>%
#   mutate(ypos = cumsum(Ratio)- 0.5*Ratio )
# ggplot(allele_new, aes(x="", y=Ratio, fill=Site.Name)) +
#   geom_bar(stat="identity", width=1, color="white") +
#   coord_polar("y", start=0) +
#   theme_void()  +
#   geom_text(aes(y = ypos, label = Site.Name), color = "white", size=5)+
#   geom_text(aes(y = ypos+1, label = round(Ratio,2)), color = "black", size=4)+
# scale_fill_manual(values = c("#0055AA","#FF9900","#00C19B","#EAC862","#7FD2FF","#99FFCC","#B2DF8A","#FFACAA","#FFFF66",
#                              "#C3EF00","#CAB2D6","#894FC6","#2C3E50","#E31A1C","#bdbdbd","#CC9966","#A6CEE3","#1F78B4",
#                              "#33FF66","#FB9A99","#CC3399","#8491B4FF"))
# ggsave(paste0('allele',j,'.pdf'),width = 12,height = 9)



##-----------------------------------------------------
##2020.11.13 对1/2/3/4个等位基因上的Ks做redas mapping图
##-----------------------------------------------------

library(ggalt)
library(ggplot2)
library(readxl)
library(reshape2)
library(tidyquant)

df <- read_excel('E:\\潘浩然\\调控元件生信任务\\数据\\Ks_reads_mapping\\Ks_reads_mapping.xlsx')
colnames(df)[1] <- 'reads_mapping'
df <- melt(df,id.vars = 'reads_mapping')
x1 <- palette_light()
ggplot(df,aes(x=reads_mapping,y=value,color=variable))+
  geom_xspline(spline_shape = 1,size=1)+
  scale_colour_manual(values = matrix(x1)[,1])+
  theme_classic()+
  theme(legend.position=c(0.8,0.8))+
  theme(text=element_text(size=12))+
  theme(legend.title  =  element_blank())+
  xlab('Synonymous substitution rate (Ks)')+
  ylab('Percentage of gene pairs')+
  scale_y_continuous(expand=c(0,0))+
  scale_x_continuous(expand=c(0,0))
  


##--------------------------------------------------------------------
##2020.11.25 对1/2/3/4个等位的基因在昼夜节律2h上的表达水平Kmeans聚类##
##--------------------------------------------------------------------
library(xlsx)
library(Mfuzz)
df <- read.xlsx('E:\\潘浩然\\调控元件生信任务\\数据\\allele_CRE_compare\\allele_gene\\allele1.Rep.daynight2h_expression.xlsx',sheetIndex=1)
rownames(df) <- df[,1]
df <- df[,-1]
count <- data.matrix(df) 
eset <- new("ExpressionSet",exprs = count)
count <- count[rowMeans(count)>0,]
pseduo <- 0.01
count <- count+pseduo
eset <- filter.std(eset,min.std=1)
eset <- standardise(eset)
c <- 10
#  评估出最佳的m值
m <- mestimate(eset)
# 聚类
cl <- mfuzz(eset, c = c, m = m)
cl$size
mfuzz.plot(eset,cl,mfrow=c(2,5),new.window= FALSE)



##--------------------------------------------------------------------
##2020.11.25 对1/2/3/4个等位的基因的CRE在上游2k水平上的分布##
##这里是将2K分成了四个bin，每个bin500kb                    ##
##--------------------------------------------------------------------
library(ggplot2)
library(reshape2)
df <- read.delim('E:\\潘浩然\\调控元件生信任务\\数据\\allele_CRE_compare\\allele_gene\\ALLallele.Rep.distribution2plot.txt')
dd <- melt(df,id.vars = 'Allele')
ggplot(dd,aes(x=variable,y=value,fill=Allele))+geom_bar(position = 'dodge',stat="identity",color='black')+scale_colour_manual(values = matrix(x1)[,1])+theme_bw()



##------------------------------------------------------------------------
##2020.12.03 对Bru1 1MB片段进行reads maping 统计每个片段的coverage后作图##
##------------------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(ggridges)
Sb <- read.delim('F:\\Bru1\\Sb.coverage',sep = '\t',col.names = c('ctgID','reads','coverage'))
Sr <- read.delim('F:\\Bru1\\Sr.coverage',sep = '\t',col.names = c('ctgID','reads','coverage'))
LA <- read.delim('F:\\Bru1\\Soffic.coverage',sep = '\t',col.names = c('ctgID','reads','coverage'))
Ss <- read.delim('F:\\Bru1\\Sspon.coverage',sep = '\t',col.names = c('ctgID','reads','coverage'))

Sb %<>% mutate(Sb,ID='Sb') %>% select(ID,reads,coverage)
Sr %<>% mutate(Sr,ID='Sr') %>% select(ID,reads,coverage)
LA %<>% mutate(LA,ID='LA') %>% select(ID,reads,coverage)
Ss %<>% mutate(Ss,ID='Ss') %>% select(ID,reads,coverage)
all <- rbind(Sr,LA,Ss,Sb)
ggplot(all_test,aes(x=reads,y=coverage,color=ID))+geom_line()+facet_wrap(ID~.,ncol = 1)+theme_bw()

#分成几个bin进行计算，每个bin包括5wbp。
Sb_bin1 <- Sb[c(1:50000),]
Sr_bin1 <- Sr[c(1:50000),]
LA_bin1 <- LA[c(1:50000),]
Ss_bin1 <- Ss[c(1:50000),]
all_bin1 <- rbind(Sb_bin1,Sr_bin1,LA_bin1,Ss_bin1)
ggplot(all_bin1,aes(x=reads,y=coverage,color=ID))+geom_area()+facet_wrap(ID~.,ncol = 1)+
  scale_colour_manual(values= c("#E31A1C", "#18BC9C", "#A6CEE3","#FF7F00"))+theme_bw()+
  theme(legend.title=element_blank())+theme(legend.position="none")

## 着重看527500-727500这20w bp 的区间
Sb_querybin1 <- Sb[c(527500:577500),]
Sr_querybin1 <- Sr[c(577501:627500),]
LA_querybin1 <- LA[c(627501:677500),]
Ss_querybin1 <- Ss[c(677501:727500),]
all_querybin1 <- rbind(Sb_querybin1,Sr_querybin1,LA_querybin1,Ss_querybin1)
ggplot(all_querybin1,aes(x=reads,y=coverage,color=ID))+geom_area()+facet_wrap(ID~.,ncol = 1,strip.position = 'right')+
  scale_colour_manual(values= c("#E31A1C", "#18BC9C", "#A6CEE3","#FF7F00"))+theme_bw()+
  theme(legend.title=element_blank())+theme(legend.position="none")+theme(axis.text=element_text(size=16),axis.title=element_text(size=20,face="bold"))


#-----------------1206给表达量数据进行log10处理------------------------
library(readxl)
library(dplyr)
library(tibble)
library(magrittr)
# df <- read_excel('C:\\Users\\admin\\Desktop\\allele1.Rep.leaf.xlsx',sheet = 1)
# df %<>% column_to_rownames('gene_id') 
setwd('E:\\潘浩然\\调控元件生信任务\\数据\\allele_CRE_compare\\1206')
df2 <- read_excel('allele4.Rep.leaf.xlsx')
df2 %<>% column_to_rownames('gene_id')
dd <- log10(df2)
write.csv(dd,'allele4.Rep.leaf.log10.csv')




#------------------1207对1/2/3/4 等位基因上保守不保守分布的统计-------------
library(ggplot2)
setwd('C:\\Users\\admin\\Desktop\\20201205\\diff_allele_all_gene')
df <- read.table('genenum_count2plot.txt',header = TRUE)
ggplot(df,aes(x=Allele,y=Num))+geom_bar(aes(fill=Type),stat="identity",position = 'dodge',width = 0.5,color='black')+scale_fill_manual(values=c("#69b3a2","orange"))+theme(legend.position="top")+labs(fill = '')+theme_bw()+theme(axis.text=element_text(size=16),axis.title=element_text(size=20,face="bold"))
#------------1207按照表达量值的高低对全基因组范围内的基因进行聚类-------------
library(pheatmap)
library(readxl)
library(reshape2)
library(dplyr)
library(tibble)
library(magrittr) #为了能使用%<>%
library(ggplot2)
df <- read_excel('E:\\潘浩然\\Result\\Transit\\fpkm_caculate\\cleaned\\全套\\SES_leaf.clean.xlsx')
df %<>% column_to_rownames('gene_id')
# dd <- round(scale(df),2)
dd <- round(log10(df),2)
p <- pheatmap(dd, show_rownames = F, cellwidth =40, cluster_cols = F, 
              cutree_rows = 6,gaps_col = c(2,4,6,8,10,12,14), angle_col = 45,fontsize = 12)
p
ggsave('cluster_base_exp.pdf',width = 50,height = 100)





library(stringr)
library(dplyr)
df1 <- read.table('C:\\Users\\admin\\Desktop\\20201205\\diff_allele_all_gene\\allele1.all_CRE_position.txt')
df1 <-  mutate(df1,Type='A1') %>% select(Type,V1,V2)
df1 <- mutate(df1,bin = str_split_fixed(df1$V2, "bin", 2)[,2]) %>% select(Type,V1,bin)
df1 <- mutate(df1,pos=as.numeric(bin)*100)

df2 <- read.table('C:\\Users\\admin\\Desktop\\20201205\\diff_allele_all_gene\\allele2.all_CRE_position.txt')
df2 <-  mutate(df2,Type='A2') %>% select(Type,V1,V2)
df2 <- mutate(df2,bin = str_split_fixed(df2$V2, "bin", 2)[,2]) %>% select(Type,V1,bin)
df2 <- mutate(df2,pos=as.numeric(bin)*100)

df3 <- read.table('C:\\Users\\admin\\Desktop\\20201205\\diff_allele_all_gene\\allele3.all_CRE_position.txt')
df3 <-  mutate(df3,Type='A3') %>% select(Type,V1,V2)
df3 <- mutate(df3,bin = str_split_fixed(df3$V2, "bin", 2)[,2]) %>% select(Type,V1,bin)
df3 <- mutate(df3,pos=as.numeric(bin)*100)

df4 <- read.table('C:\\Users\\admin\\Desktop\\20201205\\diff_allele_all_gene\\allele4.all_CRE_position.txt')
df4 <-  mutate(df4,Type='A4') %>% select(Type,V1,V2)
df4 <- mutate(df4,bin = str_split_fixed(df4$V2, "bin", 2)[,2]) %>% select(Type,V1,bin)
df4 <- mutate(df4,pos=as.numeric(bin)*100)
all_df <- rbind(df1,df2,df3,df4)
ggplot(all_df,aes(x=Type,y=pos,fill=Type))+geom_boxplot()



#--------------20201209 对叶段基因表达量分组K1/2/3/4 画小提琴图-----------
library(ggplot2)
library(dplyr)
library(readxl)
library(tibble)
library(magrittr)

setwd('E:\\潘浩然\\Result\\Transit\\fpkm_caculate\\cleaned\\全套')
for (i in (1:4)) {
  print(i)
}
df1 <- read_excel('SES_leaf.clean.xlsx',sheet = 2)
df1 %<>% column_to_rownames('gene_id')
df1 %<>% select(Average) 
df1 %<>% log10() 
df1 %<>% mutate(Type='K1') 

df2 <- read_excel('SES_leaf.clean.xlsx',sheet = 3)
df2 %<>% column_to_rownames('gene_id')
df2 %<>% select(Average) 
df2 %<>% log10() 
df2 %<>% mutate(Type='K2')

df3 <- read_excel('SES_leaf.clean.xlsx',sheet = 4)
df3 %<>% column_to_rownames('gene_id')
df3 %<>% select(Average) 
df3 %<>% log10() 
df3 %<>% mutate(Type='K3')

df4 <- read_excel('SES_leaf.clean.xlsx',sheet = 5)
df4 %<>% column_to_rownames('gene_id')
df4 %<>% select(Average) 
df4 %<>% log10() 
df4 %<>% mutate(Type='K4')

df <- rbind(df1,df2,df3,df4)
ggplot(df,aes(x=Type,y=Average,fill=Type))+geom_boxplot()+coord_flip()


#--------20201211对Np-X、LA中各个染色体上不同等位与Sb的近源关系画箱线图---------------
library(readxl)
library(reshape2)
library(dplyr)
library(tibble)
library(magrittr)
library(ggplot2)
setwd('C:\\Users\\admin\\Desktop\\LA_allele_ks')
df <- read_excel('Np-X_allele.ks.xlsx',sheet = 2)
dfch1 <- melt(df)
dfch1 %<>% separate(variable,c('Chr','allele','allele2'),'_')
dfch1 %<>% mutate(Allele=paste(allele,allele2,sep = '_')) %>% select(Chr,Allele,value)
head(dfch1)

for(i in 1:10) { 
  nam <- paste("chr", i, sep = "")
  # nam[i] <- read_excel('Np-X_allele.ks.xlsx',sheet = i+1)
  # head(nam[i])
  # assign(nam,read_excel('Np-X_allele.ks.xlsx',sheet = i+1))
  # assign(nam,melt(nam))
}



#--------对不同等位数目基因上的CRE进行统计------------
library(dplyr)
library(ggplot2)

setwd('C:\\Users\\admin\\Desktop\\20201205\\diff_allele_Rep_gene\\gene2CREnum')
df1con <- read.table('allele1.all.con.gene2freq.bed')
df1con %<>% mutate(Type='A1-Con') %>% select(Type,V2)
df1non <- read.table('allele1.all.non.gene2freq.bed')
df1non %<>% mutate(Type='A1-Non') %>% select(Type,V2)

df2con <- read.table('allele2.all.con.gene2freq.bed')
df2con %<>% mutate(Type='A2-Con') %>% select(Type,V2)
df2non <- read.table('allele2.all.non.gene2freq.bed')
df2non %<>% mutate(Type='A2-Non') %>% select(Type,V2)

df3con <- read.table('allele3.all.con.gene2freq.bed')
df3con %<>% mutate(Type='A3-Con') %>% select(Type,V2)
df3non <- read.table('allele3.all.non.gene2freq.bed')
df3non %<>% mutate(Type='A3-Non') %>% select(Type,V2)

df4con <- read.table('allele4.all.con.gene2freq.bed')
df4con %<>% mutate(Type='A4-Con') %>% select(Type,V2)
df4non <- read.table('allele4.all.non.gene2freq.bed')
df4non %<>% mutate(Type='A4-Non') %>% select(Type,V2)

alldf <- rbind(df1con,df1non,df2con,df2non,df3con,df3non,df4con,df4non)
head(alldf)
alldf %<>% separate(Type,c('type','type2'),'-')
head(alldf)
ggplot(alldf,aes(Type,V2,fill=Type))+geom_boxplot(outlier.colour = NA)+ylim(0,200)+
  stat_boxplot(geom = 'errorbar',width=0.3)+scale_fill_manual(values=c("orange", "#69b3a2"))
