#! /usr/bin/env Rscript
# @Author : Haoran Pan
# date: 2021/04/03

library(ggplot2)
library(HiTC)
options(scipen = 999)

hic <- importC('merged_500000.matrix','merged_500000_abs.bed')
head(summary(hic))
detail(hic$Chr02Chr02) #TODO

Chr_query <- extractRegion(hic$Chr5AChr5A,chr = 'Chr5A',from=1,to=68253030) #TODO
Chr_query_norm <- normPerExpected(Chr_query, method="loess")
pc_norm <- pca.hic(Chr_query_norm, normPerExpected=T, method="loess", npc=1)

ggplot(as.data.frame(pc_norm),aes(x=start,y=score,fill = ifelse(score>0,'blue','yellow')))+
  geom_bar(stat = 'identity')+theme_bw()+
  theme(legend.position="none")+
  scale_fill_manual(values=c("orange", "#69b3a2"))+
  xlab('Chr5A') #TODO

ggsave('hitc_callAB.pdf',height = 4,width = 12, units = "in", dpi = 300) #TODO