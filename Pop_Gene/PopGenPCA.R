#! /usr/bin/env Rscript
# @Author : Haoran Pan
# date: 2021/11/01

library(ggplot2)
library(optparse)
library(ggsci)
library(patchwork)

# option_list = list(make_option(c("-e", "--eigenvec"), type = "character", 
#                                default = NULL, help  = "Input pop eigenvec file"),
#                    make_option(c("-l", "--eigenval"), type = "character",
#                                default = NULL, help = "Input Plink output eigenvec file"),
#                   )

# args <- parse_args(OptionParser(option_list=option_list), 
#                   usage = "This script is used for population genetics PCA analysis")
setwd('E:\\中转站\\NpX群体')

myData <- read.table('NpX_final.snp.pca3.match.eigenvec', sep=' ',
                     header=TRUE, col.names=c('SampleID','group','pc1','pc2','pc3'))
                                              # 'pc4','pc5','pc6','pc7','pc8','pc9',
                                              # 'pc10'))
eigenval <- read.table('NpX_final.snp.pca3.eigenval')

total <- sum(eigenval)

xlabel <- paste('PC1(', round((eigenval$V1[1]/total)*100, 2), '%)', sep='')
ylabel <- paste('PC2(', round((eigenval$V1[2]/total)*100, 2), '%)', sep='')
y1label <- paste('PC3(', round((eigenval$V1[3]/total)*100, 2), '%)', sep='')
p1 <- ggplot(myData,aes(x=pc1,y=pc2,color=group))+
  geom_point(size=2)+ 
  scale_color_npg()+
  xlab(xlabel) + ylab(ylabel)+stat_ellipse(level = 0.95, show.legend = F)+
  theme_bw()
p1


p2 <- ggplot(myData,aes(x=pc1,y=pc3,color=group))+
  geom_point(size=2)+ 
  scale_color_npg()+
  xlab(xlabel) + ylab(y1label)+stat_ellipse(level = 0.90, show.legend = F)+
  theme_bw()

p2

ggsave(paste('NpX_final.snp.pca3.match.eigenvec.PC12', '.pdf', sep=''),
       p1,dpi=300,width = 10,height = 8)
ggsave(paste('NpX_final.snp.pca3.match.eigenvec.PC13', '.pdf', sep=''),
       p2,dpi=300,width = 10,height = 8)
  
p3 <- p1+p2
p3
ggsave(paste('NpX_final.snp.pca3.match.eigenvec', '.pdf', sep=''),
       p3,dpi=300,width = 18,height = 8)
