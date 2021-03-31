#! /usr/bin/env Rscript
# @Author : Haoran Pan
# date: 2021/3/31


#=================================================
#ATAC-Seq InsertFragmentSize stats
#Input：samtools view -f 0x40 ${i}.sorted.dedup.bam | perl -ane 'print abs($F[8]);print "\n";' | sort -n > ${i}.insertSize.txt
#=================================================

library(ggplot2)
library(dplyr)
library(optparse)
library(ggprism)

rm(list = ls())
options(scipen = 999)
option_list <- list(make_option(c("-i","--input"),type="character", default = NULL,help  = "Input ATAC insertSize file"),
                    make_option(c("-o","--output"),type="character",default = NULL,help  = "Output the PDF file"),
                    make_option(c("-n","--name"),type="character", default = NULL,help="Experiment name, which will be used to generate output plot title names"))

args <- parse_args(OptionParser(option_list=option_list))

df <- read.table(args$input)
dd <- df %>% group_by(V1) %>% summarise(count=n()) %>% mutate(count = count/1000)

ggplot(dd,aes(x=V1,y=count))+
  geom_bar(stat = 'identity', fill="#FF7F00", alpha=.8,width=1)+
  theme_bw()+xlab('Fragment length(bp)')+ylab('Read Count (x 10^3)')+
  scale_y_continuous(guide = guide_prism_minor())+
  theme_prism(base_size = 10)+
  theme_prism(base_fontface = "bold")+
  theme_prism(border = TRUE) + coord_cartesian(clip = "off")+
  ggtitle(args$name) +
  ggsave(args$output,width = 10,height=8)
