library(clusterProfiler)
library(optparse)
library(ggplot2)

# 描述参数的解析方式
option_list = list(make_option(c("-g", "--goanno"), type = "character", default = NULL, help  = "Input your GO annotation file"),
                   make_option(c("-l", "--genelist"), type = "character", default = NULL, help = "Input the list of genes you want to enrich"),
                   make_option(c("-k", "--koanno"), type = "character", default = NULL, help = "Input your GO annotation file"),
                   make_option(c("-p", "--pvalueCutoff"), type = "integer", default = NULL, help = "Input the pvalueCutoff value"))

opt <- parse_args(OptionParser(option_list=option_list), usage = "This Script is a test for arguments!")

go_anno <- read.delim(opt$goanno,sep='\t',header = TRUE)
ko_anno <- read.delim(opt$koanno,sep='\t',header = TRUE)

gene_select <- read.delim(opt$genelist,sep='\t',header = FALSE)
gene_select <-as.vector(gene_select$V1)

go_enrich <- enricher(gene = gene_select,
                      TERM2GENE = go_anno[c('GO', 'gene')], 
                      TERM2NAME = go_anno[c('GO', 'description')], 
                      pvalueCutoff = opt$pvalueCutoff, 
                      pAdjustMethod = 'BH')
barplot(go_enrich)
ggsave('GO_barplot.pdf',height = 10,width = 12)
dotplot(go_enrich)
ggsave('GO_dotplot.pdf',height = 10,width = 12)

ko_enrich <- enricher(gene = gene_select,
                      TERM2GENE = ko_anno[c('KO', 'gene')], 
                      TERM2NAME = ko_anno[c('KO', 'description')], 
                      pvalueCutoff = opt$pvalueCutoff, 
                      pAdjustMethod = 'BH')
barplot(ko_enrich)
ggsave('KEGG_barplot.pdf',height = 10,width = 12)
dotplot(ko_enrich)
ggsave('KEGG_dotplot.pdf',height = 10,width = 12)
