#usage : "This script is used to use edgeR to perform differential expression analysis without biological duplication"

library(optparse)
option_list = list(make_option(c("-c", "--readscount"), type = "character", default = NULL, help  = "Input your reads count file"),
                   make_option(c("-f", "--firstcol"), type = "integer", default = NULL, help = "Input the column index of the first column to be differentially expressed"),
                   make_option(c("-s", "--secondcol"), type = "integer", default = NULL, help = "Input the column index of the second column to be differentially expressed"),
                   make_option(c("-b", "--bcv"), type = "integer", default = NULL, help = "Input the bcv value"))
                   #make_option(c("-o", "--output"), type = "character", default = NULL, help = "Output the DEG result(.xlsx)"))
opt <- parse_args(OptionParser(option_list=option_list))

library(edgeR)
library(openxlsx)
library(dplyr)
df <- read.delim(opt$readscount,row.names = 1)
df <- df[,-(1:5)]
head(df)
#
df1 <- df[,c(as.numeric(opt$firstcol),as.numeric(opt$secondcol))]
head(df1)
#构建DEGlist
counts <- as.matrix(df1)
group <- 1:2
y <- DGEList(counts=counts, group = group)
#数据过滤
keep <- rowSums(cpm(y)>1) >= 1
y <- y[keep, , keep.lib.sizes=FALSE]
#标准化
y <- calcNormFactors(y)
#差异表达分析
y_bcv <- y
bcv <- as.numeric(opt$bcv)
et <- exactTest(y_bcv, dispersion = bcv ^ 2)
#查看并输出结果
gene1 <- decideTestsDGE(et, p.value = 0.05, lfc = 0)
summary(gene1)
colnames(gene1) <- "Signifi"
results <- cbind(y$counts,et$table,gene1)

DEG <- results[results$Signifi!=0,]
up <- results[results$Signifi==1,]
down <- results[results$Signifi==-1,]

result_name = paste(colnames(df1)[1],colnames(df1)[2],sep = '_vs_')
sheets <- list('DEG_overview_result'=DEG,'up_regulated_gene'=up,'down_regulated_gene'=down)
write.xlsx(sheets,paste(result_name,'.xlsx',sep = ''),row.names = T)