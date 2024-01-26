


###UVM 验证队列
#GSE22138  included 63 UVM samples obtained by enucleation in untreated patients.
#GPL570

#GSE84976  28sample  预后数据

#single-cell data set GSE139829

#GSE176345 差异基因  高通量测序数据6个样本 舍弃


####

setwd("/home/aim/Desktop/TCGA_work/TCGA_UVM_project/validation_UVM")
#BiocManager::install("GEOquery")
library(GEOquery)
options(timeout = 10000)
eSet <- getGEO(GEO = 'GSE22138', 
               destdir = '.', 
               getGPL = F)

eSet2 <- getGEO(GEO = 'GSE84976', 
                destdir = '.', 
                getGPL = F)

# eSet3 <- getGEO(GEO = 'GSE176345', 
#                 destdir = '.', 
#                 getGPL = F)



# 提取表达矩阵exp
exp1 <- exprs(eSet[[1]])   #"GPL570"
exp2 <- exprs(eSet2[[1]])   #"GPL10558"
#exp3 <- exprs(eSet3[[1]])   #GPL6244
#exp4 <- exprs(eSet4[[1]])   #GPL6244

exp1[1:4,1:4]
exp2[1:4,1:4]


range(exp1)
#[1]  1.061966 15.625932
range(exp2)
#[1]  4.444015 15.122919

dim(exp1)
dim(exp2)
#exp = log2(exp+1)

# 提取芯片平台编号
gpl1 <- eSet[[1]]@annotation
gpl2 <- eSet2[[1]]@annotation

gpl1 
gpl2
#gpl3
## GPL注释
library(devtools)
#install_github("jmzeng1314/idmap3")
## 下载后本地安装
#devtools::install_local("../idmap3-master.zip")
library(idmap3)
ids_GPL10558=idmap3::get_pipe_IDs('GPL10558')
head(ids_GPL10558)

load("/t8a/database/GEO/GPL570.RData")
ids_GPL570 <- probe_id
head(ids_GPL570)
#方法1 BioconductorR包
gpl 
#http://www.bio-info-trainee.com/1399.html
if(!require(hugene10sttranscriptcluster.db))BiocManager::install("hugene10sttranscriptcluster.db")
library(hugene10sttranscriptcluster.db)
ls("package:hugene10sttranscriptcluster.db")
ids_GPL6244 <- toTable(hugene10sttranscriptclusterSYMBOL)
head(ids_GPL6244)


library(dplyr)

exp1 <- data.frame(exp1)
exp1$probe_id = row.names(exp1)
exp1 <- exp1 %>% 
  inner_join(ids_GPL570,by="probe_id") %>% ##合并探针信息
  dplyr::select(-probe_id) %>% ##去掉多余信息
  dplyr::select(symbol, everything()) %>% #重新排列
  mutate(rowMean =rowMeans(.[grep("GSM", names(.))])) %>% #求出平均数
  arrange(desc(rowMean)) %>% #把表达量的平均值按从大到小排序
  distinct(symbol,.keep_all = T) %>% # 留下第一个symbol
  dplyr::select(-rowMean)  #去除rowMean这一列
exp1[1:4,1:4]

exp2 <- data.frame(exp2)
exp2$probe_id = row.names(exp2)

exp2 <- exp2 %>% 
  inner_join(ids_GPL10558,by="probe_id") %>% ##合并探针信息
  dplyr::select(-probe_id) %>% ##去掉多余信息
  dplyr::select(symbol, everything()) %>% #重新排列
  mutate(rowMean =rowMeans(.[grep("GSM", names(.))])) %>% #求出平均数
  arrange(desc(rowMean)) %>% #把表达量的平均值按从大到小排序
  distinct(symbol,.keep_all = T) %>% # 留下第一个symbol
  dplyr::select(-rowMean)  #去除rowMean这一列
exp2[1:4,1:4]





library(tidyverse)
str(exp1)
exp = exp1 %>% 
  inner_join(exp2,by="symbol") %>% ##合并探针信息
  inner_join(exp3,by="symbol") %>% ##合并探针信息
  inner_join(exp4,by="symbol") %>% ##合并探针信息
  tibble::column_to_rownames(colnames(.)[1]) # 把第一列变成行名并删除
exp[1:4,1:4]
# 先保存一下
# save(exp, eSet, file = "GSE108886.Rdata")
# load('GSE108886.Rdata')
# install.packages("devtools")

# 提取临床信息
library(stringr)
pd1 <- pData(eSet[[1]])
pd2 <- pData(eSet2[[1]])

GSE22138_cli <- pd1
GSE84976_cli <- pd2
GSE22138_cli$fustat <- ifelse(GSE22138_cli$`metastasis:ch1`=="no",0,1)
GSE22138_cli$futime <- as.numeric(GSE22138_cli$`months to endpoint:ch1`)
GSE22138_exp <- exp1


GSE84976_cli$fustat <- if_else(str_detect(GSE84976_cli$`death due to metastasis:ch1`,"alive"),0,1)
GSE84976_cli$futime <- as.numeric(GSE84976_cli$`follow-up (months):ch1`)  ##month
GSE84976_exp <- exp2

save(GSE22138_cli,GSE22138_exp,GSE84976_cli,GSE84976_exp,file = "~/Desktop/TCGA_work/TCGA_UVM_project/validation_UVM/UVM_GSE_validation.rds")
## 筛选诊断为IPF的样本

# pd1 = subset(pd1,characteristics_ch1.1 == 'diagnosis: IPF')
# pd2 = subset(pd2,characteristics_ch1.1 == 'diagnosis: IPF')
# exp_idp = exp[,c(pd1$geo_accession,pd2$geo_accession)]
## 批次校正
#BiocManager::install("sva")
library('sva')
library(TCGAbiolinks)
library(sva)
library(cluster)
library(oompaBase)
if(!require(paletteer))install.packages("paletteer")
if(!require(scico))install.packages('scico')
if(!require(nord))install.packages('nord')
library(paletteer)
paletteer_d("RColorBrewer::Paired")
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

# 加载自定义函数，用于绘制PCA图
source("/home/aim/Desktop/xiaoyahuitu/FigureYa203ComBat/batchPCA.R")
###去除批次并进行差异分析
#FigureYa203
####合并数据集、检查批次效应
# 设置颜色
blue <- "#2874C5"
yellow <- "#EABF00"
green <- "#008B8A"
red <- "#E31A1CFF"
# 提取在三套数据中都出现的基因
comgene <- intersect(intersect(rownames(tcga.expr), rownames(gse41613.expr)), rownames(gse65858.expr))

# 合并三套数据
exp1[1:4,1:4]
exp2[1:4,1:4]
exp3[1:4,1:4]
exp4[1:4,1:4]
GSE108886.exp <- exp1%>%
  column_to_rownames("symbol")
GSE108886.exp[1:4,1:4]
GSE145467.exp <- exp2%>%
  column_to_rownames("symbol")
GSE45885.exp <- exp3%>%
  column_to_rownames("symbol")
GSE45887.exp <- exp4%>%
  column_to_rownames("symbol")

comgene <- intersect(intersect(rownames(GSE108886.exp), rownames(GSE145467.exp)), rownames(GSE45885.exp))
combined.expr <- cbind.data.frame(GSE108886.exp[comgene,],
                                  GSE145467.exp[comgene,],
                                  GSE45885.exp[comgene,],
                                  GSE45887.exp[comgene,]
)

# 绘制PCA散点图，检查批次效应
batchPCA(indata = t(scale(t(combined.expr))),
         batch = rep(c("GSE108886","GSE145467","GSE45885","GSE45887"), times = c(ncol(GSE108886.exp),ncol(GSE145467.exp),ncol(GSE45885.exp),ncol(GSE45887.exp))),
         fig.dir = ".",
         PCA.fig.title = "Raw PCA for combined expression profile",
         cols = c(blue, yellow, green,red),
         showID = F,
         cex = 0.7,
         showLegend = T) # 可以看到三个数据集(批次)完全分开，说明有很严重的批次效应
range(combined.expr)


# 去除批次效应
batch <- data.frame(batch = rep(c("GSE108886","GSE145467","GSE45885","GSE45887"), times = c(ncol(GSE108886.exp),ncol(GSE145467.exp),ncol(GSE45885.exp),ncol(GSE45887.exp))))
modcombat = model.matrix(~1, data = batch)
combined.expr.combat <- as.data.frame(ComBat(dat=as.matrix(combined.expr), batch=batch$batch, mod=modcombat))
# 输出到文件
# 这里保存前三个基因，便于传输和查看
write.csv(combined.expr.combat[1:3,], "output_combined_expr.csv", quote = F)
# 实际使用时请运行下面这行，获得全部基因在三个数据集中的表达矩阵
#write.csv(combined.expr.combat, "output_combined_expr.csv", quote = F)

# 绘制PCA散点图，检查批次效应
batchPCA(indata = t(scale(t(combined.expr.combat))),
         batch = rep(c("GSE108886","GSE145467","GSE45885","GSE45887"),times = c(ncol(GSE108886.exp),ncol(GSE145467.exp),ncol(GSE45885.exp),ncol(GSE45887.exp))),
         fig.dir = ".",
         PCA.fig.title = "Combat PCA for combined expression profile",
         cols = c(blue, yellow, green,red),
         showID = F,
         cex = 0.7,
         showLegend = T) # 可以看到三个数据集(批次)混杂在了一次，说明批次效应被基本消除
range(combined.expr.combat)
combined.expr.combat[1:4,1:4]
# 
# 
# 
# ## 批次信息
# batch = data.frame(sample = c(pd1$geo_accession,pd2$geo_accession),
#                    batch = c(pd1$platform_id,pd2$platform_id))
# 
# ## 未批次校正前PCA
# #install.packages('FactoMineR')
# #install.packages('factoextra')
# library("FactoMineR")
# library("factoextra")
# pca.plot = function(dat,col){
#   
#   df.pca <- PCA(t(dat), graph = FALSE)
#   fviz_pca_ind(df.pca,
#                geom.ind = "point",
#                col.ind = col ,
#                addEllipses = TRUE,
#                legend.title = "Groups"
#   )
# }
# pca.plot(exp_idp,factor(batch$batch))
# ## sva 批次校正
# 
# combat_exp <- ComBat(dat = as.matrix(log2(exp_idp+1)),
#                      batch = batch$batch)
# pca.plot(combat_exp,factor(batch$batch))

combat_exp = data.frame(combined.expr.combat)
save(combat_exp, eSet,eSet2,eSet3,eSet4, file = "NOA4datasets.Rdata")

# 载入保存的数据
load("NOA4datasets.Rdata")

# 提取临床信息
pd1 <- pData(eSet[[1]])
pd2 <- pData(eSet2[[1]])
pd3 <- pData(eSet3[[1]])
pd4 <- pData(eSet4[[1]])
## 筛选诊断为IPF的样本
library(stringr)
pd1 <- data.frame(GEO_id=pd1$geo_accession,
                  Group=ifelse(str_detect(pd1$characteristics_ch1.2,"Non"),"NOA","Control"))
pd2 <- data.frame(GEO_id=pd2$geo_accession,
                  Group=ifelse(str_detect(pd2$title,"impaired"),"NOA","Control"))
pd3 <- data.frame(GEO_id=pd3$geo_accession,
                  Group=ifelse(str_detect(pd3$title,"NOA"),"NOA","Control"))
pd4 <- data.frame(GEO_id=pd4$geo_accession,
                  Group=ifelse(str_detect(pd4$title,"NOA"),"NOA","Control"))
pd <- rbind(pd1,pd2,pd3,pd4)
table(pd$GEO_id==colnames(combat_exp))

pd5 <- data.frame(GEO_id=pd5$geo_accession,
                  Group=ifelse(str_detect(pd5$title,"NOA"),"NOA","Control"))