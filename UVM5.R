


###其他周边 
#fig255免疫热图
###免疫图TIMEnewplot####


setwd("/home/aim/Desktop/TCGA_work/TCGA_UVM_project/immnue_heatmap")
library(utils)
library(GSVA)
library(ComplexHeatmap) # 用到里面的pheatmap画图然后拼图，需安装2.8以上版本的ComplexHeatmap
source("/t8a/2022backup/ssy088_202210/ssy088/tutulaile/FigureYa255TIME/pheatmap_translate.R") # 如果不想放弃较老版本的R及其对应的老ComplexHeatmap，可以只加载新版ComplexHeatmap里的这个函数，该脚本出自2.9.4版ComplexHeatmap
library(circlize)
library(viridis)
library(gplots)
library(data.table)
library(estimate)
source("/t8a/2022backup/ssy088_202210/ssy088/tutulaile/FigureYa255TIME/annTrackScale.R") # 加载函数，对数值进行标准化并对极端值做截断
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor


# 加载自定义函数
standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
  outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
  if (!is.null(halfwidth)) {
    outdata[outdata>halfwidth]=halfwidth
    outdata[outdata<(-halfwidth)]= -halfwidth
  }
  return(outdata)
}


# 读入用MOVICS获取的muscle-invasive bladder cancer (MIBC)亚型
sub <- annCol
annCol.tcga <- data.frame(CMOIC=paste0("CS",cmoic.brca$clust.res$clust))
rownames(annCol.tcga) <- cmoic.brca$clust.res$samID

# 读取450K甲基化数据
# 甲基化数据只要用到5个探针的marker就可以，所以这里我的输入数据是简化的，方便传输
library(data.table)
meth <- fread("/t8a/vip39database/database/UCSC_TCGA/TCGA_methylation450/TCGA-UVM.methylation450.tsv.gz",data.table=F,header=T)
meth <- meth%>%
  column_to_rownames("Composite Element REF")

meth[1:4,1:4]

#load("/home/data/vip39/TCGA/Immune_MOVICS_sub/BLCA_meth_metil.rds")
colnames(meth) <- substr(colnames(meth), start = 1,stop = 16)
meth <- meth[,which(substr(colnames(meth), 14, 15) == "01")]
colnames(meth) <- str_sub(colnames(meth),1,15)
meth <- meth[,rownames(annCol.tcga)] 

#meth <- meth.metil

# 匹配亚型
#colnames(meth) <- substr(colnames(meth), start = 1,stop = 16)
#meth <- meth[,which(substr(colnames(meth), 14, 15) == "01")]
#meth <- meth[,rownames(annCol.tcga)] 

MeTIL.marker <- c("cg20792833","cg20425130","cg23642747","cg12069309","cg21554552")
#meth.metil <- meth[MeTIL.marker,]
meth.meti <- beta[MeTIL.marker,]
MeTIL <- beta[MeTIL.marker,rownames(annCol.tcga)]
MeTIL <- t(scale(t(MeTIL)))

# 计算MeTIL得分
MeTIL[1:4,1:4]
#MeTIL1 <- na.omit(MeTIL)
range(MeTIL)

pca.MeTIL <- prcomp(MeTIL,center = F,scale. = F)

MeTIL <- pca.MeTIL$rotation[,1]
annCol.tcga$MeTIL <- MeTIL



# 加载表达数据
colnames(UVM_mRNA_tpm) <- str_sub(colnames(UVM_mRNA_tpm),1,15)
tpm <- UVM_mRNA_tpm[,rownames(annCol.tcga)]
immune.signature <- read.table("/t8a/2022backup/ssy088_202210/ssy088/tutulaile/FigureYa255TIME/Curated_Immune_Cell_Signature.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)

# 构建计算GSVA的列表
cell.type <- unique(immune.signature$CellType)
immune.sig.ccr <- list()
for (i in cell.type) {
  immune.sig.ccr[[i]] <- immune.signature[which(immune.signature$CellType == i),"Symbol"]
}

# 免疫检查点相关基因
imm.targets <- c("CD274","PDCD1","CD247","PDCD1LG2","CTLA4","TNFRSF9","TNFRSF4","TLR9") 

# 免疫细胞的排序
immune.sig.ccr.order <- c("T.cells.CD8", 
                          "T.cells.regulatory..Tregs.",
                          "T.cells.CD4.naive",
                          "T.cells.follicular.helper",
                          "B.cells.naive",
                          "B.cells.memory",
                          "T.cells.gamma.delta",
                          "Dendritic.cells.activated",
                          "Macrophages.M1",
                          "NK.cells.activated",
                          "Plasma.cells",
                          "T.cells.CD4.memory.resting",
                          "T.cells.CD4.memory.activated",
                          "Mast.cells.activated",
                          "NK.cells.resting",
                          "Macrophages.M0",
                          "Macrophages.M2",
                          "Eosinophils",
                          "Monocytes",
                          "Dendritic.cells.resting",
                          "Mast.cells.resting",
                          "Neutrophils",
                          "Endothelial cells",
                          "Fibroblasts")

# 计算immune/stromal得分
indata <- tpm
# 保存到文件
write.table(indata,file = "TCGA_log2TPM_hugo.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

filterCommonGenes(input.f = "TCGA_log2TPM_hugo.txt" , output.f="TCGA_log2TPM_hugo_ESTIMATE.txt", id="GeneSymbol")

estimateScore("TCGA_log2TPM_hugo_ESTIMATE.txt","TCGA_log2TPM_hugo_estimate_score.txt", platform="affymetrix")

est.tcga <- read.table(file = "TCGA_log2TPM_hugo_estimate_score.txt",header = T,row.names = NULL,check.names = F,stringsAsFactors = F,sep = "\t")
rownames(est.tcga) <- est.tcga[,2]; colnames(est.tcga) <- est.tcga[1,]; est.tcga <- est.tcga[-1,c(-1,-2)];
est.tcga <- sapply(est.tcga, as.numeric); rownames(est.tcga) <- c("StromalScore","ImmuneScore","ESTIMATEScore","TumorPurity"); est.tcga.backup = as.data.frame(est.tcga); colnames(est.tcga.backup) <- colnames(indata)

# 对数值进行标准化并对极端值做截断
est.tcga <- annTrackScale(indata = est.tcga, halfwidth = 2, poolsd = F); est.tcga <- as.data.frame(t(est.tcga)) 
rownames(est.tcga) <- colnames(tpm)

tcga.immune.gsva <- gsva(as.matrix(tpm),
                         immune.sig.ccr,
                         method = "gsva")

# 设置颜色
clust.col <- c("#DD492E","#40548A","#32A087","#EC7D21")
heatmap.BlBkRd <- c("#54FEFF","#32ABAA","#125456","#000000","#510000","#A20000","#F30000")
blue <- "#5bc0eb"
gold <- "#ECE700"
cyan <- "#00B3D0"

annCol.tcga$ImmuneScore <- as.numeric(est.tcga[rownames(annCol.tcga),"ImmuneScore"])
annCol.tcga$StromalScore <- as.numeric(est.tcga[rownames(annCol.tcga),"StromalScore"])
annCol.tcga <- annCol.tcga[order(annCol.tcga$CMOIC),]
annColors.tcga <- list() # 构建热图的图例颜色

annColors.tcga[["CMOIC"]] <- c("CS1" = clust.col[1],
                               "CS2" = clust.col[2]
                               #"CS3" = clust.col[3],
                               #"CS4" = clust.col[4]
)
annColors.tcga[["ImmuneScore"]] <- inferno(64)
annColors.tcga[["StromalScore"]] <- viridis(64)

## 热图1：免疫检查点基因表达
indata <- log2(tpm[intersect(rownames(tpm),imm.targets),] + 1)
hm1 <- pheatmap(standarize.fun(indata[imm.targets,rownames(annCol.tcga)],halfwidth = 2), # 表达谱数据标准化
                border_color = NA, # 热图单元格无边框
                annotation_col = annCol.tcga[,c("CMOIC","StromalScore","ImmuneScore")],
                annotation_colors = annColors.tcga[c("CMOIC","StromalScore","ImmuneScore")],
                color = colorpanel(64,low=blue,mid = "black",high=gold),
                show_rownames = T, # 显示行名
                show_colnames = F, # 不显示列名
                cellheight = 12, # 热图高度固定
                cellwidth = 1.5, # 热图宽度固定
                name = "ICI", # 图例名字
                cluster_rows = F, # 行不聚类
                cluster_cols = F) # 列不聚类

#pdf("CheckPoint_heatmap.pdf",width = 10,height = 10)
hm1

dev.off()

## 热图2：肿瘤免疫微环境富集得分
hm2 <- pheatmap(standarize.fun(tcga.immune.gsva[immune.sig.ccr.order,rownames(annCol.tcga)],halfwidth = 1), # 富集得分标准化并排序
                border_color = NA,
                color = NMF:::ccRamp(x = heatmap.BlBkRd,n = 64),
                gaps_row = c(14,22), # 根据不同类别的免疫细胞分割
                show_rownames = T,
                show_colnames = F,
                cellheight = 12,
                cellwidth = 1.5,
                name = "TIME",
                cluster_rows = F,
                cluster_cols = F)

#pdf("TIMEheatmap.pdf",width = 10,height = 10)
hm2

#dev.off()

## 热图3：MeTIL得分
hm3 <- pheatmap(standarize.fun(t(annCol.tcga[,"MeTIL",drop = F]),halfwidth = 1),
                border_color = NA,
                color = NMF:::ccRamp(c(cyan,"black","#F12AFE"),64),
                show_rownames = T,
                show_colnames = F,
                cellheight = 12,
                cellwidth = 1.5,
                name = "MeTIL",
                cluster_rows = F,
                cluster_cols = F)

#pdf("MeTILheatmap.pdf",width = 10,height = 10)
hm3
#dev.off()

# 合并热图并输出
pdf("TIME.pdf",width = 10,height = 10)
draw(hm1 %v% hm2 %v% hm3, 
     heatmap_legend_side = "right", 
     annotation_legend_side = "right")
invisible(dev.off())

##保存分组信息去vip39进行补充绘图
sub_vum <- annCol
sub_vum$clust <- cmoic.brca$clust.res[rownames(annCol),"clust"]
save(sub_vum,file = "~/Desktop/TCGA_work/TCGA_UVM_project/UVM_MOVICS_sub.rds")
##


#fig238 cor risk mut


###药敏fig 282 #
library(PharmacoGx)
library(parallel)
library(dplyr)
library(stringr)
library(tidyverse)
library(tibble)
library(clusterProfiler)
library(ggplot2)
library(ggrepel)
library(ggthemes)
library(gridExtra)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

##UVM药敏单独分析
##CS1预后较差 视为tumor  CS2 为 normal 

setwd("~/Desktop/xiaoyahuitu/FigureYa282CMAP_XSum/")

load("/home/aim/Desktop/TCGA_work/TCGA_UVM_project/UVM_figure282input.rds")
# 读入drug signature
camp_sig <- readRDS("/home/aim/Desktop/xiaoyahuitu/FigureYa282CMAP_XSum/camp_sig.rds")

# 读入disease signature
#dis_sig <- read.csv('dis_sig.csv', sep=',', header=TRUE)

# 读入XSum函数
source("/home/aim/Desktop/xiaoyahuitu/FigureYa282CMAP_XSum/Core_function.R")

# 选择XSum的topN（Yang et al的研究提到topN选择200效果可能比较好，但这个结论可能不适用与cmap的数据，这里我们选择topN = 500）
XLogFC <- eXtremeLogFC(camp_sig, N = 500)

up_gene <- dis_sig$id[dis_sig$fc == 1]
dn_gene <- dis_sig$id[dis_sig$fc == -1]

xsum <- data.frame(score=XSum(XLogFC, up_gene, dn_gene))
xsum <- rownames_to_column(xsum, var = "id")

# 把结果标准化至-1到1（这步也可不做）
xsum_pos <- xsum[xsum$score>0,]
xsum_pos$score <- xsum_pos$score/max(xsum_pos$score)

xsum_neg <- xsum[xsum$score<0,]
xsum_neg$score <- xsum_neg$score/min(xsum_neg$score) * -1

xsum <- rbind(xsum_pos, xsum_neg)

# 将结果从低到高排序
xsum <- xsum[order(xsum$score),]
head(xsum)
write_csv2(xsum,file = "UVM_cs1_drug.csv")

###开始绘图
xsum$number <- 1:nrow(xsum)

# 突出显示top5的药物，标出药物名
select <- xsum[1:5,]

# 开始画图
ggplot(xsum, aes(number,score))+
  geom_point(size=3, color="grey50") + 
  geom_point(data = select, alpha = 1, 
             size = 5, color = "#5ec7dd") + 
  geom_label_repel(data = select, aes(label=id), 
                   color = "white",
                   alpha = 1, point.padding = 1, 
                   size = 5, fill = "#009bc7",
                   segment.size = 1, nudge_x=-0.5, 
                   segment.color = "grey50",
                   direction = "x",
                   hjust = 1) + 
  theme_classic()

####fig 105####

library(limma) # 芯片差异表达
library(impute) # 芯片缺失值多重填补
library(dplyr)
library(pheatmap)
library(gplots)
library(pROC)
library(ggplot2)
display.progress = function ( index, totalN, breakN=20) {
  if ( index %% ceiling(totalN/breakN)  ==0  ) {
    cat(paste(round(index*100/totalN), "% ", sep=""))
  }
} 

# 自定义函数标准化表达谱
standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
  outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
  if (!is.null(halfwidth)) {
    outdata[outdata>halfwidth]=halfwidth
    outdata[outdata<(-halfwidth)]= -halfwidth
  }
  return(outdata)
}


library(rlang)
library(car)
library(pRRophetic)
library(ggplot2)
library(cowplot)

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor


ann_UVM <- sub
ann_UVM$ImmClust <- ifelse(ann_UVM$Cluster=="C1","C1","C2")

dat_UVM <- tumor_exp[,rownames(ann_UVM)]

dat <- dat_UVM
ann <- ann_UVM
#药物名字
GCP.drug <- read.table("/t8a/2022backup/ssy088/tutulaile/FigureYa105GDSC/drug.txt") #如果要例文的两种药物，就换成drug_eg.txt
GCP.drug <- GCP.drug$V1
GCP.drug


sub <- cmoic.brca$clust.res
colnames(sub) <- c("Sample","Cluster")
sub$Cluster <- paste0("C",sub$Cluster)


range(UVM_mRNA_tpm)
#[1]  0.00000 16.60746


ann_UVM <- sub
ann_UVM$ImmClust <- ifelse(ann_UVM$Cluster=="C1","C1","C2")


dat_UVM <- UVM_mRNA_tpm[,sub$Sample]

dat <- dat_UVM
ann <- ann_UVM
GCP.drug

setwd("~/Desktop/TCGA_work/TCGA_UVM_project/UVM_fig105/")
# 自定义足够多的box的颜色，颜色数量至少等于分组数量
jco <- c("#EABF00", "#2874C5", "red")

### 药物敏感性预测 ###、
GCPinfo <- GCP.IC50 <- GCP.expr <- cvOut <- predictedPtype <- predictedBoxdat <- list() # 初始化列表
plotp <- list()

for (drug in GCP.drug) {
  set.seed(1248103) # 因为预测过程默认10-fold CV，所以设置种子以便结果可重复
  cat(drug," starts!\n") # 提示当前药物已开始分析
  
  # 预测IC50值，采用默认参数，详细可参考??pRRopheticPredict参数列表
  predictedPtype[[drug]] <- pRRopheticPredict(testMatrix = as.matrix(dat[,rownames(ann)]),
                                              drug = drug,
                                              tissueType = "allSolidTumors",
                                              selection = 1) # 1表示若有重复基因取均值处理
  
  if(!all(names(predictedPtype[[drug]])==rownames(ann))) {stop("Name mismatched!\n")} # 若名字不匹配则报错退出
  
  predictedBoxdat[[drug]] <- data.frame("est.ic50"=predictedPtype[[drug]],
                                        #"ImmClust"=ifelse(ann$ImmClust == "C3","C3","C12"), # 这里我修改了C1和C2的名
                                        "ImmClust"=ifelse(ann$ImmClust == "C1","C1",ifelse(ann_UVM$Cluster=="C2","C2","C3")), # 这里我修改了C1和C2的名
                                        row.names = names(predictedPtype[[drug]])) 
  #predictedBoxdat[[drug]]$ImmClust <- factor(predictedBoxdat[[drug]]$ImmClust,levels = c("C12","C3"),ordered = T) # 把类改成因子变量
  predictedBoxdat[[drug]]$ImmClust <- factor(predictedBoxdat[[drug]]$ImmClust,levels = c("C1","C2","C3"),ordered = T) # 把类改成因子变量
  # 绘图
  p <- ggplot(data = predictedBoxdat[[drug]], aes(x=ImmClust, y=est.ic50))
  p <- p + geom_boxplot(aes(fill = ImmClust)) + 
    scale_fill_manual(values = jco[1:length(unique(ann$ImmClust))]) + #自定义box的配色
    theme(legend.position="none") + # 倾斜字体
    theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 12),plot.title = element_text(size = 12, hjust = 0.5)) + 
    xlab("") + ylab("Estimated IC50") + 
    ggtitle(drug) # 补上title
  
  plotp[[drug]] <- p # 保存在列表里供合并图片用
  cat(drug," has been finished!\n") # 提示当前药物已分析结束
}


# 合并图片
#适合展示两种药物
p1 <- plot_grid(plotp[[1]],plotp[[2]],labels = c("A","B"),nrow = 1) # title可以AI下拉到合适位置，就如例文所示
ggsave("boxplot of predicted IC50.pdf", width = 6, height = 5)

# 适合展示多种药物
p2 <- plot_grid(plotlist=plotp, ncol=11)
ggsave("./jiaowang/boxplot of predicted IC50_multiple.pdf", width = 25, height = 15)

## 检验组间差异
p <- vector()
for (drug in GCP.drug) {
  tmp <- wilcox.test(as.numeric(predictedBoxdat[[drug]][which(predictedBoxdat[[drug]]$ImmClust %in% "C2"),"est.ic50"]),
                     as.numeric(predictedBoxdat[[drug]][which(predictedBoxdat[[drug]]$ImmClust %in% "C1"),"est.ic50"]),alternative = "less")$p.value
  p <- append(p,tmp) # 两组样本秩和检验p值
}


for (drug in GCP.drug) {
  tmp <- kruskal.test(list(as.numeric(predictedBoxdat[[drug]][which(predictedBoxdat[[drug]]$ImmClust %in% "C1"),"est.ic50"]),
                           as.numeric(predictedBoxdat[[drug]][which(predictedBoxdat[[drug]]$ImmClust %in% "C2"),"est.ic50"]),
                           as.numeric(predictedBoxdat[[drug]][which(predictedBoxdat[[drug]]$ImmClust %in% "C3"),"est.ic50"])))$p.value
  # tmp <- wilcox.test(as.numeric(predictedBoxdat[[drug]][which(predictedBoxdat[[drug]]$ImmClust %in% "C2"),"est.ic50"]),
  #                    as.numeric(predictedBoxdat[[drug]][which(predictedBoxdat[[drug]]$ImmClust %in% "C1"),"est.ic50"]),alternative = "less")$p.value
  p <- append(p,tmp) # 两组样本秩和检验p值
}

names(p) <- GCP.drug
print(p) #打印p值，因为这里有一个不显著，所以当时没有放在boxplot上，有需要的也可以加在ggplot的title里，或参考FigureYa12box直接画在图上。


sort(p,decreasing = F)
names(sort(p,decreasing = F))[1:10]


#保存到文件
write.table(p,"./output_pvalueC2.txt", quote = F, sep = "\t")
#save.image(file = "./all.rds")
#load(file = "./all.rds")


###CS1 CS2作图

plot_grid(plotp[["PF.4708671"]],plotp[["Temsirolimus"]],plotp[["BMS.509744"]],
          plotp[["SB.216763"]],plotp[["PD.173074"]],plotp[["ABT.888"]],
          plotp[["Methotrexate"]],plotp[["KU.55933"]],plotp[["GW843682X"]],plotp[["MK.2206"]],
          nrow = 2)


## 6 X 10比例合适

## 检验组间差异 另一组##
p <- vector()
for (drug in GCP.drug) {
  tmp <- wilcox.test(as.numeric(predictedBoxdat[[drug]][which(predictedBoxdat[[drug]]$ImmClust %in% "C1"),"est.ic50"]),
                     as.numeric(predictedBoxdat[[drug]][which(predictedBoxdat[[drug]]$ImmClust %in% "C2"),"est.ic50"]),alternative = "less")$p.value
  p <- append(p,tmp) # 两组样本秩和检验p值
}
names(p) <- GCP.drug
print(p) #打印p值，因为这里有一个不显著，所以当时没有放在boxplot上，有需要的也可以加在ggplot的title里，或参考FigureYa12box直接画在图上。


sort(p,decreasing = F)
names(sort(p,decreasing = F))[1:10]

#保存到文件
write.table(p,"./output_pvalue_c1.txt", quote = F, sep = "\t")

###top10 p值排序
plot_grid(plotp[["Docetaxel"]],plotp[["BX.795"]],plotp[["Bicalutamide"]],
          plotp[["Cytarabine"]],plotp[["DMOG"]],plotp[["Thapsigargin"]],
          plotp[["Bexarotene"]],plotp[["PF.562271"]],plotp[["X681640"]],
          plotp[["Tipifarnib"]],
          nrow = 2)
ggsave(filename = "./C1_sensitive.pdf",height = 6,width = 10)





###临床信息比较####
library(compareGroups)
##临床三线表

rownames(allclin) <- allclin$sample
class(allclin)
allclin <- as.data.frame(allclin)
allclin[1:4,1:4]
pd <- as.data.frame(pd)
pd[1:4,1:4]

pd$submitter_id.samples <- str_sub(pd$submitter_id.samples,1,15)

comparedata <- data.frame(ID=pd$submitter_id.samples,
                          T=pd$pathologic_T,
                          N=pd$pathologic_N,
                          M=pd$pathologic_M,
                          #grade=pd$neoplasm_histologic_grade,
                          stage=pd$tumor_stage.diagnoses,
                          sex=pd$gender.demographic,
                          age=pd$age_at_initial_pathologic_diagnosis,
                          OS=allclin[pd$submitter_id.samples,"OS"],
                          OS.time=allclin[pd$submitter_id.samples,"OS.time"],
                          PFI=allclin[pd$submitter_id.samples,"PFI"],
                          PFI.time=allclin[pd$submitter_id.samples,"PFI.time"],
                          DSS=allclin[pd$submitter_id.samples,"DSS"],
                          DSS.time=allclin[pd$submitter_id.samples,"DSS.time"]
)
rownames(comparedata) <- comparedata$ID
comparedata <- comparedata[rownames(sub),]
comparedata$group <- sub$Cluster
str(comparedata) # 查看数据集结构
dim(comparedata)
comparedata <- na.omit(comparedata)
comparedata <- comparedata[comparedata$stage!="not reported",]

dim(comparedata)
#comparedata$group <- factor(comparedata$group,levels = c("C1","C2","C3"))
comparedata$group <- factor(comparedata$group,levels = c("C1","C2"))
comparedata$T <- str_sub(comparedata$T,1,2)
table(comparedata$T)
comparedata$T <- factor(comparedata$T,levels = c("T2","T3","T4"))
comparedata$N <- str_sub(comparedata$N,1,2)
table(comparedata$N)
comparedata$N <- factor(comparedata$N,levels = c("NO","NX"))
comparedata$M <- str_sub(comparedata$M,1,2)
table(comparedata$M)
comparedata$M <- factor(comparedata$M,levels = c("MO","M1","MX"))
#comparedata$grade <- str_sub(comparedata$grade,1,2)
#comparedata$grade <- factor(comparedata$grade,levels = c("G1","G2","G3","G4","GX"))
comparedata$stage <- substring(comparedata$stage,7)
table(comparedata$stage)
comparedata$stage <- gsub("iia","ii",comparedata$stage)
comparedata$stage <- gsub("iib","ii",comparedata$stage)
comparedata$stage <- gsub("iiia","iii",comparedata$stage)
comparedata$stage <- gsub("iiib","iii",comparedata$stage)
comparedata$stage <- gsub("iiic","iii",comparedata$stage)


comparedata$stage <- factor(comparedata$stage,levels = c("ii","iii","iv"))
#comparedata$sex <- str_sub(comparedata$M,1,2)
comparedata$sex <- factor(comparedata$sex)
comparedata$age <- as.numeric(comparedata$age)
#comparedata$M <- factor(comparedata$M,levels = c("MO","M1","MX"))
#comparedata$OS <- str_sub(comparedata$N,1,2)
comparedata$OS <- factor(comparedata$OS,levels = c("0","1"))
comparedata$OS.time <- as.numeric(comparedata$OS.time)
comparedata$PFI <- factor(comparedata$PFI,levels = c("0","1"))
comparedata$PFI.time <- as.numeric(comparedata$PFI.time)
comparedata$DSS <- factor(comparedata$DSS,levels = c("0","1"))
comparedata$DSS.time <- as.numeric(comparedata$DSS.time)

descrTable( ~ ., data = comparedata)
descrTable(group ~ ., data = comparedata)


## 先绘制一个基线特征表
restab1 <- descrTable( ~ ., data = comparedata)
restab <- descrTable(group ~ ., data = comparedata)  
restab

export2csv(restab, file='table1.csv')
export2xls(restab, file='table1.xlsx')
export2word(restab, file='table1.docx')
export2pdf(restab, file='table1.pdf')

export2csv(restab1, file='table_all.csv')
export2xls(restab1, file='table_all.xlsx')
export2word(restab1, file='table_all.docx')
export2pdf(restab1, file='table_all.pdf')


####40例SNRPA1验证##
comparedata40 <- comparedata[sample(1:503,40),]
comparedata40 <- comparedata40[-c(1,9:13)]
restab1 <- descrTable( ~ ., data = comparedata40)
export2xls(restab1, file='/home/data/vip39/TCGA/SNRPA1/SNRPA1table1.xlsx')