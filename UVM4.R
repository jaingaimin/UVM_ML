
#save.image(file = "/t8a/allimage_data/UVMalldata202211.rds")
load("/t8a/allimage_data/UVMalldata202211.rds")
###UVM外部验证 及临床信息的作图 ROC 等

load("~/Desktop/TCGA_work/TCGA_UVM_project/validation_UVM/UVM_GSE_validation.rds")
GSE22138_exp[1:4,1:4]
library(tidyverse)
GSE22138_exp <- GSE22138_exp%>%
  as.data.frame()%>%
  column_to_rownames("symbol")
table(rownames(GSE22138_cli)==colnames(GSE22138_exp))
table(rownames(GSE22138_cli)%in%colnames(GSE22138_exp))
GSE22138_exp <- GSE22138_exp[,rownames(GSE22138_cli)]
range(GSE22138_cli$futime)
GSE22138_cli$futime <- GSE22138_cli$futime*30


GSE22138 <- list(mRNA.expr=GSE22138_exp,
                 clin.info=GSE22138_cli)


# run NTP in Yau cohort by using up-regulated biomarkers
GSE22138.ntp.pred <- runNTP(expr       = GSE22138$mRNA.expr,
                       templates  = marker.up$templates, # the template has been already prepared in runMarker()
                       scaleFlag  = TRUE, # scale input data (by default)
                       centerFlag = TRUE, # center input data (by default)
                       doPlot     = TRUE, # to generate heatmap
                       fig.name   = "NTP HEATMAP FOR GSE22138") 
head(GSE22138.ntp.pred$ntp.res)
##外部数据生存分析
# compare survival outcome in Yau cohort

surv.GSE22138 <- compSurv(moic.res         = GSE22138.ntp.pred,
                     surv.info        = GSE22138$clin.info,
                     convt.time       = "m", # switch to month
                     surv.median.line = "hv", # switch to both
                     xyrs.est         = c(3,5),
                     fig.name         = "KAPLAN-MEIER CURVE OF NTP FOR GSE22138") 




GSE84976_exp[1:4,1:4]
library(tidyverse)
GSE84976_exp <- GSE84976_exp%>%
  as.data.frame()%>%
  column_to_rownames("symbol")
table(rownames(GSE84976_cli)==colnames(GSE84976_exp))
table(rownames(GSE84976_cli)%in%colnames(GSE84976_exp))
GSE84976_exp <- GSE84976_exp[,rownames(GSE84976_cli)]
range(GSE84976_cli$futime)

GSE84976_cli$futime <- GSE84976_cli$futime*30


GSE84976 <- list(mRNA.expr=GSE84976_exp,
                 clin.info=GSE84976_cli)


# run NTP in Yau cohort by using up-regulated biomarkers
GSE84976.ntp.pred <- runNTP(expr       = GSE84976$mRNA.expr,
                            templates  = marker.up$templates, # the template has been already prepared in runMarker()
                            scaleFlag  = TRUE, # scale input data (by default)
                            centerFlag = TRUE, # center input data (by default)
                            doPlot     = TRUE, # to generate heatmap
                            fig.name   = "NTP HEATMAP FOR GSE84976") 
head(GSE84976.ntp.pred$ntp.res)
##外部数据生存分析
# compare survival outcome in Yau cohort

surv.GSE84976 <- compSurv(moic.res         = GSE84976.ntp.pred,
                          surv.info        = GSE84976$clin.info,
                          convt.time       = "m", # switch to month
                          surv.median.line = "hv", # switch to both
                          xyrs.est         = c(3,5),
                          fig.name         = "KAPLAN-MEIER CURVE OF NTP FOR GSE84976") 


###其他临床信息比较  柱状图加点图

table(rownames(GSE22138.ntp.pred$clust.res)==rownames(GSE22138_cli))
# TRUE 
# 63 
GSE22138_cli$Group <- paste0("CS",GSE22138.ntp.pred$clust.res$clust)

options(stringsAsFactors = F)
library(ggplot2)

data1<-GSE22138_cli
table(data1$Group)
data1$Group<- factor(data1$Group,
                              levels=c("CS1","CS2")) ##将组别转换为factor，方便后续作图

cbPalette <- c("#E69F00", "#CC79A7", "#56B4E9", "#009E73", "#CC79A7", "#F0E442", "#999999","#0072B2","#D55E00")

colnames(data1)
data1$tumor_thickness <- as.numeric(data1$`tumor thickness (mm):ch1`)
p <- ggplot(data=data1,aes(x=Group,
                          y=tumor_thickness,
                          color=Group))+
  geom_jitter(alpha=1,
              position=position_jitterdodge(jitter.width = 0.35, 
                                            jitter.height = 0, 
                                            dodge.width = 0.8))+
  geom_boxplot(alpha=0.2,width=0.45,
               position=position_dodge(width=0.8),
               size=0.75,outlier.colour = NA)+
  geom_violin(alpha=0.2,width=0.9,
              position=position_dodge(width=0.8),
              size=0.75)+
  scale_color_manual(values = cbPalette)+
  theme_classic() +
  theme(legend.position="none") + ##图例设置为没有‘none’
  theme(text = element_text(size=16)) + 
  #ylim(0.0,1.3)+ ##用来规定y轴的标尺
  ylab("tumor thickness (mm)")+
  #scale_x_discrete(labels=c("A","B","C","D"))+ ##可以重置x轴你想要展示的名称
  annotate("segment", x = 1-0.01, y = 17, xend = 2.01,lineend = "round", 
           yend = 17,size=1,colour="black",arrow = arrow(length = unit(0.02, "npc")))+##添加辅助线的位置
  annotate("segment", x = 2.01, y = 17, xend = 0.99,lineend = "round", 
           yend = 17,size=1,colour="black",arrow = arrow(length = unit(0.02, "npc")))+##添加辅助线的位置
  annotate("text", x=1.5,y=17.5, 
           label=expression("**"~"FDR"~2.41%*%10^-10),vjust=0)##添加符号和文字
p

##通过ggpubr计算p值
library(ggpubr)


###构建分型score 模型分析 及评估 ####
#参考脚本 fig204 201
library(Boruta)
library(ggplot2)
library(ggpubr)
library(pheatmap)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor


# 加载表达谱
#expr <- read.table("easy_input_expr.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
#tumsam <- colnames(expr)[substr(colnames(expr),11,12) == "01"]
load("/t8a/vip39database/database/tcga_counts_fpkm_tpm/TCGA-UVM_tpm_gene_symbol.Rdata")
tpms <- as.data.frame(tpms)
range(tpms)
colnames(tpms) <- str_sub(colnames(tpms),1,15)
tpms[1:4,1:4]
tpms <- na.omit(tpms)
#提取mRNA矩阵
expr <- tpms[mRNA_anno$gene_name,cmoic.brca$clust.res$samID]
range(expr)
#[1]     0.0 94062.9
indata <- t(scale(t(log2(expr + 1)))) # z-score表达谱
table(is.na(indata))
range(indata)
#[1] -6.081861 10.431476

# 加载亚型结果
annCol <- read.table("easy_input_cluster.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
annCol <- cmoic.brca$clust.res
annCol$`Gene cluster` <- ifelse(annCol$clust==1,"A","B")
table(annCol$`Gene cluster`)

# 加载Figure201ClusterCorrelation的ICI signature gene结果
##开始计算相关性 fig201
geneclust <- cmoic.brca$clust.res

table(geneclust$clust)
samorder <- geneclust$clust
names(samorder) <- cmoic.brca$clust.res$samID
samorder <- sort(samorder)


names(samorder)

table(names(samorder)%in%colnames(indata))
##只需要用于聚类的gene 


indata <- as.data.frame(indata)
indata <- na.omit(indata)
table(c(unique(marker.up$templates$probe,marker.up$templates$probe))%in%rownames(indata))
gene_need <- intersect(rownames(indata),c(unique(marker.up$templates$probe,marker.up$templates$probe)))
indata <- indata[gene_need,]
outTab <- NULL

for (i in rownames(indata)) {
  tmp <- as.numeric(indata[i,names(samorder)])
  cor.res <- cor.test(tmp, as.numeric(samorder), method = "pearson")
  outTab <- rbind.data.frame(outTab,
                             data.frame(gene = i,
                                        r = cor.res$estimate,
                                        p = cor.res$p.value,
                                        stringsAsFactors = F),
                             stringsAsFactors = F)
}


# 按相关性正负来分类
outTab$direct <- ifelse(outTab$r > 0, "A","B") # 正相关标为A，否则标为B
outTab <- outTab[order(outTab$r, decreasing = T),]
table(outTab$direct) # 事实上最终用于展示热图的基因数目更少，因为作者进一步采用了随机森林降维
# A  B 
# 96 91
# 把基因与cluster的相关系数r、p值及其所属的基因集（A、B）保存到文件
write.table(outTab,"ouput_ICIsignatureGene.txt",sep = "\t", row.names = F, col.names = T, quote = F)

annCol <- data.frame("Gene cluster" = ifelse(samorder == 1,"A", ifelse(samorder == 2, "B", "C")),
                     row.names = names(samorder),
                     check.names = F,
                     stringsAsFactors = F)
annCol <- cbind.data.frame(annCol,clin.info[rownames(annCol),])
colnames(allclin)
table(rownames(annCol)%in%rownames(allclin))
allclin <- as.data.frame(allclin)
annCol$Gender <- allclin[rownames(annCol),"gender"]
rownames(allclin) <- allclin$sample
annCol$Grade <- annCol$pstage
annCol$Status <- ifelse(annCol$fustat==0,"Dead","Alive")
annCol$Age <- annCol$age

annRow <- data.frame("ICI signature gene" = outTab$direct,
                     row.names = outTab$gene,
                     check.names = F,
                     stringsAsFactors = F)

annColors <- list("Gene cluster" = c("A" = "#008ECB", "B" = "#EA921D"),
                  "Gender" = c("MALE" = "#79B789", "FEMALE" = "#B5262A"),
                  "Status" = c("Alive" = "#79B789", "Dead" = "#B5262A"),
                  "Grade" = c( "T2" = "#53B1E7", "T3" = "#C78157", "T4" = "#A54D48"),
                  "Age" = colorRampPalette(c("#F7D202", "#96862A"))(64),
                  "ICI signature gene" = c("A" = "#D14039", "B" = "#008ECB"))

plotdata <- t(scale(t(indata[rownames(annRow), rownames(annCol)])))
plotdata[plotdata > 3] <- 3 # 截断极端值
plotdata[plotdata < -3] <- -3 # 截断极端值

pheatmap(plotdata,
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = F,
         show_colnames = F,
         annotation_row = annRow,
         annotation_col = annCol,
         annotation_colors = annColors,
         color = colorRampPalette(c("#343493", "white", "#C24A45"))(64))


#outTab <- read.table("ouput_ICIsignatureGene.txt",sep = "\t",row.names = NULL,header = T,stringsAsFactors = F,check.names = F)
table(outTab$direct) # AB两种签名按相关性正负来分类，正相关标为A，否则标为B


set.seed(20221111)
dat.boruta <- as.data.frame(t(indata[outTab$gene,rownames(annCol)]))
borutafit <- Boruta(x = as.matrix(dat.boruta), 
                    y = as.factor(annCol$`Gene cluster`), # multiclassification
                    doTrace = 2,
                    maxRuns = 100,
                    ntree = 500)

boruta_fea <- attStats(borutafit)
boruta_fea <- rownames(boruta_fea[which(boruta_fea$decision == "Confirmed"),])
boruta.all <- outTab[which(outTab$gene %in% boruta_fea),]
table(boruta.all$direct)


# 提取A签名
set.seed(20221111)
dat.boruta <- as.data.frame(t(indata[outTab[which(outTab$direct == "A"),"gene"],rownames(annCol)]))
borutafit <- Boruta(x = as.matrix(dat.boruta), 
                    y = as.factor(annCol$`Gene cluster`), # multiclassification
                    doTrace = 2,
                    maxRuns = 100,
                    ntree = 500)

boruta_fea <- attStats(borutafit)
boruta_fea <- rownames(boruta_fea[which(boruta_fea$decision == "Confirmed"),])
boruta.A <- outTab[which(outTab$gene %in% boruta_fea),]
# 降维前后基因数量
ncol(dat.boruta)
#[1] 96
nrow(boruta.A)
#[1] 21
# 提取B签名
set.seed(20221111)
dat.boruta <- as.data.frame(t(indata[outTab[which(outTab$direct == "B"),"gene"],rownames(annCol)]))
borutafit <- Boruta(x = as.matrix(dat.boruta), 
                    y = as.factor(annCol$`Gene cluster`), # multiclassification
                    doTrace = 2,
                    maxRuns = 100,
                    ntree = 500)
boruta_fea <- attStats(borutafit)
boruta_fea <- rownames(boruta_fea[which(boruta_fea$decision == "Confirmed"),])
boruta.B <- outTab[which(outTab$gene %in% boruta_fea),]
# 降维前后基因数量
ncol(dat.boruta)
#[1] 91
nrow(boruta.B)
#[1] 35

##降维后绘图
annRow <- data.frame("ICI signature gene" = rep(c("A","B"), c(nrow(boruta.A), nrow(boruta.B))),
                     row.names = c(boruta.A$gene,boruta.B$gene),
                     check.names = F,
                     stringsAsFactors = F)
annColors <- list("Gene cluster" = c("A" = "#008ECB", "B" = "#EA921D"),
                  "Gender" = c("MALE" = "#79B789", "FEMALE" = "#B5262A"),
                  "Status" = c("Alive" = "#79B789", "Dead" = "#B5262A"),
                  "Grade" = c("T2" = "#53B1E7", "T3" = "#C78157", "T4" = "#A54D48"),
                  "Age" = colorRampPalette(c("#F7D202", "#96862A"))(64),
                  "ICI signature gene" = c("A" = "#D14039", "B" = "#008ECB"))
plotdata <- indata[rownames(annRow), rownames(annCol)]
plotdata[plotdata > 3] <- 3 # 截断极端值
plotdata[plotdata < -3] <- -3 # 截断极端值
pheatmap(plotdata,
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = F,
         show_colnames = F,
         annotation_row = annRow,
         annotation_col = annCol,
         annotation_colors = annColors,
         color = colorRampPalette(c("#343493", "white", "#C24A45"))(64))
dev.copy2pdf(file = "heatmap_dimension_reduced.pdf", width = 10, height = 8) # 保存图像

# 在“输入文件”那里，已经做过z-score处理，因此这里不再进行标准化
expr.A <- indata[boruta.A$gene,rownames(annCol)]
pca.A <- prcomp(t(expr.A), scale = F, center = F) # 如果数据没有标准化，这里都要设置为TRUE
pca1.A <- pca.A$x[,1] # 取出第一主成分

expr.B <- indata[boruta.B$gene,rownames(annCol)]
pca.B <- prcomp(t(expr.B), scale = F, center = F) #如果数据没有标准化，这里都要设置为TRUE
pca1.B <- pca.B$x[,1] # 取出第一主成分

ICI.score <- pca1.A - pca1.B # 主成份相减得到ICI得分
ICI.outtab <- data.frame(samID = rownames(annCol),
                         pca1.A = pca1.A[rownames(annCol)],
                         pca1.B = pca1.B[rownames(annCol)],
                         ICI.score = ICI.score[rownames(annCol)],
                         subtype = annCol$`Gene cluster`,
                         stringsAsFactors = F)
# 输出到文件
write.table(ICI.outtab,"~/Desktop/TCGA_work/TCGA_UVM_project/模型分数/output_ICI_score_TCGA.txt",sep = "\t",row.names = F,col.names = T,quote = F)


ggplot(data = ICI.outtab,aes(x = subtype, y = ICI.score, fill = subtype))+
  scale_fill_manual(values = c("#008ECB", "#EA921D")) + 
  geom_violin(alpha=0.4, position = position_dodge(width = .75),
              size=0.8, color="black") + # 边框线黑色
  geom_boxplot(notch = TRUE, outlier.size = -1, 
               color="black", lwd=0.8, alpha = 0.7)+ # 背景色透明化
  geom_point(shape = 21, size=2, 
             position = position_jitterdodge(), 
             color="black", alpha=1)+ # 边框线黑色
  theme_classic() +
  ylab(expression("ICI score")) +
  xlab("Subtype")  +
  theme(axis.ticks = element_line(size=0.2,color="black"),
        axis.ticks.length = unit(0.2,"cm"),
        legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)) +
  stat_compare_means(method = "kruskal.test", label.y = max(ICI.score))




##ROC time  
#三套数据 平滑ROC
##在ICI.outtab上进行绘制
table(rownames(ICI.outtab)==rownames(annCol))
annCol$riskscore <- ICI.outtab$ICI.score


library(survival)
library(survminer)
### 中位值划分
Type=ifelse(annCol[,'riskscore']<= median(annCol$riskscore), "Low", "High")
data=annCol
data$group=Type
data$group=factor(data$group, levels=c("Low", "High"))
diff=survdiff(Surv(futime, fustat) ~ group, data = data)
length=length(levels(factor(data[,"group"])))
pValue=1-pchisq(diff$chisq, df=length-1)

fit <- survfit(Surv(futime, fustat) ~ group, data = data)
bioCol=c("#0073C2","#EFC000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]

p=ggsurvplot(fit, 
             data=data,
             conf.int=F,
             pval=pValue,
             pval.size=6,
             legend.title='Risk_score',
             legend.labs=levels(factor(data[,"group"])),
             legend = c(0.88, 0.9),
             font.legend=12,
             xlab="Time(Days)",
             palette = bioCol,
             surv.median.line = "hv",
             risk.table=T,
             cumevents=F,
             risk.table.height=.25)

p    


rt9 <- data
class(rt9)
class(rt9$futime)
range(rt9$futime)

library(timeROC)
riskRoc <- timeROC(T = rt9$futime/365,delta = rt9$fustat,
                   marker = -rt9$riskscore,cause = 1,
                   weighting="marginal",
                   times = c(0.5,1,2,3,5))
multiTimeplot <- function(ROC,time,cex,xlab,ylab,title){
  library(ggsci)
  color <- pal_lancet()(length(time))
  plot(ROC$FP[,1], ROC$TP[,1], type="l", xlim=c(0,1), ylim=c(0,1), 
       col=color[1],
       xlab=xlab, 
       ylab=ylab,main=title)
  #如果直接plot roc对象，无法修改标题和坐标轴标签
  for(i in 2:length(time)){
    plot(ROC,time=time[i],add=T,col=color[i])
  }
  legend("bottomright",
         legend =paste("AUC at",time,"year:",round(ROC$AUC,digits = 4)),
         col = color,lwd = 1,
         bty = "n",cex = cex,text.col = color
  )
}

multiTimeplot(riskRoc,time = c(0.5,1,2,3,5),
              title="Time dependent ROC curve",
              xlab="False positive rate",
              ylab="True positive rate",
              cex=0.7)
riskRoc               #验证绘图结果




ROC <- timeROC(T = rt9$futime/365,   
               delta = rt9$fustat,   
               marker = -rt9$riskscore,   
               cause = 1,                
               weighting = "marginal",   
               times = c(1,2,3,5),       
               iid = TRUE)
ROC

df_plot <- data.frame(tpr = as.numeric(ROC$TP),
                      fpr = as.numeric(ROC$FP),
                      year = rep(c("1-year","2-year","3-year","5-year"),each = nrow(ROC$TP)))

head(df_plot)

library(ggplot2)

p <- ggplot(df_plot, aes(fpr, tpr, color = year)) +
  geom_smooth(se=FALSE, size=1.2)+ # 这就是平滑曲线的关键
  geom_abline(slope = 1, intercept = 0, color = "grey10",linetype = 2) +
  scale_color_manual(values = c("#E41A1C","#377EB8","#4DAF4A","#6A3D9AFF"),
                     name = NULL, 
                     labels = c(paste0("AUC at 1 year: ",round(ROC[["AUC"]][1],2)), 
                                paste0("AUC at 2 year: ",round(ROC[["AUC"]][2],2)), 
                                paste0("AUC at 3 year: ",round(ROC[["AUC"]][3],2)),
                                paste0("AUC at 5 year: ",round(ROC[["AUC"]][4],2)))
  ) + 
  coord_fixed(ratio = 1) +
  labs(x = "1 - Specificity", y = "Sensitivity") +
  theme_minimal(base_size = 14, base_family = "sans") +
  theme(legend.position = c(0.7,0.15), 
        panel.border = element_rect(fill = NA),
        axis.text = element_text(color = "black"))

p


##DCA
#rt9$Stage <- as.numeric(factor(str_split_fixed(KIRC_cli[rownames(rt9),"pathologic_stage"]," ",2)[,2]))

rt9$Stage <- as.numeric(factor(rt9$pstage,levels = c("T2","T3","T4")))
rt9$Age
rt9$cancer <- rt9$fustat==1
rt9$ttcancer <- rt9$futime/365 ##年份

rt9$riskscore
# 建立多个模型
Risk_score <- coxph(Surv(ttcancer, cancer) ~ riskscore + Stage, 
                    data = rt9)
Stage <- coxph(Surv(ttcancer, cancer) ~ Stage, data = rt9)
Age <- coxph(Surv(ttcancer, cancer) ~ Age, data = rt9)

df3 <- ggDCA::dca(Risk_score,Stage,Age,
                  times = c(1,2,3,5)
)
ggplot(df3,linetype = F)+
  scale_color_jama(name="Model Type")+
  theme_bw()+
  facet_wrap(~time)

##校准曲线

library(survival)
library(rms)
library(dplyr)
library(tidyr)
library(paletteer)
paletteer_d("RColorBrewer::Paired")

colnames(rt9)
rt9$futime
rt9$time <- rt9$futime
table(rt9$time>1835)

##1年
f1 <- cph(formula = Surv(time, fustat) ~  riskscore,data=rt9,
          x=T,y=T,surv = T,na.action=na.delete,time.inc = 365) 

#参数m=50表示每组50个样本进行重复计算
cal1 <- calibrate(f1, cmethod="KM", method="boot",u=365,m=25,B=1000) 

##2年
f2 <- cph(formula = Surv(time, fustat) ~  riskscore,data=rt9,
          x=T,y=T,surv = T,na.action=na.delete,time.inc = 730) 

#参数m=50表示每组50个样本进行重复计算
cal2 <- calibrate(f2, cmethod="KM", method="boot",u=730,m=25,B=1000) 

##3年
f3 <- cph(formula = Surv(time, fustat) ~  riskscore,data=rt9,
          x=T,y=T,surv = T,na.action=na.delete,time.inc = 1095) 

#参数m=50表示每组50个样本进行重复计算
cal3 <- calibrate(f3, cmethod="KM", method="boot",u=1095,m=25,B=1000) 


##5年
f5 <- cph(formula = Surv(time, fustat) ~  riskscore,data=rt9,
          x=T,y=T,surv = T,na.action=na.delete,time.inc = 1500) 

#参数m=50表示每组50个样本进行重复计算
cal5 <- calibrate(f5, cmethod="KM", method="boot",u=1825,m=25,B=1000) 

cal5 <- calibrate(f5, cmethod="KM", method="boot",u=1500) 
# #8年
# f8 <- cph(formula = Surv(time, fustat) ~  gpr_score,
#           data=rt9,x=T,y=T,surv = T,na.action=na.delete,time.inc = 2920)
# cal8 <- calibrate(f8, cmethod="KM", method="boot",u=2920,m=25,B=1000)


#pdf("calibration_compare.pdf",width = 8,height = 8)
plot(cal1,lwd = 2,lty = 0,errbar.col = c("#2166AC"),
     bty = "l", #只画左边和下边框
     xlim = c(0.2,1),ylim= c(0.2,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c("#2166AC"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal1[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#2166AC"), pch = 16)
mtext("")

plot(cal2,lwd = 2,lty = 0,errbar.col = c("#B2182B"),
     xlim = c(0.2,1),ylim= c(0.2,1),col = c("#B2182B"),add = T)
lines(cal2[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#B2182B"), pch = 16)

plot(cal3,lwd = 2,lty = 0,errbar.col = c("#FF7F00FF"),
     xlim = c(0.2,1),ylim= c(0.2,1),col = c("#FF7F00FF"),add = T)
lines(cal3[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#FF7F00FF"), pch = 16)

plot(cal5,lwd = 2,lty = 0,errbar.col = c("#6A3D9AFF"),
     xlim = c(0.2,1),ylim= c(0.2,1),col = c("#6A3D9AFF"),add = T)
lines(cal5[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#6A3D9AFF"), pch = 16)

abline(0,1, lwd = 2, lty = 3, col = c("#224444"))

legend("topleft", #图例的位置
       legend = c("1-year","2-year","3-year","5-year"), #图例文字
       col =c("#2166AC","#B2182B","#FF7F00FF","#6A3D9AFF"), #图例线的颜色，与文字对应
       lwd = 2,#图例中线的粗细
       cex = 1.2,#图例字体大小
       bty = "n")#不显示图例边框
#dev.off()


##time c-index  fig85


##其他两个队列进行验证分析

####GSE22138####

expr.A <- GSE22138_exp[boruta.A$gene,rownames(GSE22138_cli)]
pca.A <- prcomp(t(expr.A), scale = F, center = F) # 如果数据没有标准化，这里都要设置为TRUE
pca1.A <- pca.A$x[,1] # 取出第一主成分

expr.B <- GSE22138_exp[boruta.B$gene,rownames(GSE22138_cli)]
expr.B <- na.omit(expr.B)
pca.B <- prcomp(t(expr.B), scale = F, center = F) #如果数据没有标准化，这里都要设置为TRUE
pca1.B <- pca.B$x[,1] # 取出第一主成分

ICI.score <- pca1.A - pca1.B # 主成份相减得到ICI得分
# ICI.outtab <- data.frame(samID = rownames(annCol),
#                          pca1.A = pca1.A[rownames(annCol)],
#                          pca1.B = pca1.B[rownames(annCol)],
#                          ICI.score = ICI.score[rownames(annCol)],
#                          subtype = annCol$`Gene cluster`,
#                          stringsAsFactors = F)
# 输出到文件
table(rownames(GSE22138_cli)==names(ICI.score))
GSE22138_cli$riskscore <- ICI.score
#write.table(ICI.outtab,"output_ICI_score.txt",sep = "\t",row.names = F,col.names = T,quote = F)

library(survival)
library(survminer)
### 中位值划分
Type=ifelse(GSE22138_cli[,'riskscore']<= median(GSE22138_cli$riskscore), "Low", "High")
data=GSE22138_cli
data$group=Type
data$group=factor(data$group, levels=c("Low", "High"))
diff=survdiff(Surv(futime, fustat) ~ group, data = data)
length=length(levels(factor(data[,"group"])))
pValue=1-pchisq(diff$chisq, df=length-1)

fit <- survfit(Surv(futime, fustat) ~ group, data = data)
bioCol=c("#0073C2","#EFC000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]

p=ggsurvplot(fit, 
             data=data,
             conf.int=F,
             pval=pValue,
             pval.size=6,
             legend.title='Risk_score',
             legend.labs=levels(factor(data[,"group"])),
             legend = c(0.88, 0.9),
             font.legend=12,
             xlab="Time(Months)",
             palette = bioCol,
             surv.median.line = "hv",
             risk.table=T,
             cumevents=F,
             risk.table.height=.25)

p    


rt9 <- data
class(rt9)
class(rt9$futime)
range(rt9$futime)

library(timeROC)
riskRoc <- timeROC(T = rt9$futime*30/365,delta = rt9$fustat,
                   marker = -rt9$riskscore,cause = 1,
                   weighting="marginal",
                   times = c(0.5,1,2,3,5))
multiTimeplot <- function(ROC,time,cex,xlab,ylab,title){
  library(ggsci)
  color <- pal_lancet()(length(time))
  plot(ROC$FP[,1], ROC$TP[,1], type="l", xlim=c(0,1), ylim=c(0,1), 
       col=color[1],
       xlab=xlab, 
       ylab=ylab,main=title)
  #如果直接plot roc对象，无法修改标题和坐标轴标签
  for(i in 2:length(time)){
    plot(ROC,time=time[i],add=T,col=color[i])
  }
  legend("bottomright",
         legend =paste("AUC at",time,"year:",round(ROC$AUC,digits = 4)),
         col = color,lwd = 1,
         bty = "n",cex = cex,text.col = color
  )
}

multiTimeplot(riskRoc,time = c(0.5,1,2,3,5),
              title="Time dependent ROC curve",
              xlab="False positive rate",
              ylab="True positive rate",
              cex=0.7)
riskRoc               #验证绘图结果




ROC <- timeROC(T = rt9$futime*30/365,   
               delta = rt9$fustat,   
               marker = -rt9$riskscore,   
               cause = 1,                
               weighting = "marginal",   
               times = c(1,2,3,5),       
               iid = TRUE)
ROC

df_plot <- data.frame(tpr = as.numeric(ROC$TP),
                      fpr = as.numeric(ROC$FP),
                      year = rep(c("1-year","2-year","3-year","5-year"),each = nrow(ROC$TP)))

head(df_plot)

library(ggplot2)

p <- ggplot(df_plot, aes(fpr, tpr, color = year)) +
  geom_smooth(se=FALSE, size=1.2)+ # 这就是平滑曲线的关键
  geom_abline(slope = 1, intercept = 0, color = "grey10",linetype = 2) +
  scale_color_manual(values = c("#E41A1C","#377EB8","#4DAF4A","#6A3D9AFF"),
                     name = NULL, 
                     labels = c(paste0("AUC at 1 year: ",round(ROC[["AUC"]][1],2)), 
                                paste0("AUC at 2 year: ",round(ROC[["AUC"]][2],2)), 
                                paste0("AUC at 3 year: ",round(ROC[["AUC"]][3],2)),
                                paste0("AUC at 5 year: ",round(ROC[["AUC"]][4],2)))
  ) + 
  coord_fixed(ratio = 1) +
  labs(x = "1 - Specificity", y = "Sensitivity") +
  theme_minimal(base_size = 14, base_family = "sans") +
  theme(legend.position = c(0.7,0.15), 
        panel.border = element_rect(fill = NA),
        axis.text = element_text(color = "black"))

p


##DCA
#rt9$Stage <- as.numeric(factor(str_split_fixed(KIRC_cli[rownames(rt9),"pathologic_stage"]," ",2)[,2]))

rt9$Stage <- as.numeric(factor(rt9$pstage,levels = c("T2","T3","T4")))
rt9$Age
rt9$cancer <- rt9$fustat==1
rt9$ttcancer <- rt9$futime/365 ##年份

rt9$riskscore
# 建立多个模型
Risk_score <- coxph(Surv(ttcancer, cancer) ~ riskscore + Stage, 
                    data = rt9)
Stage <- coxph(Surv(ttcancer, cancer) ~ Stage, data = rt9)
Age <- coxph(Surv(ttcancer, cancer) ~ Age, data = rt9)

df3 <- ggDCA::dca(Risk_score,Stage,Age,
                  times = c(1,2,3,5)
)
ggplot(df3,linetype = F)+
  scale_color_jama(name="Model Type")+
  theme_bw()+
  facet_wrap(~time)

##校准曲线

library(survival)
library(rms)
library(dplyr)
library(tidyr)
library(paletteer)
paletteer_d("RColorBrewer::Paired")

colnames(rt9)
rt9$futime
rt9$time <- rt9$futime*30
table(rt9$time>1835)

##1年
f1 <- cph(formula = Surv(time, fustat) ~  riskscore,data=rt9,
          x=T,y=T,surv = T,na.action=na.delete,time.inc = 365) 

#参数m=50表示每组50个样本进行重复计算
cal1 <- calibrate(f1, cmethod="KM", method="boot",u=365,m=25,B=1000) 

##2年
f2 <- cph(formula = Surv(time, fustat) ~  riskscore,data=rt9,
          x=T,y=T,surv = T,na.action=na.delete,time.inc = 730) 

#参数m=50表示每组50个样本进行重复计算
cal2 <- calibrate(f2, cmethod="KM", method="boot",u=730,m=25,B=1000) 

##3年
f3 <- cph(formula = Surv(time, fustat) ~  riskscore,data=rt9,
          x=T,y=T,surv = T,na.action=na.delete,time.inc = 1095) 

#参数m=50表示每组50个样本进行重复计算
cal3 <- calibrate(f3, cmethod="KM", method="boot",u=1095,m=25,B=1000) 


##5年
f5 <- cph(formula = Surv(time, fustat) ~  riskscore,data=rt9,
          x=T,y=T,surv = T,na.action=na.delete,time.inc = 1825) 

#参数m=50表示每组50个样本进行重复计算
cal5 <- calibrate(f5, cmethod="KM", method="boot",u=1825,m=25,B=1000) 

#cal5 <- calibrate(f5, cmethod="KM", method="boot",u=1500) 
# #8年
# f8 <- cph(formula = Surv(time, fustat) ~  gpr_score,
#           data=rt9,x=T,y=T,surv = T,na.action=na.delete,time.inc = 2920)
# cal8 <- calibrate(f8, cmethod="KM", method="boot",u=2920,m=25,B=1000)


#pdf("calibration_compare.pdf",width = 8,height = 8)
plot(cal1,lwd = 2,lty = 0,errbar.col = c("#2166AC"),
     bty = "l", #只画左边和下边框
     xlim = c(0.2,1),ylim= c(0.2,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c("#2166AC"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal1[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#2166AC"), pch = 16)
mtext("")

plot(cal2,lwd = 2,lty = 0,errbar.col = c("#B2182B"),
     xlim = c(0.2,1),ylim= c(0.2,1),col = c("#B2182B"),add = T)
lines(cal2[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#B2182B"), pch = 16)

plot(cal3,lwd = 2,lty = 0,errbar.col = c("#FF7F00FF"),
     xlim = c(0.2,1),ylim= c(0.2,1),col = c("#FF7F00FF"),add = T)
lines(cal3[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#FF7F00FF"), pch = 16)

plot(cal5,lwd = 2,lty = 0,errbar.col = c("#6A3D9AFF"),
     xlim = c(0.2,1),ylim= c(0.2,1),col = c("#6A3D9AFF"),add = T)
lines(cal5[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#6A3D9AFF"), pch = 16)

abline(0,1, lwd = 2, lty = 3, col = c("#224444"))

legend("topleft", #图例的位置
       legend = c("1-year","2-year","3-year","5-year"), #图例文字
       col =c("#2166AC","#B2182B","#FF7F00FF","#6A3D9AFF"), #图例线的颜色，与文字对应
       lwd = 2,#图例中线的粗细
       cex = 1.2,#图例字体大小
       bty = "n")#不显示图例边框




###GSE84976#####
expr.A <- GSE84976_exp[boruta.A$gene,rownames(GSE84976_cli)]
pca.A <- prcomp(t(expr.A), scale = F, center = F) # 如果数据没有标准化，这里都要设置为TRUE
pca1.A <- pca.A$x[,1] # 取出第一主成分

expr.B <- GSE84976_exp[boruta.B$gene,rownames(GSE84976_cli)]
expr.B <- na.omit(expr.B)
pca.B <- prcomp(t(expr.B), scale = F, center = F) #如果数据没有标准化，这里都要设置为TRUE
pca1.B <- pca.B$x[,1] # 取出第一主成分

ICI.score <- -(pca1.A - pca1.B) # 主成份相减得到ICI得分
# ICI.outtab <- data.frame(samID = rownames(annCol),
#                          pca1.A = pca1.A[rownames(annCol)],
#                          pca1.B = pca1.B[rownames(annCol)],
#                          ICI.score = ICI.score[rownames(annCol)],
#                          subtype = annCol$`Gene cluster`,
#                          stringsAsFactors = F)
# 输出到文件
table(rownames(GSE84976_cli)==names(ICI.score))
GSE84976_cli$riskscore <- ICI.score
#write.table(ICI.outtab,"output_ICI_score.txt",sep = "\t",row.names = F,col.names = T,quote = F)

library(survival)
library(survminer)
### 中位值划分
Type=ifelse(GSE84976_cli[,'riskscore']<= median(GSE84976_cli$riskscore), "Low", "High")
data=GSE84976_cli
data$group=Type
data$group=factor(data$group, levels=c("Low", "High"))
diff=survdiff(Surv(futime, fustat) ~ group, data = data)
length=length(levels(factor(data[,"group"])))
pValue=1-pchisq(diff$chisq, df=length-1)

fit <- survfit(Surv(futime, fustat) ~ group, data = data)
bioCol=c("#0073C2","#EFC000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]

p=ggsurvplot(fit, 
             data=data,
             conf.int=F,
             pval=pValue,
             pval.size=6,
             legend.title='Risk_score',
             legend.labs=levels(factor(data[,"group"])),
             legend = c(0.88, 0.9),
             font.legend=12,
             xlab="Time(Months)",
             palette = bioCol,
             surv.median.line = "hv",
             risk.table=T,
             cumevents=F,
             risk.table.height=.25)

p    


rt9 <- data
class(rt9)
class(rt9$futime)
range(rt9$futime)

library(timeROC)
riskRoc <- timeROC(T = rt9$futime*30/365,delta = rt9$fustat,
                   marker = rt9$riskscore,cause = 1,
                   weighting="marginal",
                   times = c(0.5,1,2,3,5))
multiTimeplot <- function(ROC,time,cex,xlab,ylab,title){
  library(ggsci)
  color <- pal_lancet()(length(time))
  plot(ROC$FP[,1], ROC$TP[,1], type="l", xlim=c(0,1), ylim=c(0,1), 
       col=color[1],
       xlab=xlab, 
       ylab=ylab,main=title)
  #如果直接plot roc对象，无法修改标题和坐标轴标签
  for(i in 2:length(time)){
    plot(ROC,time=time[i],add=T,col=color[i])
  }
  legend("bottomright",
         legend =paste("AUC at",time,"year:",round(ROC$AUC,digits = 4)),
         col = color,lwd = 1,
         bty = "n",cex = cex,text.col = color
  )
}

multiTimeplot(riskRoc,time = c(0.5,1,2,3,5),
              title="Time dependent ROC curve",
              xlab="False positive rate",
              ylab="True positive rate",
              cex=0.7)
riskRoc               #验证绘图结果




ROC <- timeROC(T = rt9$futime*30/365,   
               delta = rt9$fustat,   
               marker = rt9$riskscore,   
               cause = 1,                
               weighting = "marginal",   
               times = c(2,3,5,10),       
               iid = TRUE)
ROC

df_plot <- data.frame(tpr = as.numeric(ROC$TP),
                      fpr = as.numeric(ROC$FP),
                      year = rep(c("2-year","3-year","5-year","10-year"),each = nrow(ROC$TP)))

head(df_plot)

library(ggplot2)

p <- ggplot(df_plot, aes(fpr, tpr, color = year)) +
  geom_smooth(se=FALSE, size=1.2)+ # 这就是平滑曲线的关键
  geom_abline(slope = 1, intercept = 0, color = "grey10",linetype = 2) +
  scale_color_manual(values = c("#E41A1C","#377EB8","#4DAF4A","#6A3D9AFF"),
                     name = NULL, 
                     labels = c(paste0("AUC at 2 year: ",round(ROC[["AUC"]][1],2)), 
                                paste0("AUC at 3 year: ",round(ROC[["AUC"]][2],2)), 
                                paste0("AUC at 5 year: ",round(ROC[["AUC"]][3],2)),
                                paste0("AUC at 10 year: ",round(ROC[["AUC"]][4],2)))
  ) + 
  coord_fixed(ratio = 1) +
  labs(x = "1 - Specificity", y = "Sensitivity") +
  theme_minimal(base_size = 14, base_family = "sans") +
  theme(legend.position = c(0.7,0.15), 
        panel.border = element_rect(fill = NA),
        axis.text = element_text(color = "black"))

p


##DCA
#rt9$Stage <- as.numeric(factor(str_split_fixed(KIRC_cli[rownames(rt9),"pathologic_stage"]," ",2)[,2]))

rt9$Stage <- as.numeric(factor(rt9$pstage,levels = c("T2","T3","T4")))
rt9$Age
rt9$cancer <- rt9$fustat==1
rt9$ttcancer <- rt9$futime/365 ##年份

rt9$riskscore
# 建立多个模型
Risk_score <- coxph(Surv(ttcancer, cancer) ~ riskscore + Stage, 
                    data = rt9)
Stage <- coxph(Surv(ttcancer, cancer) ~ Stage, data = rt9)
Age <- coxph(Surv(ttcancer, cancer) ~ Age, data = rt9)

df3 <- ggDCA::dca(Risk_score,Stage,Age,
                  times = c(1,2,3,5)
)
ggplot(df3,linetype = F)+
  scale_color_jama(name="Model Type")+
  theme_bw()+
  facet_wrap(~time)

##校准曲线

library(survival)
library(rms)
library(dplyr)
library(tidyr)
library(paletteer)
paletteer_d("RColorBrewer::Paired")

colnames(rt9)
rt9$futime
rt9$time <- rt9$futime*30
table(rt9$time>1835)

# ##1年
# f1 <- cph(formula = Surv(time, fustat) ~  riskscore,data=rt9,
#           x=T,y=T,surv = T,na.action=na.delete,time.inc = 365) 
# 
# #参数m=50表示每组50个样本进行重复计算
# cal1 <- calibrate(f1, cmethod="KM", method="boot",u=365,m=25,B=1000) 

##2年
f2 <- cph(formula = Surv(time, fustat) ~  riskscore,data=rt9,
          x=T,y=T,surv = T,na.action=na.delete,time.inc = 730) 

#参数m=50表示每组50个样本进行重复计算
cal2 <- calibrate(f2, cmethod="KM", method="boot",u=730,m=25,B=1000) 

##3年
f3 <- cph(formula = Surv(time, fustat) ~  riskscore,data=rt9,
          x=T,y=T,surv = T,na.action=na.delete,time.inc = 1095) 

#参数m=50表示每组50个样本进行重复计算
cal3 <- calibrate(f3, cmethod="KM", method="boot",u=1095,m=25,B=1000) 


##5年
f5 <- cph(formula = Surv(time, fustat) ~  riskscore,data=rt9,
          x=T,y=T,surv = T,na.action=na.delete,time.inc = 1825) 

#参数m=50表示每组50个样本进行重复计算
cal5 <- calibrate(f5, cmethod="KM", method="boot",u=1825,m=25,B=1000) 

#cal5 <- calibrate(f5, cmethod="KM", method="boot",u=1500) 
#8年
f8 <- cph(formula = Surv(time, fustat) ~  riskscore,
          data=rt9,x=T,y=T,surv = T,na.action=na.delete,time.inc = 3000)
cal8 <- calibrate(f8, cmethod="KM", method="boot",u=3000,m=25,B=1000)


#pdf("calibration_compare.pdf",width = 8,height = 8)
plot(cal1,lwd = 2,lty = 0,errbar.col = c("#2166AC"),
     bty = "l", #只画左边和下边框
     xlim = c(0.2,1),ylim= c(0.2,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c("#2166AC"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal1[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#2166AC"), pch = 16)
mtext("")

plot(cal2,lwd = 2,lty = 0,errbar.col = c("#B2182B"),
     xlim = c(0.2,1),ylim= c(0.2,1),col = c("#B2182B"),add = T)
lines(cal2[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#B2182B"), pch = 16)

plot(cal3,lwd = 2,lty = 0,errbar.col = c("#FF7F00FF"),
     xlim = c(0.2,1),ylim= c(0.2,1),col = c("#FF7F00FF"),add = T)
lines(cal3[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#FF7F00FF"), pch = 16)

plot(cal5,lwd = 2,lty = 0,errbar.col = c("#6A3D9AFF"),
     xlim = c(0.2,1),ylim= c(0.2,1),col = c("#6A3D9AFF"),add = T)
lines(cal5[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#6A3D9AFF"), pch = 16)

plot(cal8,lwd = 2,lty = 0,errbar.col = c("#6A3D9AFF"),
     xlim = c(0.2,1),ylim= c(0.2,1),col = c("#6A3D9AFF"),add = T)
lines(cal8[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#6A3D9AFF"), pch = 16)

abline(0,1, lwd = 2, lty = 3, col = c("#224444"))

legend("topleft", #图例的位置
       legend = c("2-year","3-year","5-year","10-year"), #图例文字
       col =c("#2166AC","#B2182B","#FF7F00FF","#6A3D9AFF"), #图例线的颜色，与文字对应
       lwd = 2,#图例中线的粗细
       cex = 1.2,#图例字体大小
       bty = "n")#不显示图例边框

