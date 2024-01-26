#save.image(file = "alldata.Rdata")
#load(file="alldata.Rdata")
setwd("/home/aim/Desktop/TCGA_work/TCGA_UVM_project/")

#load("/t8a/allimage_data/UVMall.rds")
#save.image(file = "/t8a/allimage_data/UVM_MOVICS_202208alldata.rds")

#load(file = "all.Rdata")
library(CMScaller)
library(MOVICS)
library(ggplot2)
library(devtools)
#devtools::install_local("/home/aim/Desktop/sc_sequnce/CMScaller-master.zip")
library(CMScaller)
###数据准备
load(file = "UVM_tcga.Rdata")
mo.data   <- UVM.tcga[1:5]# extract multi-omics data需要的多组学数据
count     <- UVM.tcga$count# extract raw count data for downstream analyses下游分析（差异分析）的count矩阵
fpkm      <- UVM.tcga$fpkm#extract fpkm data for downstream analyses下游分析的FPKM
maf       <- UVM.tcga$maf# extract maf for downstream analysis下游分析的maf文件
segment   <- UVM.tcga$segment#extract segmented copy number for downstream analyses下游分析的cnv数据
surv.info <- UVM.tcga$clin.info# extract survival information下游分析的生存数据

######开始聚类#####
# identify optimal clustering number (may take a while)
optk.uvm <- getClustNum(data        = mo.data,
                         is.binary   = c(F,F,F,F,T), # 注意第四个矩阵是一个二进制矩阵note: the 4th data is somatic mutation which is a binary matrix
                         try.N.clust = 2:8, # 推荐2:8 try cluster number from 2 to 8
                         fig.name    = "CLUSTER NUMBER OF TCGA-UVM")

optk.brca <- getClustNum(data        = mo.data,
                        is.binary   = c(F,F,F,F,T), # 注意第四个矩阵是一个二进制矩阵note: the 4th data is somatic mutation which is a binary matrix
                        try.N.clust = 2:8, # 推荐2:8 try cluster number from 2 to 8
                        fig.name    = "CLUSTER NUMBER OF TCGA-UVM")

# perform iClusterBayes (may take a while)使用getiClusterBayes，从单个算法中提取答案
iClusterBayes.res <- getiClusterBayes(data        = mo.data,
                                      N.clust     = 2,
                                      type        = c("gaussian","gaussian","gaussian","binomial"),
                                      n.burnin    = 1800,
                                      n.draw      = 1200,
                                      prior.gamma = c(0.5, 0.5, 0.5, 0.5),
                                      sdev        = 0.05,
                                      thin        = 3)
save(iClusterBayes.res,file = "iCluster.Rdata")
load(file = "iCluster.Rdata")
#> clustering done...
#> feature selection done...

###使用getMOIC函数，指定特定methodlist,也可以选择多种方法
iClusterBayes.res <- getMOIC(data        = mo.data,
                             N.clust     = 2,
                             methodslist = "iClusterBayes", # specify only ONE algorithm here
                             type        = c("gaussian","gaussian","gaussian","gaussian","binomial"), # data type corresponding to the list
                             n.burnin    = 1800,
                             n.draw      = 1200,
                             prior.gamma = c(0.5, 0.5, 0.5, 0.5),
                             sdev        = 0.05,
                             thin        = 3)
##一步到位 直接使用10中算法

# 使用九种算法 perform multi-omics integrative clustering with the rest of 10 algorithms
moic.res.list <- getMOIC(data        = mo.data,
                         methodslist = list("iClusterBayes","SNF", "PINSPlus", "NEMO", "COCA", "LRAcluster", "ConsensusClustering", "IntNMF", "CIMLR", "MoCluster"),
                         N.clust     = 2,
                         type        = c("gaussian", "gaussian","gaussian", "gaussian", "binomial"))

##先用
moic.res.list <- getMOIC(data        = mo.data,
                         methodslist = list("SNF", "PINSPlus", "NEMO", "COCA", "LRAcluster", "ConsensusClustering", "IntNMF", "CIMLR", "MoCluster"),
                         N.clust     = 2,
                         type        = c("gaussian", "gaussian","gaussian", "gaussian", "binomial"))
# 保存结果 attach iClusterBayes.res as a list using append() to moic.res.list with 9 results already
moic.res.list <- append(moic.res.list, 
                        list("iClusterBayes" = iClusterBayes.res))

# save moic.res.list to local path
save(moic.res.list, file = "moic.res.list_2types_10methods.rda")

load(file = "moic.res.list.rda")
##绘制共识矩阵 一个list和一张对称热图
cmoic.brca <- getConsensusMOIC(moic.res.list = moic.res.list,
                               fig.name      = "CONSENSUS HEATMAP",
                               distance      = "euclidean",
                               linkage       = "average")
##看簇相似性
getSilhouette(sil      = cmoic.brca$sil, # a sil object returned by getConsensusMOIC()
              fig.path = getwd(),
              fig.name = "SILHOUETTE",
              height   = 5.5,
              width    = 5)

# convert beta value to M value for stronger signal
indata <- mo.data
#indata$meth.beta <- log2(indata$meth.beta / (1 - indata$meth.beta))
#methtmp <- indata$meth.beta

# data normalization for heatmap四个组学一起绘制热图
plotdata <- getStdiz(data       = indata,
                     halfwidth  = c(2,2,2,2,NA), # no truncation for mutation
                     centerFlag = c(T,T,T,T,F), # no center for mutation
                     scaleFlag  = c(T,T,T,T,F)) # no scale for mutation
feat   <- iClusterBayes.res$feat.res
feat1  <- feat[which(feat$dataset == "mRNA.expr"),][1:10,"feature"] 
feat2  <- feat[which(feat$dataset == "lnc.expr"),][1:10,"feature"]
feat3  <- feat[which(feat$dataset == "mi.expr"),][1:10,"feature"]
feat4  <- feat[which(feat$dataset == "meth.beta"),][1:10,"feature"]
feat5  <- feat[which(feat$dataset == "mut.status"),][1:10,"feature"]
annRow <- list(feat1, feat2, feat3, feat4,feat5)
# set color for each omics data
# if no color list specified all subheatmaps will be unified to green and red color pattern
mRNA.col   <- c("#00FF00", "#008000", "#000000", "#800000", "#FF0000")
lncRNA.col <- c("#6699CC", "white"  , "#FF3C38")
mi.col <- c("#5AA3DAFF","#194155FF","#8A3F2AFF")
meth.col   <- c("#0074FE", "#96EBF9", "#FEE900", "#F00003")
mut.col    <- c("grey90" , "black")
col.list   <- list(mRNA.col, lncRNA.col,mi.col, meth.col, mut.col)

# comprehensive heatmap (may take a while)
getMoHeatmap(data          = plotdata,
             row.title     = c("mRNA","lncRNA","miRNA","Methylation","Mutation"),
             is.binary     = c(F,F,F,F,T), # the 4th data is mutation which is binary
             legend.name   = c("mRNA.FPKM","lncRNA.FPKM","miRNA.FPKM","M value","Mutated"),
             clust.res     = iClusterBayes.res$clust.res, # cluster results
             clust.dend    = NULL, # no dendrogram
             show.rownames = c(F,F,F,F,F), # specify for each omics data
             show.colnames = FALSE, # show no sample names
             annRow        = annRow, # mark selected features
             color         = col.list,
             annCol        = NULL, # no annotation for samples
             annColors     = NULL, # no annotation color
             width         = 10, # width of each subheatmap
             height        = 7, # height of each subheatmap
             fig.path = "./",
             fig.name      = "COMPREHENSIVE HEATMAP OF ICLUSTERBAYES")

# comprehensive heatmap (may take a while)这里悬着COCA聚类做热图
getMoHeatmap(data          = plotdata,
             row.title     = c("mRNA","lncRNA","Methylation","Mutation"),
             is.binary     = c(F,F,F,T), # the 4th data is mutation which is binary
             legend.name   = c("mRNA.FPKM","lncRNA.FPKM","M value","Mutated"),
             clust.res     = moic.res.list$COCA$clust.res, # cluster results
             clust.dend    = moic.res.list$COCA$clust.dend, # show dendrogram for samples
             color         = col.list,
             width         = 10, # width of each subheatmap
             height        = 7, # height of each subheatmap
             fig.name      = "COMPREHENSIVE HEATMAP OF COCA")

###映射其他临床信息增加列注释
# extract PAM50, pathologic stage and age for sample annotation
clin.info <- surv.info
table(is.na(clin.info$age))
mean(na.omit(clin.info$age))
clin.info[is.na(clin.info$age),"age"] <- 62
#gsub 更改AJCC
table(clin.info$AJCC)
AJCC <- clin.info$AJCC
AJCC <- gsub("IIA","II",AJCC)
AJCC <- gsub("IIB","II",AJCC)
AJCC <- gsub("IIIB","III",AJCC)
AJCC <- gsub("IIIA","III",AJCC)
AJCC <- gsub("IIIC","III",AJCC)
table(AJCC)
clin.info$AJCC <- AJCC

surv.info <- clin.info
annCol    <- surv.info[,c("AJCC", "pstage", "age"), drop = FALSE]
table(clin.info$pstage)
# generate corresponding colors for sample annotation
annColors <- list(age    = circlize::colorRamp2(breaks = c(min(annCol$age),
                                                           median(annCol$age),
                                                           max(annCol$age)), 
                                                colors = c("#0000AA", "#555555", "#AAAA00")),
                  AJCC  = c(#"I" = "blue",
                            "II"   = "red",
                            "III"   = "yellow",
                            "IV"   = "green",
                            "X"="black"),
                  
                  pstage = c(
                    #"T1"    = "green",
                    "T2"    = "blue",
                    "T3"    = "red",
                    "T4"    = "yellow"
                    #"TX"="black
                    ))


# comprehensive heatmap (may take a while)

getMoHeatmap(data          = plotdata,
             row.title     = c("mRNA","lncRNA","miRNA","Methylation","Mutation"),
             is.binary     = c(F,F,F,F,T), # the 4th data is mutation which is binary
             legend.name   = c("mRNA.FPKM","lncRNA.FPKM","miRNA.FPKM","M value","Mutated"),
             clust.res     = cmoic.brca$clust.res, # consensusMOIC results
             clust.dend    = NULL, # show no dendrogram for samples
             show.rownames = c(F,F,F,F,F), # specify for each omics data
             show.colnames = FALSE, # show no sample names
             show.row.dend = c(F,F,F,F,F), # show no dendrogram for features
             #annRow        = NULL, # no selected features
             annRow        = annRow,
             color         = col.list,
             annCol        = annCol, # annotation for samples
             annColors     = annColors, # annotation color
             width         = 10, # width of each subheatmap
             height        = 7, # height of each subheatmap
             fig.name      = "clincial_COMPREHENSIVE HEATMAP OF CONSENSUSMOIC")



####探究亚型绘制图片#####
####生存比较
# survival comparison
cmoic.UVM <- cmoic.brca
surv.info$futime

surv.brca <- compSurv(moic.res         = cmoic.UVM,
                      surv.info        = surv.info,
                      convt.time       = "m", # convert day unit to month
                      surv.median.line = "h", # draw horizontal line at median survival
                      xyrs.est         = c(3,5), # estimate 5 and 10-year survival
                      fig.name         = "KAPLAN-MEIER CURVE OF CONSENSUSMOIC")
#> --a total of 643 samples are identified.
#> --removed missing values.
#> --leaving 642 observations.
#> --cut survival curve up to 10 years.

###临床信息比较
pseudo.moic.res                 <- list("clust.res" = surv.info,
                                        "mo.method" = "PAM50")

# make pseudo samID
pseudo.moic.res$clust.res$samID <- rownames(pseudo.moic.res$clust.res)
clin.brca <- compClinvar(moic.res      = pseudo.moic.res,
                         var2comp      = surv.info, # data.frame needs to summarize (must has row names of samples)
                         strata        = "Subtype", # stratifying variable (e.g., Subtype in this example)
                         factorVars    = c("AJCC","pstage","fustat"), # features that are considered categorical variables
                         nonnormalVars = "futime", # feature(s) that are considered using nonparametric test
                         exactVars     = "pstage", # feature(s) that are considered using exact test
                         doWord        = TRUE, # generate .docx file in local path
                         tab.name      = "SUMMARIZATION OF CLINICAL FEATURES")
#> --all samples matched.
print(clin.brca$compTab)

###突变频率比较
# mutational frequency comparison
mut.brca <- compMut(moic.res     = cmoic.brca,
                    mut.matrix   = mo.data$mut.status, # binary somatic mutation matrix
                    doWord       = TRUE, # generate table in .docx format
                    doPlot       = TRUE, # draw OncoPrint
                    freq.cutoff  = 0.05, # keep those genes that mutated in at least 5% of samples
                    p.adj.cutoff = 0.25, # keep those genes with adjusted p value < 0.05 to draw OncoPrint
                    innerclust   = TRUE, # perform clustering within each subtype
                    annCol       = annCol, # same annotation for heatmap
                    annColors    = annColors, # same annotation color for heatmap
                    width        = 6, 
                    height       = 2,
                    fig.name     = "ONCOPRINT FOR SIGNIFICANT MUTATIONS",
                    tab.name     = "INDEPENDENT TEST BETWEEN SUBTYPE AND MUTATION")
#> --all samples matched.
print(mut.brca)
###肿瘤突变负荷比较
head(maf)
# compare TMB
tmb.brca <- compTMB(moic.res     = cmoic.brca,
                    maf          = maf,
                    rmDup        = TRUE, # remove duplicated variants per sample
                    rmFLAGS      = FALSE, # keep FLAGS mutations
                    exome.size   = 38, # estimated exome size
                    test.method  = "nonparametric", # statistical testing method
                    width        = 6, 
                    height       = 6,
                    fig.name     = "DISTRIBUTION OF TMB AND TITV")
head(tmb.brca$TMB.dat)
###基于拷贝数变异计算
# change column names of segment data
colnames(segment) <- c("sample","chrom","start","end","value")
head(segment)
segment <- as.data.frame(segment)
# compare FGA, FGG, and FGL
fga.brca <- compFGA(moic.res     = cmoic.brca,
                    segment      = segment,
                    iscopynumber = FALSE, # this is a segmented copy number file
                    cnathreshold = 0.2, # threshold to determine CNA gain or loss
                    test.method  = "nonparametric", # statistical testing method
                    fig.name     = "BARPLOT OF FGA")
fga.brca
head(fga.brca$summary)
####比较药物敏感性
# drug sensitivity comparison
colnames(fpkms) <- str_sub(colnames(fpkms),1,15)
fpkm <-UVM_mRNA_fpkm

drug.brca <- compDrugsen(moic.res    = cmoic.brca,
                         norm.expr   = fpkm[,cmoic.brca$clust.res$samID], # double guarantee sample order
                         drugs       = c("Sunitinib", "Pembrolizumab"), # a vector of names of drug in GDSC
                         tissueType  = "skin", # choose specific tissue type to construct ridge regression model
                         test.method = "nonparametric", # statistical testing method
                         prefix      = "BOXVIOLIN OF ESTIMATED IC50") 
head(drug.brca$Sunitinib)

###使用最全版本
drug.brca <- compDrugsen(moic.res    = cmoic.brca,
                         norm.expr   = fpkm[,cmoic.brca$clust.res$samID], # double guarantee sample order
                         drugs       = c("Saracatinib","Crizotinib","Axitinib","Sunitinib","Erlotinib","Pazopanib","Temsirolimus"), # a vector of names of drug in GDSC
                         tissueType  = "skin", # choose specific tissue type to construct ridge regression model
                         test.method = "nonparametric", # statistical testing method
                         prefix      = "BOXVIOLIN OF ESTIMATED  IC50 ")



drug.brca <- compDrugsen(moic.res    = cmoic.brca,
                         norm.expr   = fpkm[,cmoic.brca$clust.res$samID], # double guarantee sample order
                         drugs       = c("A.443654", "A.770041", "ABT.263", "ABT.888", "AG.014699", 
                                         "AICAR", "AKT.inhibitor.VIII", "AMG.706"," AP.24534", 
                                         "AS601245","ATRA", "AUY922", "Axitinib", "AZ628", "AZD.0530", 
                                         "AZD.2281", "AZD6244", "AZD6482", "AZD7762", "AZD8055", "BAY.61.3606", 
                                         "Bexarotene", "BI.2536", "BIBW2992", "Bicalutamide", "BI.D1870", 
                                         "BIRB.0796", "Bleomycin", "BMS.509744", "BMS.536924", "BMS.708163", 
                                         "BMS.754807", "Bortezomib", "Bosutinib", "Bryostatin.1", "BX.795", 
                                         "Camptothecin", "CCT007093", "CCT018159", "CEP.701", "CGP.082996", 
                                         "CGP.60474"," CHIR.99021", "CI.1040", "Cisplatin", "CMK", "Cyclopamine", 
                                         "Cytarabine", "Dasatinib", "DMOG", "Docetaxel", "Doxorubicin", "EHT.1864", 
                                         "Elesclomol", "Embelin", "Epothilone.B", "Erlotinib", "Etoposide", "FH535", 
                                         "FTI.277", "GDC.0449", "GDC0941", "Gefitinib", "Gemcitabine", "GNF.2", "GSK269962A", 
                                         "GSK.650394", "GW.441756", "GW843682X", "Imatinib", "IPA.3", "JNJ.26854165", "JNK.9L", 
                                         "JNK.Inhibitor.VIII", "JW.7.52.1", "KIN001.135", "KU.55933", "Lapatinib", "Lenalidomide", 
                                         "LFM.A13", "Metformin", "Methotrexate", "MG.132"," Midostaurin", "Mitomycin.C", "MK.2206", 
                                         "MS.275", "Nilotinib", "NSC.87877", "NU.7441", "Nutlin.3a", "NVP.BEZ235", "NVP.TAE684", 
                                         "Obatoclax.Mesylate", "OSI.906", "PAC.1", "Paclitaxel", "Parthenolide", "Pazopanib", 
                                         "PD.0325901", "PD.0332991", "PD.173074", "PF.02341066", "PF.4708671", "PF.562271", "PHA.665752", 
                                         "PLX4720", "Pyrimethamine", "QS11", "Rapamycin", "RDEA119", "RO.3306", "Roscovitine", "Salubrinal", 
                                         "SB.216763", "SB590885", "Shikonin", "SL.0101.1", "Sorafenib", "S.Trityl.L.cysteine", "Sunitinib", 
                                         "Temsirolimus", "Thapsigargin", "Tipifarnib", "TW.37", "Vinblastine", "Vinorelbine", "Vorinostat", 
                                         "VX.680", "VX.702", "WH.4.023", "WO2009093972", "WZ.1.84", "X17.AAG", "X681640", "XMD8.85", 
                                         "Z.LLNle.CHO", "ZM.447439"), # a vector of names of drug in GDSC
                         tissueType  = "all", # choose specific tissue type to construct ridge regression model
                         test.method = "nonparametric", # statistical testing method
                         prefix      = "BOXVIOLIN OF ESTIMATED  IC50 ")


###皮肤癌专属药物组
drug.brca <- compDrugsen(moic.res    = cmoic.brca,
                         norm.expr   = fpkm[,cmoic.brca$clust.res$samID], # double guarantee sample order
                         drugs       = c("Erlotinib", "Rapamycin", "Sunitinib", "PHA.665752", "MG.132", "Paclitaxel", 
                                         "Cyclopamine", "AZ628", "Sorafenib", 
                                         "VX.680", "Imatinib", "TAE684", "Crizotinib", "Saracatinib", "S.Trityl.L.cysteine", 
                                         "Z.LLNle.CHO", "Dasatinib", 
                                         "GNF.2", "CGP.60474", "CGP.082996", "A.770041", "WH.4.023", "WZ.1.84", 
                                         "BI.2536", "BMS.536924", "BMS.509744", "CMK", 
                                         "Pyrimethamine", "JW.7.52.1", "A.443654", "GW843682X", 
                                         "MS.275","Parthenolide", "KIN001.135", "TGX221", "Bortezomib", 
                                         "XMD8.85", "Roscovitine", "Salubrinal", "Lapatinib", "GSK269962A", "Doxorubicin", 
                                         "Etoposide", "Gemcitabine", "Mitomycin.C", 
                                         "Vinorelbine", "NSC.87877", "Bicalutamide"," QS11", "CP466722", "Midostaurin", 
                                         "CHIR.99021", "AP.24534", "AZD6482", "JNK.9L", 
                                         "PF.562271","HG.6.64.1", "JQ1", "JQ12", "DMOG", "FTI.277", "OSU.03012", "Shikonin", 
                                         "AKT.inhibitor.VIII", "Embelin", "FH535", 
                                         "PAC.1", "IPA.3", "GSK.650394", "BAY.61.3606.5","Fluorouracil", "Thapsigargin",
                                         "Obatoclax.Mesylate", "BMS.754807", "Lisitinib", "Bexarotene", "Bleomycin", 
                                         "LFM.A13", "GW.2580", 
                                         "AUY922", "Phenformin", "Bryostatin" ), # a vector of names of drug in GDSC
                         tissueType  = "nervous_system", # choose specific tissue type to construct ridge regression model
                         test.method = "nonparametric", # statistical testing method
                         prefix      = "BOXVIOLIN OF ESTIMATED  IC50 ")

###和传统分类比较
# customize the factor level for pstage
table(surv.info$pstage)

surv.info$pstage <- factor(surv.info$pstage, levels = c("T2","T3","T4"))

# agreement comparison (support up to 6 classifications include current subtype)
agree.brca <- compAgree(moic.res  = cmoic.brca,
                        subt2comp = surv.info[,c("AJCC","pstage")],
                        doPlot    = TRUE,
                        box.width = 0.2,
                        fig.name  = "AGREEMENT OF CONSENSUSMOIC WITH PAM50 AND PSTAGE")
#> --all samples matched.
print(agree.brca)



####亚型之间比较#####
##差异基因
# run DEA with edgeR
runDEA(dea.method = "edger",
       expr       = count, # raw count data
       moic.res   = cmoic.brca,
       prefix     = "TCGA-UVM") # prefix of figure name
colnames(count)
count_int <- ceiling(count)
runDEA(dea.method = "deseq2",
       expr       = count_int,
       moic.res   = cmoic.brca,
       prefix     = "TCGA-UVM")

runDEA(dea.method = "limma",
       expr       = fpkm, # normalized expression data
       moic.res   = cmoic.brca,
       prefix     = "TCGA-UVM")

###每组的marker基因
UVM_mRNA_fpkm <- fpkm
# choose edgeR result to identify subtype-specific up-regulated biomarkers
# choose edgeR result to identify subtype-specific up-regulated biomarkers
marker.up <- runMarker(moic.res      = cmoic.brca,
                       dea.method    = "edger", # name of DEA method
                       prefix        = "TCGA-UVM", # MUST be the same of argument in runDEA()
                       dat.path      = getwd(), # path of DEA files
                       res.path      = getwd(), # path to save marker files
                       p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                       p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                       dirct         = "up", # direction of dysregulation in expression
                       n.marker      = 100, # number of biomarkers for each subtype
                       doplot        = TRUE, # generate diagonal heatmap
                       norm.expr     = UVM_mRNA_fpkm, # use normalized expression as heatmap input
                       annCol        = annCol, # sample annotation in heatmap
                       annColors     = annColors, # colors for sample annotation
                       show_rownames = FALSE, # show no rownames (biomarker name)
                       fig.name      = "UPREGULATED BIOMARKER HEATMAP")
#> --all samples matched.
#> --log2 transformation done for expression data.
head(marker.up$templates)

# choose edgeR result to identify subtype-specific up-regulated biomarkers
marker.down <- runMarker(moic.res      = cmoic.brca,
                         dea.method    = "edger", # name of DEA method
                         prefix        = "TCGA-UVM", # MUST be the same of argument in runDEA()
                         dat.path      = getwd(), # path of DEA files
                         res.path      = getwd(), # path to save marker files
                         p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                         p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                         dirct         = "down", # direction of dysregulation in expression
                         n.marker      = 100, # number of biomarkers for each subtype
                         doplot        = TRUE, # generate diagonal heatmap
                         norm.expr     = UVM_mRNA_fpkm, # use normalized expression as heatmap input
                         annCol        = annCol, # sample annotation in heatmap
                         annColors     = annColors, # colors for sample annotation
                         show_rownames = FALSE, # show no rownames (biomarker name)
                         fig.name      = "UPREGULATED BIOMARKER HEATMAP")
#> --all samples matched.
#> --log2 transformation done for expression data.
head(marker.down$templates)
####GSEA分析
# MUST locate ABSOLUTE path of msigdb file
MSIGDB.FILE <- system.file("extdata", "c5.bp.v7.1.symbols.xls", package = "MOVICS", mustWork = TRUE)
# run GSEA to identify up-regulated GO pathways using results from edgeR
gsea.up <- runGSEA(moic.res     = cmoic.brca,
                   dea.method   = "deseq2", # name of DEA method
                   prefix       = "TCGA-UVM", # MUST be the same of argument in runDEA()
                   dat.path     = getwd(), # path of DEA files
                   res.path     = getwd(), # path to save GSEA files
                   msigdb.path  = MSIGDB.FILE, # MUST be the ABSOLUTE path of msigdb file
                   norm.expr    = UVM_mRNA_fpkm, # use normalized expression to calculate enrichment score
                   dirct        = "up", # direction of dysregulation in pathway
                   p.cutoff     = 0.05, # p cutoff to identify significant pathways
                   p.adj.cutoff = 0.25, # padj cutoff to identify significant pathways
                   gsva.method  = "gsva", # method to calculate single sample enrichment score
                   norm.method  = "mean", # normalization method to calculate subtype-specific enrichment score
                   fig.name     = "UPREGULATED PATHWAY HEATMAP")

print(gsea.up$gsea.list$CS1[1:6,3:6])##查看CS1的GSEA
head(round(gsea.up$grouped.es,3))##查看五组的GSEA



# run GSEA to identify down-regulated GO pathways using results from DESeq2
gsea.dn <- runGSEA(moic.res     = cmoic.brca,
                   dea.method   = "deseq2",
                   prefix       = "TCGA-UVM",
                   dat.path     = getwd(), # path of DEA files
                   res.path     = getwd(), # path to save GSEA files
                   msigdb.path  = MSIGDB.FILE,
                   norm.expr    = UVM_mRNA_fpkm,
                   dirct        = "down",
                   p.cutoff     = 0.05,
                   p.adj.cutoff = 0.25,
                   gsva.method  = "gsva", # switch to ssgsea
                   norm.method  = "median", # switch to median
                   fig.name     = "DOWNREGULATED PATHWAY HEATMAP") 

##cc分析

# run GSEA to identify up-regulated GO pathways using results from edgeR
gsea.up <- runGSEA(moic.res     = cmoic.brca,
                   dea.method   = "deseq2", # name of DEA method
                   prefix       = "TCGA-UVM", # MUST be the same of argument in runDEA()
                   dat.path     = getwd(), # path of DEA files
                   res.path     = getwd(), # path to save GSEA files
                   msigdb.path  = "../c5.go.cc.v7.2.symbols.gmt", # MUST be the ABSOLUTE path of msigdb file
                   norm.expr    = UVM_mRNA_fpkm, # use normalized expression to calculate enrichment score
                   dirct        = "up", # direction of dysregulation in pathway
                   p.cutoff     = 0.05, # p cutoff to identify significant pathways
                   p.adj.cutoff = 0.25, # padj cutoff to identify significant pathways
                   gsva.method  = "gsva", # method to calculate single sample enrichment score
                   norm.method  = "mean", # normalization method to calculate subtype-specific enrichment score
                   fig.name     = "UPREGULATED GO-CC HEATMAP")

print(gsea.up$gsea.list$CS1[1:6,3:6])##查看CS1的GSEA
head(round(gsea.up$grouped.es,3))##查看五组的GSEA

# run GSEA to identify down-regulated GO pathways using results from DESeq2
gsea.dn <- runGSEA(moic.res     = cmoic.brca,
                   dea.method   = "deseq2",
                   prefix       = "TCGA-UVM",
                   dat.path     = getwd(), # path of DEA files
                   res.path     = getwd(), # path to save GSEA files
                   msigdb.path  = "../c5.go.cc.v7.2.symbols.gmt",
                   norm.expr    = UVM_mRNA_fpkm,
                   dirct        = "down",
                   p.cutoff     = 0.05,
                   p.adj.cutoff = 0.25,
                   gsva.method  = "gsva", # switch to ssgsea
                   norm.method  = "mean", # switch to median
                   fig.name     = "DOWNREGULATED GO-CC HEATMAP") 


##MF分析

# run GSEA to identify up-regulated GO pathways using results from edgeR
gsea.up <- runGSEA(moic.res     = cmoic.brca,
                   dea.method   = "deseq2", # name of DEA method
                   prefix       = "TCGA-UVM", # MUST be the same of argument in runDEA()
                   dat.path     = getwd(), # path of DEA files
                   res.path     = getwd(), # path to save GSEA files
                   msigdb.path  = "../c5.go.mf.v7.2.symbols.gmt", # MUST be the ABSOLUTE path of msigdb file
                   norm.expr    = UVM_mRNA_fpkm, # use normalized expression to calculate enrichment score
                   dirct        = "up", # direction of dysregulation in pathway
                   p.cutoff     = 0.05, # p cutoff to identify significant pathways
                   p.adj.cutoff = 0.25, # padj cutoff to identify significant pathways
                   gsva.method  = "gsva", # method to calculate single sample enrichment score
                   norm.method  = "mean", # normalization method to calculate subtype-specific enrichment score
                   fig.name     = "UPREGULATED GO-MF HEATMAP")

print(gsea.up$gsea.list$CS1[1:6,3:6])##查看CS1的GSEA
head(round(gsea.up$grouped.es,3))##查看五组的GSEA

# run GSEA to identify down-regulated GO pathways using results from DESeq2
gsea.dn <- runGSEA(moic.res     = cmoic.brca,
                   dea.method   = "deseq2",
                   prefix       = "TCGA-UVM",
                   dat.path     = getwd(), # path of DEA files
                   res.path     = getwd(), # path to save GSEA files
                   msigdb.path  = "../c5.go.mf.v7.2.symbols.gmt",
                   norm.expr    = UVM_mRNA_fpkm,
                   dirct        = "down",
                   p.cutoff     = 0.05,
                   p.adj.cutoff = 0.25,
                   gsva.method  = "gsva", # switch to ssgsea
                   norm.method  = "mean", # switch to median
                   fig.name     = "DOWNREGULATED GO-MF HEATMAP") 

##GSVA这里作者内置数据集其实是免疫浸润相关的
# MUST locate ABSOLUTE path of gene set file
GSET.FILE <- 
  system.file("extdata", "gene sets of interest.gmt", package = "MOVICS", mustWork = TRUE)
# run GSVA to estimate single sample enrichment score based on given gene set of interest
gsva.res <- 
  runGSVA(moic.res      = cmoic.brca,
          norm.expr     = UVM_mRNA_fpkm,
          gset.gmt.path = GSET.FILE, # ABSOLUTE path of gene set file
          gsva.method   = "gsva", # method to calculate single sample enrichment score
          annCol        = annCol,
          annColors     = annColors,
          fig.path      = getwd(),
          fig.name      = "GENE SETS OF INTEREST HEATMAP",
          height        = 5,
          width         = 8)
print(gsva.res$raw.es[1:3,1:3])


# run GSVA to estimate single sample enrichment score based on given gene set of interest
gsva.res2 <- 
  runGSVA(moic.res      = cmoic.brca,
          norm.expr     = UVM_mRNA_fpkm,
          gset.gmt.path = "./c7.all.v7.2.symbols.gmt", # ABSOLUTE path of gene set file
          gsva.method   = "gsva", # method to calculate single sample enrichment score
          annCol        = annCol,
          annColors     = annColors,
          show_rownames = F,
          fig.path      = getwd(),
          fig.name      = "GENE SETS OF IMMUNE HEATMAP",
          height        = 5,
          width         = 8)
print(gsva.res$raw.es[1:3,1:3])



##代谢分数
GSET.FILE2 <- "/t8a/GEO_work/gmtdata/try.gmt"
gsva.res <- 
  runGSVA(moic.res      = cmoic.brca,
          norm.expr     = UVM_mRNA_fpkm,
          gset.gmt.path = GSET.FILE2, # ABSOLUTE path of gene set file
          gsva.method   = "gsva", # method to calculate single sample enrichment score
          annCol        = annCol,
          annColors     = annColors,
          fig.path      = getwd(),
          show_rownames = T,
          fig.name      = "GENE SETS OF metabolism HEATMAP",
          height        = 20,
          width         = 8)


##TME 免疫相关富集分数
GSET.FILE3 <- "/t8a/GEO_work/gmtdata/tme.gmt"
gsva.res <- 
  runGSVA(moic.res      = cmoic.brca,
          norm.expr     = UVM_mRNA_fpkm,
          gset.gmt.path = GSET.FILE3, # ABSOLUTE path of gene set file
          gsva.method   = "gsva", # method to calculate single sample enrichment score
          annCol        = annCol,
          annColors     = annColors,
          fig.path      = getwd(),
          show_rownames = T,
          fig.name      = "GENE SETS OF TME HEATMAP",
          height        = 20,
          width         = 8)

###进行差异分析

###肿瘤相关富集分数
GSET.FILE4 <- "/t8a/GEO_work/gmtdata/tumor.gmt"
gsva.res <- 
  runGSVA(moic.res      = cmoic.brca,
          norm.expr     = UVM_mRNA_fpkm,
          gset.gmt.path = GSET.FILE4, # ABSOLUTE path of gene set file
          gsva.method   = "gsva", # method to calculate single sample enrichment score
          annCol        = annCol,
          annColors     = annColors,
          fig.path      = getwd(),
          show_rownames = T,
          fig.name      = "GENE SETS OF tumor HEATMAP",
          height        = 5,
          width         = 8)


####截止在这里
save.image(file = "all.Rdata")

###提取分组信息 用于其他组学维度差异分析
###拿外部数据验证
##使用GSE2748 ##cancer research文章 带有预后信息   28个样本
load("/t8a/GEO/UVM_validation/GSE2748input.rds")


# run NTP in Yau cohort by using up-regulated biomarkers
yau.ntp.pred <- runNTP(expr       = brca.yau$mRNA.expr,
                       templates  = marker.up$templates, # the template has been already prepared in runMarker()
                       scaleFlag  = TRUE, # scale input data (by default)
                       centerFlag = TRUE, # center input data (by default)
                       doPlot     = TRUE, # to generate heatmap
                       fig.name   = "NTP HEATMAP FOR YAU") 
head(yau.ntp.pred$ntp.res)
##外部数据生存分析
# compare survival outcome in Yau cohort
surv.yau <- compSurv(moic.res         = yau.ntp.pred,
                     surv.info        = brca.yau$clin.info,
                     convt.time       = "m", # switch to month
                     surv.median.line = "hv", # switch to both
                     fig.name         = "KAPLAN-MEIER CURVE OF NTP FOR YAU") 
#> --a total of 682 samples are identified.
#> --cut survival curve up to 10 years.
print(surv.yau)
# compare agreement in Yau cohort
agree.yau <- compAgree(moic.res  = yau.ntp.pred,
                       subt2comp = brca.yau$clin.info[, "PAM50", drop = FALSE],
                       doPlot    = TRUE,
                       fig.name  = "YAU PREDICTEDMOIC WITH PAM50")
#> --all samples matched.
print(agree.yau)

###预测相似性
yau.pam.pred <- runPAM(train.expr  = fpkm,
                       moic.res    = cmoic.brca,
                       test.expr   = brca.yau$mRNA.expr)
print(yau.pam.pred$IGP)

###其他相似性
# predict subtype in discovery cohort using NTP
tcga.ntp.pred <- runNTP(expr      = fpkm,
                        templates = marker.up$templates,
                        doPlot    = FALSE) 
# predict subtype in discovery cohort using PAM
tcga.pam.pred <- runPAM(train.expr  = fpkm,
                        moic.res    = cmoic.brca,
                        test.expr   = fpkm)
# check consistency between current and NTP-predicted subtype in discovery TCGA-BRCA
runKappa(subt1     = cmoic.brca$clust.res$clust,
         subt2     = tcga.ntp.pred$clust.res$clust,
         subt1.lab = "CMOIC",
         subt2.lab = "NTP",
         fig.name  = "CONSISTENCY HEATMAP FOR TCGA between CMOIC and NTP")
# check consistency between current and PAM-predicted subtype in discovery TCGA-BRCA
runKappa(subt1     = cmoic.brca$clust.res$clust,
         subt2     = tcga.pam.pred$clust.res$clust,
         subt1.lab = "CMOIC",
         subt2.lab = "PAM",
         fig.name  = "CONSISTENCY HEATMAP FOR TCGA between CMOIC and PAM")
# check consistency between NTP and PAM-predicted subtype in validation Yau-BRCA
runKappa(subt1     = yau.ntp.pred$clust.res$clust,
         subt2     = yau.pam.pred$clust.res$clust,
         subt1.lab = "NTP",
         subt2.lab = "PAM",
         fig.name  = "CONSISTENCY HEATMAP FOR YAU")




###一些小技巧####
# include original clinical information as `clust.res` and a string value for `mo.method` to a list
pseudo.moic.res                 <- list("clust.res" = surv.info,
                                        "mo.method" = "PAM50")

# make pseudo samID
pseudo.moic.res$clust.res$samID <- rownames(pseudo.moic.res$clust.res)

# make pseudo clust using a mapping relationship
pseudo.moic.res$clust.res$clust <- sapply(pseudo.moic.res$clust.res$PAM50,
                                          switch,
                                          "Basal"   = 1, # relabel Basal as 1
                                          "Her2"    = 2, # relabel Her2 as 2
                                          "LumA"    = 3, # relabel LumA as 3
                                          "LumB"    = 4, # relabel LumnB as 4
                                          "Normal"  = 5) # relabel Normal as 5
pseudo.moic.res$clust.res$clust <- cmoic.UVM$clust.res$clust
head(pseudo.moic.res$clust.res)
# survival comparison
pam50.brca <- compSurv(moic.res         = pseudo.moic.res,
                       surv.info        = surv.info,
                       convt.time       = "y", # convert day unit to year
                       surv.median.line = "h", # draw horizontal line at median survival
                       fig.name         = "KAPLAN-MEIER CURVE OF PAM50 BY PSEUDO")


save.image(file = "alldata.Rdata")
load(file = "alldata.Rdata")


tmp.moic.res <- list("clust.res"=cmoic.brca[["clust.res"]],
                     "mo.method"=cmoic.brca[["mo.method"]])

tmb.brca <- compTMB(moic.res     = tmp.moic.res,
                    maf          = maf,
                    rmDup        = TRUE, # remove duplicated variants per sample
                    rmFLAGS      = FALSE, # keep FLAGS mutations
                    exome.size   = 38, # estimated exome size
                    test.method  = "nonparametric", # statistical testing method
                    fig.name     = "DISTRIBUTION OF TMB AND TITV TMP")

clin.brca <- compClinvar(moic.res      = tmp.moic.res,
                         var2comp      = surv.info, # data.frame needs to summarize (must has row names of samples)
                         strata        = "Subtype", # stratifying variable (e.g., Subtype in this example)
                         factorVars    = c("PAM50","pstage","fustat"), # features that are considered categorical variables
                         nonnormalVars = "futime", # feature(s) that are considered using nonparametric test
                         exactVars     = "pstage", # feature(s) that are considered using exact test
                         doWord        = TRUE, # generate .docx file in local path
                         tab.name      = "SUMMARIZATION OF CLINICAL FEATURES TMP")

tmb.brca <- compTMB(moic.res     = pseudo.moic.res,
                    maf          = maf,
                    rmDup        = TRUE, # remove duplicated variants per sample
                    rmFLAGS      = FALSE, # keep FLAGS mutations
                    exome.size   = 38, # estimated exome size
                    test.method  = "nonparametric", # statistical testing method
                    fig.name     = "DISTRIBUTION OF TMB AND TITV")
head(tmb.brca$TMB.dat)

fga.brca <- compFGA(moic.res     = pseudo.moic.res,
                    segment      = segment,
                    iscopynumber = FALSE, # this is a segmented copy number file
                    cnathreshold = 0.2, # threshold to determine CNA gain or loss
                    test.method  = "nonparametric", # statistical testing method
                    fig.name     = "BARPLOT OF FGA")

clin.brca <- compClinvar(moic.res      = pseudo.moic.res,
                         var2comp      = surv.info, # data.frame needs to summarize (must has row names of samples)
                         strata        = "Subtype", # stratifying variable (e.g., Subtype in this example)
                         factorVars    = c("PAM50","pstage","fustat"), # features that are considered categorical variables
                         nonnormalVars = "futime", # feature(s) that are considered using nonparametric test
                         exactVars     = "pstage", # feature(s) that are considered using exact test
                         doWord        = TRUE, # generate .docx file in local path
                         tab.name      = "SUMMARIZATION OF CLINICAL FEATURES")
print(clin.brca$compTab)


save(marker.up,file = "makerup_UVM.rds")
