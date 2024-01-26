



setwd("~/Desktop/TCGA_work/TCGA_UVM_project/")
###UVM 脚本 多组学  葡萄膜黑色素瘤
#80 UVM sample

library(MOVICS)

###数据准备工作####
library(MOVICS)
library(data.table)
library(tidyverse)
####mRNA和lncRNA####
#rm(list = ls())
options(stringsAsFactors = F)
if(!require(tidyverse)) install.packages('tidyverse')
if(!require(data.table)) install.packages('data.table')
dev.off()
project<-'TCGA-UVM'
##数据位置

#load("/t8a/vip39database/database/UCSC_TCGA/TCGA_methylation450/TCGA-UVM.methylation450.tsv.gz")
fpkms<-read_tsv('/t8a/vip39database/database/UCSC_TCGA/TCGA_FPKM/TCGA-UVM.htseq_fpkm.tsv.gz')%>%
  column_to_rownames(var="Ensembl_ID")%>%
  as.matrix()
fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
tpms<-as.data.frame(apply(fpkms,2,fpkmToTpm))%>%
  rownames_to_column(var='gene_id')%>%
  mutate(gene_id1=substr(.$gene_id,1,15))%>%
  dplyr::select(gene_id1,everything(),-gene_id)%>%
  write_csv(paste0(project,'tpms.csv'))
#UVM_count<-read_tsv(file = "/home/aim/Desktop/urology/Urinary_data/kidney/UVM_DATA/TCGA-UVM.htseq_counts.tsv.gz")
load("/t8a/vip39database/database/tcga_counts_fpkm_tpm/TCGA-UVM_counts_gene_symbol.Rdata")
UVM_count<-as.data.frame(fpkm2)
UVM_count[1:4,1:4]
tail(rownames(UVM_count))

###名称转化
##count
count_tmp <- UVM_count
load(file = "/home/aim/database/genome/anno.Rdata")
sum(rownames(count_tmp) %in% mRNA_anno$gene_name)
#> [1] 19712
UVM_mRNA_count<-as.data.frame(UVM_count)%>%
  rownames_to_column(var='gene_name')%>%
  inner_join(mRNA_anno,.,by=c("gene_name"))%>%#.很重要 #通过gene_id 合并注释
  dplyr::select(-gene_id,-gene_type)%>%#删除gene_id
  mutate(rowmean=apply(.[,2:length(colnames(.))],1,function(x) mean(as.numeric(x),na.rm=T)))%>%
  arrange(desc(rowmean))%>%#相同的注释取最大值
  distinct(gene_name,.keep_all = T)%>%#第1列为gene name 去除重复的
  dplyr::select(-rowmean)%>%
  column_to_rownames(var='gene_name')
UVM_mRNA_count[1:5,1:5]
#同样办法拆出lncRNA
UVM_lncRNA_count<-as.data.frame(UVM_count)%>%
  rownames_to_column(var='gene_name')%>%
  inner_join(lnc_anno,.,by=c("gene_name"))%>%#.很重要 #通过gene_id 合并注释
  dplyr::select(-gene_id,-gene_type)%>%#删除gene_id
  mutate(rowmean=apply(.[,2:length(colnames(.))],1,function(x) mean(as.numeric(x),na.rm=T)))%>%
  arrange(desc(rowmean))%>%#相同的注释取最大值
  distinct(gene_name,.keep_all = T)%>%#第1列为gene name 去除重复的
  dplyr::select(-rowmean)%>%
  column_to_rownames(var='gene_name')
UVM_lncRNA_count[1:5,1:5]

save(UVM_mRNA_count,UVM_lncRNA_count,file = paste0(project,"_count.Rdata"))
load(file = "TCGA-UVM_count.Rdata")

##FPKM

UVM_fpkm <- fpkms
UVM_fpkm[1:4,1:4]
sum(rownames(fpkms) %in% mRNA_anno$gene_id)
#> [1] 19814
UVM_mRNA_fpkm<-as.data.frame(UVM_fpkm)%>%
  rownames_to_column(var='gene_id')%>%
  inner_join(mRNA_anno,.,by=c("gene_id"))%>%#.很重要 #通过gene_id 合并注释
  dplyr::select(-gene_id,-gene_type)%>%#删除gene_id
  mutate(rowmean=apply(.[,2:length(colnames(.))],1,function(x) mean(as.numeric(x),na.rm=T)))%>%
  arrange(desc(rowmean))%>%#相同的注释取最大值
  distinct(gene_name,.keep_all = T)%>%#第1列为gene name 去除重复的
  dplyr::select(-rowmean)%>%
  column_to_rownames(var='gene_name')
UVM_mRNA_fpkm[1:5,1:5]
#同样办法拆出lncRNA
UVM_lncRNA_fpkm<-as.data.frame(UVM_fpkm)%>%
  rownames_to_column(var='gene_id')%>%
  inner_join(lnc_anno,.,by=c("gene_id"))%>%#.很重要 #通过gene_id 合并注释
  dplyr::select(-gene_id,-gene_type)%>%#删除gene_id
  mutate(rowmean=apply(.[,2:length(colnames(.))],1,function(x) mean(as.numeric(x),na.rm=T)))%>%
  arrange(desc(rowmean))%>%#相同的注释取最大值
  distinct(gene_name,.keep_all = T)%>%#第1列为gene name 去除重复的
  dplyr::select(-rowmean)%>%
  column_to_rownames(var='gene_name')
UVM_lncRNA_fpkm[1:5,1:5]

save(UVM_lncRNA_fpkm,UVM_mRNA_fpkm,file = paste0(project,"_fpkm.Rdata"))
load(file = "TCGA-UVM_fpkm.Rdata")

###tpm
load(file = "/t8a/vip39database/database/tcga_counts_fpkm_tpm/TCGA-UVM_tpm_gene_symbol.Rdata")
#> [1] 19814
UVM_mRNA_tpm<-as.data.frame(tpms)%>%
  rownames_to_column(var='gene_name')%>%
  inner_join(mRNA_anno,.,by=c("gene_name"))%>%#.很重要 #通过gene_id 合并注释
  dplyr::select(-gene_id,-gene_type)%>%#删除gene_id
  mutate(rowmean=apply(.[,2:length(colnames(.))],1,function(x) mean(as.numeric(x),na.rm=T)))%>%
  arrange(desc(rowmean))%>%#相同的注释取最大值
  distinct(gene_name,.keep_all = T)%>%#第1列为gene name 去除重复的
  dplyr::select(-rowmean)%>%
  column_to_rownames(var='gene_name')
UVM_mRNA_tpm[1:5,1:5]
#log2 处理
UVM_mRNA_tpm <- log2(UVM_mRNA_tpm+1)

#同样办法拆出lncRNA
UVM_lncRNA_tpm<-as.data.frame(tpms)%>%
  rownames_to_column(var='gene_name')%>%
  inner_join(lnc_anno,.,by=c("gene_name"))%>%#.很重要 #通过gene_id 合并注释
  dplyr::select(-gene_id,-gene_type)%>%#删除gene_id
  mutate(rowmean=apply(.[,2:length(colnames(.))],1,function(x) mean(as.numeric(x),na.rm=T)))%>%
  arrange(desc(rowmean))%>%#相同的注释取最大值
  distinct(gene_name,.keep_all = T)%>%#第1列为gene name 去除重复的
  dplyr::select(-rowmean)%>%
  column_to_rownames(var='gene_name')
UVM_lncRNA_tpm[1:5,1:5]
UVM_lncRNA_tpm <- log2(UVM_lncRNA_tpm+1)
table(str_sub(colnames(UVM_mRNA_tpm),14,15))

##只有一例normal sample
# 01 
# 80
###miRNA数据####
#load("/t8a/vip39database/database/UCSC_TCGA/TCGA_mirna/TCGA-UVM.mirna.tsv.gz")
UVM_miRNA <- read_tsv("/t8a/vip39database/database/UCSC_TCGA/TCGA_mirna/TCGA-UVM.mirna.tsv.gz")
UVM_miRNA <- as.data.frame(UVM_miRNA)
rownames(UVM_miRNA) <- UVM_miRNA$miRNA_ID
UVM_miRNA <- UVM_miRNA[,-1]
UVM_miRNA <- as.data.frame(UVM_miRNA)
colnames(UVM_miRNA) <- str_sub(colnames(UVM_miRNA),1,15)

table(colnames(mo.data$mRNA.expr)%in%colnames(UVM_miRNA))
UVM_miRNA <- UVM_miRNA[,colnames(mo.data$mRNA.expr)]
dim(UVM_miRNA)
#[1] 1881  80
##去除几乎没有表达的miRNA数据
UVM_miRNA <- UVM_miRNA[apply(UVM_miRNA,1,function(x) sum(x > 1) > 0),]
dim(UVM_miRNA)
#[1] 809 80

###甲基化数据####
library(data.table)
library(impute)
library(ChAMP)
library(stringr)
library(tibble)
options(stringsAsFactors = F)
meth <- data.table::fread('/t8a/vip39database/database/UCSC_TCGA/TCGA_methylation450/TCGA-UVM.methylation450.tsv.gz')
a <- meth
a = column_to_rownames(a,"Composite Element REF")
colnames(a)= str_sub(colnames(a),1,15)
beta=as.matrix(a)
# beta信号值矩阵里面不能有NA值
beta=impute.knn(beta) 
sum(is.na(beta))
beta=beta$data
beta=beta+0.00001
dim(beta)
beta[1:5,1:5]

save(beta,file = paste0(project,"_meth.Rdata"))
load(file = "TCGA-UVM_meth.Rdata")

###mut.status文件####
options(stringsAsFactors = F) 
require(maftools) 
require(dplyr)
#laml <- data.table::fread("/home/aim/Desktop/urology/Urinary_data/bladder/UVM_DATA/TCGA-UVM.varscan2_snv.tsv.gz")
###这个 flow过于繁琐

laml = read.maf(maf = './TCGA.UVM.mutect.1ab98b62-5863-4440-84f9-3c15d476d523.DR-10.0.somatic.maf.gz')
data1 <- fread("/t8a/vip39database/database/UCSC_TCGA/TCGA_SNV/mutect2_snv/TCGA-UVM.mutect2_snv.tsv.gz")
colnames(data1)[2] <- "Hugo_Symbol"
colnames(data1)[3] <- "Chromosome"
colnames(data1)[4:5] <- c("Start_Position","End_Position")
colnames(data1)[6] <- "Reference_Allele"
colnames(data1)[7] <- "Tumor_Seq_Allele2"
colnames(data1)[9] <- "Variant_Classification"
#clinical1 <- fread("")
laml <- read.maf(maf=data1)
laml 

###使用TCGA_biolinks 进行分析
library(TCGAbiolinks)
setwd("/home/aim/Desktop/TCGA_work/TCGA_UVM/mut_data/")
query <- GDCquery(
  project = "TCGA-UVM", 
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  access = "open"
)

GDCdownload(query)
GDCprepare(query, save = T,save.filename = "TCGA-UVM_SNP.Rdata")
##总是有duplicated sample
# 
# library(maftools)
# load(file = "/home/aim/Desktop/TCGA_work/TCGA_UVM/mut_data/")
# maf.coad <- data

###使用技能树脚本分析
require(maftools)
options(stringsAsFactors = F) 
library(data.table)
#tmp=fread('TCGA-BRCA.mutect2_snv.tsv.gz')
tmp <- fread("/t8a/vip39database/database/UCSC_TCGA/TCGA_SNV/mutect2_snv/TCGA-UVM.mutect2_snv.tsv.gz")
head(tmp)   
colnames(tmp) =c( "Tumor_Sample_Barcode", "Hugo_Symbol", 
                  "Chromosome", "Start_Position", 
                  "End_Position", "Reference_Allele", "Tumor_Seq_Allele2", 
                  "HGVSp_Short" , 'effect' ,"Consequence",
                  "vaf" ) 
tmp$Entrez_Gene_Id =1
tmp$Center ='ucsc'
tmp$NCBI_Build ='GRCh38'
tmp$NCBI_Build ='GRCh38'
tmp$Strand ='+'
tmp$Variant_Classification = tmp$effect
tail(sort(table(tmp$Variant_Classification )))
tmp$Tumor_Seq_Allele1 = tmp$Reference_Allele
tmp$Variant_Type = ifelse(
  tmp$Reference_Allele %in% c('A','C','T','G') & tmp$Tumor_Seq_Allele2 %in% c('A','C','T','G'),
  'SNP','INDEL'
)
table(tmp$Variant_Type )
tcga.UVM = read.maf(maf = tmp,
                     vc_nonSyn=names(tail(sort(table(tmp$Variant_Classification )))))

oncoplot(maf = tcga.UVM, top = 10) # 高频突变的前10个基因

laml <- tcga.UVM
maf_df = laml@data
mut.status <- maftools::mutCountMatrix(maf = laml)
mut.status <- ifelse(mut.status==0,0,1)
range(mut.status)
colnames(mut.status)= str_sub(colnames(mut.status),1,15)
table(str_sub(colnames(mut.status),14,15))
# 01  06 
# 104 363
table(colnames(mut.status)%in%colnames(beta))
colnames(beta)[1:5]
colnames(mut.status)[1:5]
#选出前30个mut
mut.statu <- mut.status[c(1:30),]
save(mut.statu,file = paste0(project,"_mut.Rdata"))
load(file = "TCGA-UVM_mut.Rdata")


####获取共有数据，构建mo.data#####
table(colnames(UVM_mRNA_count)%in%colnames(UVM_lncRNA_fpkm))
colnames(UVM_lncRNA_count)[1:5]
colnames(UVM_lncRNA_count)= str_sub(colnames(UVM_lncRNA_count),1,15)
colnames(UVM_mRNA_count)= str_sub(colnames(UVM_mRNA_count),1,15)
colnames(UVM_lncRNA_fpkm)= str_sub(colnames(UVM_lncRNA_fpkm),1,15)
colnames(UVM_mRNA_fpkm)= str_sub(colnames(UVM_mRNA_fpkm),1,15)

allid <- intersect(colnames(UVM_mRNA_count),intersect(colnames(beta),colnames(mut.statu)))
length(allid)
#[1] 80
table(str_sub(allid,14,15))
##全部都是肿瘤sample
dim(UVM_mRNA_fpkm)

UVM_mRNA_fpkm <- round(UVM_mRNA_fpkm,2)
UVM_lncRNA_fpkm <- round(UVM_lncRNA_fpkm,2)
meth.beta <- round(beta,3)
#mut.status <- mut.statu
##并集
UVM_mRNA_fpkm <- UVM_mRNA_fpkm[,allid]
dim(UVM_mRNA_fpkm)
UVM_lncRNA_fpkm <- UVM_lncRNA_fpkm[,allid]
dim(UVM_lncRNA_fpkm)
meth.beta <- meth.beta[,allid]
dim(meth.beta)
mut.status <- mut.status[,allid]
mut.status <- mut.status[c(1:30),]
dim(mut.status)

###选出对于预后有影响的signature

##精简数据
tmp<- UVM_mRNA_fpkm # get expression data with 500 features
elite.tmp1 <- getElites(dat       = tmp,
                        method    = "mad",
                        elite.pct = (500/nrow(tmp))) # this time only top 10% features with high mad values are kept
dim(elite.tmp1$elite.dat) # get 50 elite left只留下0.1的数据

tmp<- UVM_lncRNA_fpkm # get expression data with 500 features
elite.tmp2 <- getElites(dat       = tmp,
                        method    = "mad",
                        elite.pct = (500/nrow(tmp))) # this time only top 10% features with high mad values are kept
dim(elite.tmp2$elite.dat) # get 50 elite left只留下0.1的数据

tmp<- meth.beta # get expression data with 500 features
elite.tmp3 <- getElites(dat       = tmp,
                        method    = "mad",
                        elite.pct = (1000/nrow(tmp))) # this time only top 10% features with high mad values are kept
dim(elite.tmp3$elite.dat) # get 50 elite left只留下0.1的数据

range(mut.status)



range(UVM_mRNA_fpkm)
###使用预后数据进行区分
##精简数据
tmp<- UVM_mRNA_fpkm[,allid]# get expression data with 500 features
surv.info <- clin.info
elite.tmp1 <- getElites(dat       = tmp,
                        method    = "cox",
                        surv.info = surv.info, # must provide survival information with 'futime' and 'fustat'
                        p.cutoff  = 0.01)
#elite.num = 100) # this time elite.num argument will be disabled because cox method refers to p.cutoff to select elites
#> --all sample matched between omics matrix and survival data.
#> 5% 10% 15% 20% 25% 30% 35% 40% 45% 50% 55% 60% 65% 70% 75% 80% 85% 90% 95% 100%
dim(elite.tmp1$elite.dat) # get 125 elites
#> [1] 3179  465
table(elite.tmp1$unicox$pvalue < 0.01) # 125 genes have nominal pvalue < 0.05 in univariate Cox regression


elite.tmp1 <- getElites(dat       = elite.tmp1$elite.dat,
                        method    = "mad",
                        elite.num = 1000) # this time only top 10% features with high mad values are kept
dim(elite.tmp1$elite.dat) # get 50 elite left只留下0.1的数据
#[1] 1000  465


tmp<- UVM_lncRNA_fpkm[,allid] # get expression data with 500 features

elite.tmp2 <- getElites(dat       = tmp,
                        method    = "cox",
                        surv.info = surv.info, # must provide survival information with 'futime' and 'fustat'
                        p.cutoff  = 0.05)
#elite.num = 100) # this time elite.num argument will be disabled because cox method refers to p.cutoff to select elites
#> --all sample matched between omics matrix and survival data.
#> 5% 10% 15% 20% 25% 30% 35% 40% 45% 50% 55% 60% 65% 70% 75% 80% 85% 90% 95% 100%
dim(elite.tmp2$elite.dat) # get 125 elites
#> [1]4413   80
table(elite.tmp2$unicox$pvalue < 0.05) # 125 genes have nominal pvalue < 0.05 in univariate Cox regression

elite.tmp2 <- getElites(dat       = elite.tmp2$elite.dat,
                        method    = "mad",
                        elite.num = 1000) # this time only top 10% features with high mad values are kept
dim(elite.tmp2$elite.dat) # get 50 elite left只留下0.1的数据
#tmp <- elite.tmp2$elite.dat


tmp<- UVM_miRNA[,allid] # get expression data with 500 features
elite.tmp6 <- getElites(dat       = tmp,
                        method    = "cox",
                        surv.info = surv.info, # must provide survival information with 'futime' and 'fustat'
                        p.cutoff  = 0.05)
#elite.num = 100) # this time elite.num argument will be disabled because cox method refers to p.cutoff to select elites
#> --all sample matched between omics matrix and survival data.
#> 5% 10% 15% 20% 25% 30% 35% 40% 45% 50% 55% 60% 65% 70% 75% 80% 85% 90% 95% 100%
dim(elite.tmp6$elite.dat) # get 125 elites
#> [1] 212  80
table(elite.tmp6$unicox$pvalue < 0.05) # 125 genes have nominal pvalue < 0.05 in univariate Cox regression
#tmp <- elite.tmp6$elite.dat

###这一步计算耗时较久 需要保存结果
tmp<- meth.beta[,allid] # get expression data with 500 features
elite.tmp4 <- getElites(dat       = tmp,
                        method    = "cox",
                        surv.info = surv.info, # must provide survival information with 'futime' and 'fustat'
                        p.cutoff  = 0.01)

#elite.num = 100) # this time elite.num argument will be disabled because cox method refers to p.cutoff to select elites
#> --all sample matched between omics matrix and survival data.
#> 5% 10% 15% 20% 25% 30% 35% 40% 45% 50% 55% 60% 65% 70% 75% 80% 85% 90% 95% 100%
dim(elite.tmp4$elite.dat) # get 125 elites
#> [1] 60433    80
table(elite.tmp4$unicox$pvalue < 0.01) # 125 genes have nominal pvalue < 0.05 in univariate Cox regression

save(elite.tmp4,file = "~/Desktop/TCGA_work/TCGA_UVM/UVM_data/meth_cox.rds")
load(file = "./meth_cox.rds")
##甲基化位点太多了
elite.tmp4 <- getElites(dat       = elite.tmp4$elite.dat,
                        method    = "mad",
                        elite.num = 1000) # this time only top 10% features with high mad values are kept
dim(elite.tmp4$elite.dat) # get 50 elite left只留下0.1的数据
tmp <- elite.tmp4$elite.dat


mo.data <- list(mRNA.expr = elite.tmp1$elite.dat,
                lnc.expr=elite.tmp2$elite.dat,
                meth.beta=elite.tmp3$elite.dat,
                mut.status=mut.status)
if(F){mo.data <- list(mRNA.expr = mRNA.expr,
                      lnc.expr=lnc.expr,
                      meth.beta=meth.beta,
                      mut.status=mut.status)}
save(mo.data,file = paste0(project,"_modata.Rdata"))


###maf文件数据
maf_df$Tumor_Sample_Barcode
maf_df$Hugo_Symbol
maf_df$Chromosome
maf_df$Start_Position
maf_df$End_Position
maf_df$Variant_Classification
maf_df$Reference_Allele
maf_df$Tumor_Seq_Allele1
maf_df$Tumor_Seq_Allele2
maf <- maf_df[,c("Tumor_Sample_Barcode","Hugo_Symbol","Chromosome",
                 "Start_Position","End_Position","Variant_Classification","Variant_Type",
                 "Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2")]
maf$Tumor_Sample_Barcode <- str_sub(maf$Tumor_Sample_Barcode,1,15)
maf <- maf%>%
  filter(Tumor_Sample_Barcode%in%allid)
table(allid%in%maf$Tumor_Sample_Barcode)
length(table(maf$Tumor_Sample_Barcode))
save(maf,file = paste0(project,"_maf.Rdata"))
load(file = "TCGA-UVM_maf.Rdata")

###segment
CNS<-read_tsv('/t8a/vip39database/database/UCSC_TCGA/TCGA_CNV/TCGA-UVM.cnv.tsv.gz')

segment <- CNS
segment$sample <- str_sub(segment$sample,1,15)
table(allid%in%segment$sample)
#save(segment,file = paste0(project,"_segment.Rdata"))
#load(file = "TCGA-UVM_segment.Rdata")


###clin.info，利用allid取出
clin<-read_tsv('/home/aim/Desktop/TCGA_work/TCGA_UVM_project/UVM_data/TCGA-UVM.survival.tsv')
pd <- read_tsv('/home/aim/Desktop/TCGA_work/TCGA_UVM_project/UVM_data/TCGA-UVM.GDC_phenotype.tsv.gz')
allclin<-read_tsv('/home/aim/Desktop/urology/Urinary_data/kidney/KIRC_DATA/Survival_SupplementalTable_S1_20171025_xena_sp.gz')

clin$sample <- str_sub(clin$sample,1,15)##只要id的前15位
table(allid%in%clin$sample)
rownames(allclin) <- allclin$sample

pd$submitter_id.samples <- str_sub(pd$submitter_id.samples,1,15)
table(allid%in%pd$submitter_id.samples)
pd <- pd[!duplicated(pd$submitter_id.samples),]
rownames(pd) <- pd$submitter_id.samples
table(allid%in%pd$submitter_id.samples)

clinneed <- allclin[allid,]
pdneed <- pd[allid,]
rownames(clinneed) <- clinneed$sample
rownames(pdneed) <- pdneed$submitter_id.samples
table(rownames(clinneed)==rownames(pdneed))
clin.info <- data.frame(fustat=clinneed$OS,
                        futime=clinneed$OS.time,
                        AJCC=clinneed$ajcc_pathologic_tumor_stage,
                        pstage=pdneed$pathologic_T,
                        age=clinneed$age_at_initial_pathologic_diagnosis
)
rownames(clin.info) <- clinneed$sample
clin.info$pstage <- str_sub(clin.info$pstage,1,2)
clin.info$pstage <- as.character(clin.info$pstage)
table(is.na(clin.info$pstage))
table(clin.info$pstage)

table(clin.info$AJCC)
table(is.na(clin.info$AJCC))
clin.info[rownames(subset(clin.info,AJCC=="[Discrepancy]")),"AJCC"] <- "Stage X"
clin.info[is.na(clin.info$AJCC),"AJCC"] <- "Stage X"
table(clin.info$AJCC)
clin.info[,c("A","B")] <- str_split_fixed(clin.info$AJCC," ",2)
clin.info <- clin.info[,-c(3,6)]
colnames(clin.info)[5] <- "AJCC"
clin.info <- clin.info[,c(1,2,5,3,4)]

save(clin.info,file = paste0(project,"_clin_info.Rdata"))
load(file="TCGA-UVM_clin_info.Rdata")

mRNA_count <- UVM_mRNA_count[,allid]
dim(mRNA_count)
mRNA_fpkm <- UVM_mRNA_fpkm[,allid]
dim(mRNA_fpkm)
muttmp <- as.data.frame(mut.status)
mRNA_count <- as.data.frame(mRNA_count)
mRNA_fpkm <- as.data.frame(mRNA_fpkm)


UVM.tcga <- list(mRNA.expr = elite.tmp1$elite.dat,
                  lnc.expr=elite.tmp2$elite.dat,
                  mi.expr=elite.tmp6$elite.dat,
                  meth.beta=elite.tmp4$elite.dat,
                  mut.status=muttmp,
                  count=mRNA_count,
                  fpkm=mRNA_fpkm,
                  maf=maf,
                  segment=segment,
                  clin.info=clin.info)
save(UVM.tcga,file = "UVM_tcga.Rdata")
