

###UVM突变数据分析
library(MOVICS)
library(data.table)
library(maftools)

laml = read.maf(maf = '/t8a/vip39database/database/GDC_TCGA/TCGA_SNV/mutec2/TCGA.UVM.mutect.6c7b01bc-b068-4e01-8b4d-0362f5959f65.DR-10.0.somatic.maf.gz')
laml 
oncoplot(maf = laml, top = 15) # 高频突变的前10个基因
sub <- moic.res.list$iClusterBayes$clust.res
sub$Sample <- sub$samID
sub$Cluster <- paste0("C",sub$clust)
head(sub)
sub$Cluster
sub$Sample
idC1 <- sub$Sample[sub$Cluster=="C1"]
idC2 <- sub$Sample[sub$Cluster=="C2"]


library(stringr)
vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)



#mafC1 <- subsetMaf(maf = laml, tsb = str_sub(idC1,1,12),  isTCGA = TRUE)
#mafC2 <- subsetMaf(maf = laml, tsb = str_sub(idC2,1,12), isTCGA = TRUE)
tmp=fread('/t8a/database/UCSC_TCGA/TCGA_SNV/mutect2_snv/TCGA-UVM.mutect2_snv.tsv.gz')
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

###读取maf 单个文件
# 设置工作地址为解压后的文件地址
setwd("/home/aim/Desktop/TCGA_work/TCGA_UVM_project/mut_data/")
# 创建一个目录保存所有结果
dir.create('0000_all_maf')
# 下载的文件名为36个字母才是所需的，其他的可以忽略
dir_all <- dir()[nchar(dir()) == 36]
for (dir_maf in dir_all) {
  #内部文件也是压缩的，需要解压出来并保存到之前创建的目录中
  maf_file <- list.files(dir_maf, pattern = ".*maf")
  if (grepl('gz$',maf_file)){
    R.utils::gunzip(paste0(dir_maf,"/",maf_file))
  }
  file_extracted <- list.files(dir_maf, pattern = ".*maf$")
  file.copy(paste0(dir_maf,"/",file_extracted),"0000_all_maf")
}
# 将工作地址设置为之前创建的目录
setwd("/home/aim/Desktop/TCGA_work/TCGA_UVM_project/mut_data/0000_all_maf")

file_extracted_maf <- list.files()
first_file <- read.delim(file_extracted_maf[1], header = T, sep = '\t', comment.char = '#',stringsAsFactors = F)
for (extracted_maf in file_extracted_maf[2:length(file_extracted_maf)]) {
  file_appended <- read.delim(extracted_maf, header = T, sep = '\t', comment.char = '#',stringsAsFactors = F)
  first_file <- rbind(first_file,file_appended)
}
UVM_maf <- first_file

UVM_maf[1:4,1:4]
#save(UVM_maf,file = "~/Desktop/TCGA_work/TCGA_UVM_project/UVM_maf.rds")
load("~/Desktop/TCGA_work/TCGA_UVM_project/UVM_maf.rds")
tmp <- UVM_maf
library(stringr)
library(maftools)
UVM_maf_C1 <- tmp[str_sub(tmp$Tumor_Sample_Barcode,1,15)%in%idC1,]
mafC1 <- read.maf(UVM_maf_C1)

UVM_maf_C2 <- UVM_maf[str_sub(UVM_maf$Tumor_Sample_Barcode,1,15)%in%idC2,]
mafC2 <- read.maf(UVM_maf_C2)


library(stringr)
vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)

mafall <- read.maf(UVM_maf)
par(oma=c(1,1,1,1), mar=c(2,2,2,2))

oncoplot(maf = mafall,colors = vc_cols,
         #genes = unique(RNA_gene),
         ##colors = paletteer_c("scico::berlin", n = 3),#给突变配色
         #annotationColor = annocolors, #给临床信息配色
         #top = 20,
         sortByAnnotation = TRUE,
         # clinicalFeatures = c("Subtype"),
         writeMatrix =T)

oncoplot(maf = mafC1,colors = vc_cols,
         ##colors = paletteer_c("scico::berlin", n = 3),#给突变配色
         #annotationColor = annocolors, #给临床信息配色
         #top = 20,
         sortByAnnotation = TRUE,
         # clinicalFeatures = c("Subtype"),
         writeMatrix =T)

oncoplot(maf = mafC2,colors = vc_cols,
         ##colors = paletteer_c("scico::berlin", n = 3),#给突变配色
         #annotationColor = annocolors, #给临床信息配色
         #top = 20,
         sortByAnnotation = TRUE,
         # clinicalFeatures = c("Subtype"),
         writeMatrix =T)



###两组比较

C1.vs.C2 <- mafCompare(m1 = mafC1, m2 = mafC2, m1Name = 'C1', m2Name = 'C2', minMut = 5)
print(C1.vs.C2)

forestPlot(mafCompareRes = C1.vs.C2, pVal = 0.05,geneFontSize = 0.5)
coBarplot(m1 = mafC1, m2 = mafC2,colors = vc_cols,m1Name = "C1", m2Name = "C2")

C1.titv = titv(maf = mafC1, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = C1.titv)
C2.titv = titv(maf = mafC2, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = C2.titv)


dgi = drugInteractions(maf = mafC1, fontSize = 0.75)

dgi = drugInteractions(maf = mafC2, fontSize = 0.75)



OncogenicPathways(maf = mafC1)
OncogenicPathways(maf = mafC2)
OncogenicPathways(maf = mafC3)

if(!require(paletteer))install.packages("paletteer")
if(!require(scico))install.packages('scico')
if(!require(nord))install.packages('nord')
library(paletteer)
somaticInteractions(maf = mafC1, top = 25,colPal = "GnBu",showSum = F,sigSymbolsSize = 3,pvalue = c(0.01, 0.05))
library(eoffice)

somaticInteractions(maf = mafC2, top = 25,colPal = "GnBu",showSum = F,sigSymbolsSize = 3,pvalue = c(0.01, 0.05))

par(oma=c(1,1,1,1), mar=c(2,2,2,2))
somaticInteractions(maf = mafC2, 
                    #genes = as.character(tmp[tmp$pval < 0.05,"Hugo_Symbol"]),# 取出上面分析得到的显著基因
                    showCounts = FALSE, # 不展示OR值
                    showSum = FALSE, # 不展示突变综述
                    
                    pvalue = c(0.05, 0.01), # p值的区间（原文这里也有问题，P<0.05显著的是.，但是这个函数好像没有修改的参数）
                    colPal = "PiYG") # 选择原文配色
par(oma=c(1,1,1,1), mar=c(2,2,2,2))
somaticInteractions(maf = mafC1, 
                    #genes = as.character(tmp[tmp$pval < 0.05,"Hugo_Symbol"]),# 取出上面分析得到的显著基因
                    showCounts = FALSE, # 不展示OR值
                    showSum = FALSE, # 不展示突变综述 
                    
                    pvalue = c(0.05, 0.01), # p值的区间（原文这里也有问题，P<0.05显著的是.，但是这个函数好像没有修改的参数）
                    colPal = "PiYG") # 选择原文配色
