library(ggplot2)
library(gplots)
library(argparse)
source("/humgen/atgu1/fs03/berylc/MuscDisease/QC/MuscleCheck/scripts/bin.R")

parser <- ArgumentParser(description='Compare patient RPKMs with GTEx tissues to check validity of tissue identity')
parser$add_argument('-tissues', action='store',help='Tissues (Current options: liver,muscle,skin_sun,skin_nosun,subc_adipose,fibroblast. Default is all except liver',default="muscle,skin_sun,skin_nosun,subc_adipose,fibroblast")
parser$add_argument('-gtex_rpkm',action='store',help='Output from MakeRPKMfile.R, default in place',default='/humgen/atgu1/fs03/berylc/MuscDisease/QC/MuscleCheck/data/GTEx_5Tissues_MeleGenes_muscle,skin,adipose,fibroblasts.txt')
parser$add_argument('-patient_rpkm',action='store',help='Patient RPKM file e.g. *.genes.rpkm.gct file from RNASeQC output')
parser$add_argument('-outfile',action='store',help='File name to create and plot PCAs to (no.pdf in the end)')
parser$add_argument('-writePCADat',action='store_true',help='Write PCA file, Boolean option',default=FALSE)

args <- parser$parse_args()


tissue_list = strsplit(args$tissues,split=",")[[1]]

nClus = length(tissue_list)
#Subctutaneous adipose&fibroblasts and sun exposed and non sun exposed skin cluster together, so if they're in the tissue list, count them each as one cluster
if(length(grep("skin",tissue_list))==2 | length(grep(paste("subc_adipose","fibroblast",sep='|'),tissue_list))==2){nClus = length(tissue_list)-1}
#11/23/2015 actually subcuataneous adipse and fibroblast shouldn't be clustering together, but I am keeping the command because in PCA space they do seem to be clustering

if(length(grep("skin",tissue_list))==2 && length(grep(paste("subc_adipose","fibroblast",sep='|'),tissue_list))==2){nClus = length(tissue_list)-2}



#gtexRPKMfile = "/humgen/atgu1/fs03/berylc/MuscDisease/QC/MuscleCheck/data/GTEx_6Tissues_MeleGenes.txt"
#patientRPKMfile="/humgen/atgu1/fs03/taru/mendelian/RNA-SeQC/genes.rpkm.gct"

all_gtex<-read.delim(args$gtex_rpkm,header=T,row.names=1,stringsAsFactors=F)
all_gtex<-all_gtex[grep(paste(tissue_list,collapse="|"),rownames(all_gtex)),]
all_patient<-read.delim(args$patient_rpkm,header=T,skip=2,stringsAsFactors=F)

colnames(all_patient)<-gsub("\\.","-",colnames(all_patient))
#This is kind of messed up since it's different in bin.R so you have to fix it each time
subsetMeleGenes(all_patient,"patient",makecolnames=F)

#Added this 01/15/16 :in because mele_useful now contains genes for colon etc. So need to subset down to the same number of columns
patient_mele<-patient_mele[,colnames(all_gtex)]


all_tissues<-rbind(all_gtex,patient_mele)

#Normalize values
all_tissues<-all_tissues+1
all_tissues<-log2(all_tissues)


pdfFile = paste(args$outfile,".pdf",sep="")
pdf(pdfFile)

thrTissPCA<-prcomp(as.matrix(all_tissues), retx=TRUE,na.action=na.omit,center=TRUE) 
plot(thrTissPCA,type="l",main="Variance Explained by PCs")


PCADat<-data.frame(thrTissPCA$x)

cols = ncol(PCADat)
#Indicate tissues for the coloring
PCADat[,cols+1]<-sapply(strsplit(rownames(PCADat),split="\\."), function(x) x[1])
PCADat[-grep(paste(tissue_list,collapse="|"),PCADat[,cols+1]),cols+1]<-"Patient muscle"

names(PCADat)[cols+1]<-"tissue"

#This is just for a pretty graph
PCADat[PCADat$tissue=="muscle",cols+1]<-"Muscle"
PCADat[PCADat$tissue=="skin_sun",cols+1]<-"Skin - lower leg"
PCADat[PCADat$tissue=="skin_nosun",cols+1]<-"Skin-suprapubic"
PCADat[PCADat$tissue=="subc_adipose",cols+1]<-"Adipose Subcutaneous"
PCADat[PCADat$tissue=="fibroblast",cols+1]<-"Transformed fibroblasts"

ggplot(PCADat,aes(x=PC1,y=PC2,col=tissue))+geom_point(size=3)+theme_bw()+theme(plot.title=element_text(size=20),axis.text=element_text(size=15),axis.title=element_text(size=14),legend.title=element_blank(),legend.text=element_text(size=10))+ggtitle("PC1 vs PC2")
 
ggplot(PCADat,aes(x=PC1,y=PC3,col=tissue))+geom_point(size=3)+theme_bw()+theme(plot.title=element_text(size=20),axis.text=element_text(size=15),axis.title=element_text(size=14),legend.title=element_blank(),legend.text=element_text(size=10))+ggtitle("PC1 vs PC3")
# 
ggplot(PCADat,aes(x=PC2,y=PC3,col=tissue))+geom_point(size=3)+theme_bw()+theme(plot.title=element_text(size=20),axis.text=element_text(size=15),axis.title=element_text(size=14),legend.title=element_blank(),legend.text=element_text(size=10))+ggtitle("PC2 vs PC3")

dev.off()

#Hiearchichal clustering
hclustfunc <- function(x) hclust(x, method="average")
distfunc <- function(x) as.dist(1-cor(t(x),method = "spearman"))
d <- distfunc(all_tissues)
fit <- hclustfunc(d)
 

pdf("/humgen/atgu1/fs03/berylc/MuscDisease/QC/MuscleCheck/scripts/plotdendogram.pdf",width=800,height=100)

par(cex=3,font=10)
plot(fit,hang=-1, main="Dendrogram of Ward's Method")

dev.off() 
 
allClusters<-cutree(fit,k=nClus)
for(clusnum in 1:nClus){
  clus<- allClusters[allClusters==clusnum]
  clusnames<-sapply(strsplit(names(clus),split="\\."), function(x) x[1])
  clusnames = unique(clusnames)
  if("muscle" %in% clusnames){
    cat(paste("\n","Patient Samples that cluster with GTEx muscle samples"))
    for(elem in sort(clusnames)){
      if(elem=="muscle"){next}
      if(substr(elem, 1, 1)=="X"){cat(paste("\n",substr(elem,2,nchar(elem))))} 
      else{cat(paste("\n",elem))}}}
   
  if(!"muscle" %in% clusnames && length(setdiff(clusnames,tissue_list))!=0){
    cat(paste("\n\n","Patient samples that do not cluster tightly with muscle:"))
     
    for(eachSamp in setdiff(clusnames,tissue_list)){cat(paste("\n",eachSamp))                                                   
    cat(paste("\n","PCA coordinates are:",paste(PCADat[eachSamp,c("PC1","PC2")],collapse=","),'\n'),"\n") 
                                                     
     }}}
   
if(args$writePCADat){
  PCADatFiletoWrite = PCADat[PCADat[,cols+1]=="Patient muscle",]
  PCADatFiletoWrite = PCADatFiletoWrite[order(rownames(PCADatFiletoWrite)),1:5]
  Name = paste(args$outfile,".PCA.txt",sep="")
  write.table(PCADatFiletoWrite, file = Name, quote=F, row.names = T, col.names=T, sep="\t") }

