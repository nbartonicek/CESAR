# Transcription factor ELF5 is essential for production of milk in mammals,
# but has been observed to contribute to resistance to endocrine therapies 
# targeting estrogen receptor (ER) in breast cancer.

# In order to delineate how ELF5 modulates ER binding we performed transcriptomic
# analysis (RNAseq) to observe expression changes which occur under induced 
# expression of ELF5, and integrated it with analysis of binding of transcription 
# factors ELF5, ER and FOXA1. We showed that elevated levels of ELF5 change locations
# of regular ER binding towards the areas that cause acquisition of resistance to 
# drug therapies. The results of the analysis have been published in Piggin et al
# (2020), where I am the second-last author and produced main figures 1-4,6 and most of 
# the supplementary figures. https://pubmed.ncbi.nlm.nih.gov/31895944/

# This script below performs the first part of the analysis and visualization involved in the study - 
# classic transcriptomic analysis.

# The key steps include
# 1. Data import
# 2. Data wrangling and integration
# 3. Statistical analysis
# 4. Visualization 
# 5. Answering specific biological questions

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("librarian")
BiocManager::install("org.Hs.eg.db")
librarian::shelf(ggplot2,edgeR,RColorBrewer,pheatmap,org.Hs.eg.db,ggrepel,kegga)


#########################################
########## 1. Directory structure
#########################################

projectDir<-"."
imageDir<-paste0(projectDir,"/figures/")
robjectsDir<-paste0(projectDir,"/Robjects/")
tableDir<-paste0(projectDir,"/tables/")
rnaseqDir<-paste0(projectDir,"/RNAseq/")
annotationDir<-paste0(projectDir,"/annotation/")
system(paste("mkdir -p",imageDir))
system(paste("mkdir -p",robjectsDir))
system(paste("mkdir -p",tableDir))

#########################################
########## 2. Load in the expression data 
#########################################

inFiles<-list.files(rnaseqDir,pattern="genes",full.names=T)

#collate data from multiple tables as outputs of NGS pipelines
results<-list()
for(inFile in inFiles){
data<-read.table(inFile,header=T)
sampleName<-basename(inFile) 
sampleName<-gsub(".*V5_","",sampleName)
sampleName<-gsub("_.*","",sampleName)
cat(sampleName)
cat("\n")
results[[sampleName]]<-data$expected_count
}
df<-do.call("cbind",results)
df<-as.data.frame(df)
row.names(df)<-data$gene_id

#MDS plot for QC of technical replicates
group <- c(rep("minus",3),rep("plus",3))
y <- DGEList(df,group=group)
y <- calcNormFactors(y) 
cols <- brewer.pal(8,"Set1")[factor(group)]
pdf(paste0(imageDir,"samples_QC.pdf"),width=12,height=8)
mds <- plotMDS(y, top=200, col=cols,labels=colnames(df))
dev.off()
plotMDS(y, top=200, col=cols,labels=colnames(df))

#annotate data with external databases to get gene names
egENSEMBL <- toTable(org.Hs.egENSEMBL)
ensemblIDs=gsub("\\..*","",row.names(df))
m <- match(ensemblIDs, egENSEMBL$ensembl_id)
df$EntrezGene<-egENSEMBL$gene_id[m]

egSYMBOL <- toTable(org.Hs.egSYMBOL)
m <- match(df$EntrezGene, egSYMBOL$gene_id)
df$symbol<-egSYMBOL$symbol[m]

#eliminate duplicated symbols
o <- order(rowSums(df[,1:6]), decreasing=TRUE)
df=df[o,]
d<-duplicated(df$symbol)
df<-df[!d,]

#filter out genes that do not have at least 10 reads in at least 2 samples
include<-apply(df[,1:6],1,function(x){sum(x>=10)>=2})
df<-df[include,]
#eliminate rows without symbol
df=df[!is.na(df$symbol),]
row.names(df)<-df$symbol

#########################################
########## 3. Create differential gene expression data 
#########################################

#build generalised linear model of the data with package edgeR
treatment <- c(rep("minus",3),rep("plus",3))
design <- model.matrix(~0+treatment)
expr <- DGEList(counts=df[,1:6])
expr <- calcNormFactors(expr)
expr <- estimateDisp(expr,design)
fit <- glmQLFit(expr,design)
lrt <- glmQLFTest(fit,contrast=c(-1,1))
ELF5_expression<-as.data.frame(topTags(lrt,n=20000))

cat(paste0("Number of differentially expressed genes (FDR<0.01) is:"),
	"\nUPREGULATED: ",length(ELF5_expression[ELF5_expression$logFC>log(2) & ELF5_expression$FDR<0.01,1]),
	"\nDOWNREGULATED: ",length(ELF5_expression[ELF5_expression$logFC<(-log(2)) & ELF5_expression$FDR<0.01,1]))


#########################################
########## 4. Visualise data and critical genes 
#########################################

#take in an external list of Estrogen associated genes and plot them on a volcano plot
ER_pathway_genes<-read.table(paste0(annotationDir,"ER_associated_genes.txt"),header=T)[,1]
ER_pathway_genes_DE<-ER_pathway_genes[ER_pathway_genes %in% row.names(ELF5_expression)[1:1000]]
ELF5_expression$symbol<-row.names(ELF5_expression)
ELF5_expression$genelabels <- ""
ELF5_expression$genelabels[ELF5_expression$symbol %in% ER_pathway_genes_DE] <- ELF5_expression$symbol[ELF5_expression$symbol %in% ER_pathway_genes_DE]
ELF5_expression$color<-"gray"
ELF5_expression$color[ELF5_expression$symbol %in% ER_pathway_genes_DE]<-"darkred"

#create PDF of a volcano plot: p-value on y axis and log fold change of differential expression on X.
pdf(paste0(imageDir,"volcano_with_ER_pathway.pdf"),width=16,height=16)
  p<-ggplot(ELF5_expression) +
    geom_point(aes(x=logFC, y=-log10(FDR)),colour=ELF5_expression$color) +
    ggtitle("ER associated") +
    xlab("log2 fold change") + xlim(c(-6,6)) +
    ylab("-log10 adjusted p-value") +
    geom_hline(yintercept=1, linetype="dashed", color = "red", size=1) +
    geom_text_repel(data=subset(ELF5_expression,logFC<0),
	  max.overlaps = Inf,
      nudge_x = -6 - subset(ELF5_expression,logFC<0)$logFC, 
      segment.size  = 0.6,
      segment.colour = "gray40",
      aes(x = logFC, y = -log10(FDR), 
        label = genelabels)) +
    geom_text_repel(data=subset(ELF5_expression,logFC>0),
      max.overlaps = Inf,
      nudge_x = 6 - subset(ELF5_expression,logFC>0)$logFC,       
      segment.size  = 0.6,
      segment.colour = "gray40",
      aes(x = logFC, y = -log10(FDR), 
        label = genelabels)) +
    theme(legend.position = "none",
          plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = rel(1.25))) 
  print(p)
dev.off()
#plot the figure
print(p)

minus<-ELF5_expression[ELF5_expression$logFC<=-1,]
plus<-ELF5_expression[ELF5_expression$logFC>=1,]
significantIDs<-c(row.names(ELF5_expression[1:50,]))

fittedValues<-fitted.values(lrt)
fittedValues50<-fittedValues[row.names(fittedValues) %in% significantIDs,]
colsums<-apply(fittedValues50[,4:6],1,median)
pdf(paste0(imageDir,"heatmap_top50.pdf"),width=6,height=10)
pheatmap(log1p(fittedValues50),cex=0.7,cluster_rows = F,cluster_cols=F)
dev.off()


#########################################
########## 5. Pathway analysis 
#########################################


df<-do.call("cbind",results)
df<-as.data.frame(df)
row.names(df)<-data$gene_id

#annotate data with external databases to get gene names
egENSEMBL <- toTable(org.Hs.egENSEMBL)
ensemblIDs=gsub("\\..*","",row.names(df))
m <- match(ensemblIDs, egENSEMBL$ensembl_id)
df$EntrezGene<-egENSEMBL$gene_id[m]
o <- order(rowSums(df[,1:(dim(df)[2]-2)]), decreasing=TRUE)
df=df[o,]
d<-duplicated(df$EntrezGene)
df<-df[!d,]
df<-df[!is.na(df$EntrezGene),]
row.names(df)<-df$EntrezGene
dfClean<-as.matrix(df[,1:(dim(df)[2]-2)])

treatment <- c(rep("minus",3),rep("plus",3))
design <- model.matrix(~0+treatment)
expr <- DGEList(counts=df[,1:6])
expr <- calcNormFactors(expr)
expr <- estimateDisp(expr,design)
fit <- glmQLFit(expr,design)
fitEntrez <- glmQLFTest(fit,contrast=c(-1,1))

keg <- kegga(fitEntrez, species="Hs")
topKEGG(keg)
write.table(topKEGG(keg,number=Inf),paste0(tableDir,"KEGG_categories.txt"),row.names=T,quote=F,sep="\t")




