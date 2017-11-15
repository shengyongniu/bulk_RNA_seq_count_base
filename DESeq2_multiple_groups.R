source("http://bioconductor.org/biocLite.R")
#biocLite('DESeq2')
library('DESeq2')
library("stringi")
directory<-'/myPath/'
setwd(directory)


sampleFiles<-grep('_count.txt',list.files(directory), value=TRUE)
sampleFiles

sampleCondition<-c(rep("GFP",3), rep("M280L",2), rep("WT",2))


sampleTable<- data.frame(sampleName=stri_sub(sampleFiles,length(sampleNames),-11), fileName=sampleFiles, condition=sampleCondition)
sampleTable
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory = directory, design=~condition)
ddsHTSeq
dds<-DESeq(ddsHTSeq)
#dds$batch<-factor(c(rep("3", 4), rep("5", 4)))
nrow(dds)
dds <- dds[rowSums(counts(dds)) > 1, ]
nrow(dds)

##Control should be the first string ! -> consifering plotting
#dds$condition <- factor(dds$condition, levels=c("ND", "T2D"))
#res<- results(dds, contrast=c("condition","M280L","GFP"))
res_M280L_WT=results(dds, contrast=c("condition","M280L","WT"),independentFiltering = FALSE)
res_M280L_WT
res_M280L_GFP=results(dds, contrast=c("condition","M280L","GFP"),independentFiltering = FALSE)
res_M280L_GFP
res_WT_GFP=results(dds, contrast=c("condition","WT","GFP"),independentFiltering = FALSE)
res_WT_GFP

resultsNames(dds)
mcols(res_M280L_GFP,use.names=TRUE)
resLFC_M280L_GFP<- lfcShrink(dds, res=res_M280L_GFP, contrast=c("condition","M280L","GFP"))
resLFC_M280L_WT<- lfcShrink(dds, res=res_M280L_WT, contrast=c("condition","M280L","WT"))
resLFC_WT_GFP<- lfcShrink(dds, res=res_WT_GFP, contrast=c("condition","WT","GFP"))
#resLFC
res_M280L_GFP_Ordered<- res_M280L_GFP[order(res_M280L_GFP$pvalue), ]
res_M280L_WT_Ordered<- res_M280L_WT[order(res_M280L_WT$pvalue), ]
res_WT_GFP_Ordered<- res_WT_GFP[order(res_WT_GFP$pvalue), ]

summary(res_M280L_GFP)
summary(res_M280L_WT)

sum(resLFC_M280L_GFP$padj < 0.1, na.rm=TRUE)
sum(resLFC_M280L_WT$padj < 0.1, na.rm=TRUE)
sum(resLFC_WT_GFP$padj < 0.1, na.rm=TRUE)

# MA plot 
#plotMA(res,ylim=c(-2,2))

# MA plot in log2 fold change and interative mode

#M280L_GFP
pdf("M280L_GFP_after_MAplot.pdf")
plotMA(res_M280L_GFP, ylim = c(-15,15), main="M280L_GFP_MA plot")
topGene <- rownames(res_M280L_GFP)[which.min(res_M280L_GFP$padj)]
with(res_M280L_GFP[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})
dev.off()

#M280L_WT
pdf("M280L_WT_after_MAplot.pdf")
plotMA(res_M280L_WT, ylim = c(-10,10), main="M280L_WT_MA plot")
topGene <- rownames(res_M280L_WT)[which.min(res_M280L_WT$padj)]
with(res_M280L_WT[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})
dev.off()


# WT_GFP
pdf("WT_GFP_after_MAplot.pdf")
plotMA(res_WT_GFP, ylim = c(-15,15), main="WT_GFP_MA plot")
topGene <- rownames(res_WT_GFP)[which.min(res_WT_GFP$padj)]
with(res_WT_GFP[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})
dev.off()


# plot counts
pdf("M280L_WT_after_topGene_countPlot.pdf")
topGene <- rownames(res_M280L_WT)[which.min(res_M280L_WT$padj)]
plotCounts(dds, gene=topGene, intgroup="condition")

d <- plotCounts(dds, gene=which.min(res_M280L_WT$padj), intgroup="condition", 
                returnData=TRUE,main="M280L_WT_after_topGene_countPlot")
dev.off()

pdf("M280L_GFP_after_topGene_countPlot.pdf")
topGene <- rownames(res_M280L_GFP)[which.min(res_M280L_GFP$padj)]
plotCounts(dds, gene=topGene, intgroup="condition")

d <- plotCounts(dds, gene=which.min(res_M280L_GFP$padj), intgroup="condition", 
                returnData=TRUE,main="M280L_GFP_after_topGene_countPlot")
dev.off()

pdf("WT_GFP_after_topGene_countPlot.pdf")
topGene <- rownames(res_WT_GFP)[which.min(res_WT_GFP$padj)]
plotCounts(dds, gene=topGene, intgroup="condition")

d <- plotCounts(dds, gene=which.min(res_WT_GFP$padj), intgroup="condition", 
                returnData=TRUE,main="WT_GFP_after_topGene_countPlot")
dev.off()



library("ggplot2")
#ggplot(d, aes(x=condition, y=count)) + 
#  geom_point(position=position_jitter(w=0.1,h=0)) + 
#  scale_y_log10(breaks=c(25,100,400))

rld <- rlog(dds, blind=TRUE)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vsd.fast <- vst(dds, blind=TRUE)
head(assay(rld), 3)

# Effect of the transformation 
# this gives log2(n + 1)
ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))
meanSdPlot(assay(rld))
meanSdPlot(assay(vsd))


# Heatmaps for comparing various transformation methods
library("pheatmap")
library("grid")
library("data.table")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,"condition"])
colnames(df)[1] <- "condition"
(setattr(df, "row.names",c("GFP_1","GFP_2","GFP_3","M280L_1","M280L_2","WT_2","WT_3")))

pdf("M280L_GFP_after_ntd_heatmap.pdf")
mat <- assay(ntd)[head(order(res_M280L_GFP$padj),60),]
mat <- mat - rowMeans(mat)
df <- as.data.frame (colData(ntd)[,c("condition")])
colnames(df)[1] <- "condition"
(setattr(df, "row.names",c("GFP_1","GFP_2","GFP_3","M280L_1","M280L_2","WT_2","WT_3")))
pheatmap(mat,  annotation_col=df)
dev.off()
pdf("M280L_GFP_after_rld_heatmap.pdf")
#pheatmap(assay(rld)[select,], cluster_rows=TRUE, show_rownames=TRUE,
#         cluster_cols=TRUE, annotation_col=df)
mat <- assay(rld)[head(order(res_M280L_GFP$padj),60),]
mat <- mat - rowMeans(mat)
df <- as.data.frame (colData(rld)[,c("condition")])
colnames(df)[1] <- "condition"
(setattr(df, "row.names",c("GFP_1","GFP_2","GFP_3","M280L_1","M280L_2","WT_2","WT_3")))
pheatmap(mat,  annotation_col=df)
dev.off()
pdf("M280L_GFP_after_vsd_heatmap.pdf")
#pheatmap(assay(vsd)[select,], cluster_rows=TRUE, show_rownames=TRUE,
#         cluster_cols=TRUE, annotation_col=df)
mat <- assay(vsd)[head(order(res_M280L_GFP$padj),60),]
mat <- mat - rowMeans(mat)
df <- as.data.frame (colData(vsd)[,c("condition")])
colnames(df)[1] <- "condition"
(setattr(df, "row.names",c("GFP_1","GFP_2","GFP_3","M280L_1","M280L_2","WT_2","WT_3")))
pheatmap(mat,  annotation_col=df)
dev.off()


# Heatmap of the sample-to-sample distances
pdf("vsd_after_sampleDistance.pdf")
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

# sampleDists <- dist(t(assay(rld)))
# library("RColorBrewer")
# sampleDistMatrix <- as.matrix(sampleDists)
# rownames(sampleDistMatrix) <- paste(rld$condition, rld$type, sep="-")
# colnames(sampleDistMatrix) <- paste(rld$condition, rld$type, sep="-")
# colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
# pheatmap(sampleDistMatrix,
#          clustering_distance_rows=sampleDists,
#          clustering_distance_cols=sampleDists,
#          col=colors)

# sampleDists <- dist(t(assay(ntd)))
# library("RColorBrewer")
# sampleDistMatrix <- as.matrix(sampleDists)
# rownames(sampleDistMatrix) <- paste(ntd$condition, ntd$type, sep="-")
# colnames(sampleDistMatrix) <- paste(rnt$condition, ntd$type, sep="-")
# colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
# pheatmap(sampleDistMatrix,
#          clustering_distance_rows=sampleDists,
#          clustering_distance_cols=sampleDists,
#          col=colors)


# PCA
pdf("vsd_after_PCA.pdf")
plotPCA(vsd, intgroup=c("condition"))
# customized PCA by ggplot2
library("ggrepel")
pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + 
  geom_text_repel(aes(label=rownames(pcaData), size=3))
dev.off()

#annotate gene names
library("AnnotationDbi")
library("org.Hs.eg.db")
columns(org.Hs.eg.db)
res_M280L_GFP$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(res_M280L_GFP),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res_M280L_GFP$entrez <- mapIds(org.Hs.eg.db,
                     keys=row.names(res_M280L_GFP),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
res_M280L_WT$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(res_M280L_WT),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res_M280L_WT$entrez <- mapIds(org.Hs.eg.db,
                     keys=row.names(res_M280L_WT),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
res_WT_GFP$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(res_WT_GFP),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res_WT_GFP$entrez <- mapIds(org.Hs.eg.db,
                     keys=row.names(res_WT_GFP),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
#M280L_GFP
resOrdered <- res_M280L_GFP[order(res_M280L_GFP$pvalue),]
head(resOrdered)
resOrderedDF <- as.data.frame(resOrdered)[1:100, ]
write.csv(resOrderedDF, file = "M280L_GFP_after_DEG_results.csv")


library("ReportingTools")
htmlRep <- HTMLReport(shortName="report", title="M280L_GFP_after_report",
                      reportDirectory="./M280L_GFP_after_report")
publish(resOrderedDF, htmlRep)
url <- finish(htmlRep)
browseURL(url)

#M280L_WT
resOrdered <- res_M280L_WT[order(res_M280L_WT$pvalue),]
head(resOrdered)
resOrderedDF <- as.data.frame(resOrdered)[1:100, ]
write.csv(resOrderedDF, file = "M280L_WT_after_DEG_results.csv")


htmlRep <- HTMLReport(shortName="report", title="M280L_WT_after_report",
                      reportDirectory="./M280L_WT_after_report")
publish(resOrderedDF, htmlRep)
url <- finish(htmlRep)
browseURL(url)

#WT_GFP
resOrdered <- res_WT_GFP[order(res_WT_GFP$pvalue),]
head(resOrdered)
resOrderedDF <- as.data.frame(resOrdered)[1:100, ]
write.csv(resOrderedDF, file = "WT_GFP_after_DEG_results.csv")

htmlRep <- HTMLReport(shortName="report", title="WT_GFP_after_report",
                      reportDirectory="./WT_GFP_after_report")
publish(resOrderedDF, htmlRep)
url <- finish(htmlRep)
browseURL(url)

library("genefilter")
pdf("Integrative_heatmap_without_outlier_afterCluster.pdf",width=15,height=15)
topVarGenes <- head(order(rowVars(assay(vsd)),decreasing=TRUE),200)
#distCor <- function(x) as.dist(1-cor(t(x)))
#hclustAvg <- function(x) hclust(x, method="average")
library(gplots)
heatmap.2( assay(vsd)[ topVarGenes, ], scale="row",
           trace="none", dendrogram="row", margin=c(6,6),
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
           ColSideColors = c( GFP="gray", M280L="darkgreen", WT="orange" )[
             colData(vsd)$condition ] )
dev.off()


## GO enrichment analysis
library("goseq")
head(res_M280L_GFP)

all_genes <- row.names(res_M280L_GFP)
all_genes
DE_genes <- all_genes[which(res_M280L_GFP$padj<0.05)]
head(DE_genes)
DE_genes <- as.integer(res_M280L_GFP$padj<=0.05)
names(DE_genes) <- rownames(res_M280L_GFP)
DE_genes <- DE_genes[!is.na(DE_genes)]
#supportedOrganisms()[supportedOrganisms()$Genome=="hg19",]

pwf=nullp(DE_genes,"hg19","ensGene")
head(pwf)
GO.wall=goseq(pwf,"hg19","ensGene",use_genes_without_cat=TRUE)
head(GO.wall)

enriched.GO=GO.wall$category[p.adjust(GO.wall$over_represented_pvalue,method="BH")<.05]
head(enriched.GO)

library(GO.db)
for(go in enriched.GO[1:45]){ 
   print(GOTERM[[go]])
   cat("--------------------------------------\n")
  }


## KEGG analysis
library(pathview)
library(gage)
library(gageData)
data(kegg.sets.hs)
data(sigmet.idx.hs)
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
head(kegg.sets.hs, 3)
foldchanges = res_M280L_GFP$log2FoldChange
names(foldchanges) = res_M280L_GFP$entrez
head(foldchanges)

# Get the results
keggres = gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)

# Look at both up, down and statistics
lapply(keggres, head)

# Get the pathway
library(magrittr)
library(dplyr)
keggrespathways = data.frame(id=rownames(keggres$greater), keggres$greater) %>%
  tbl_df() %>%
  filter(row_number() <= 5) %>%
  .$id %>%
  as.character()
keggrespathways

# Get the IDs
keggresids = substr(keggrespathways, start=1, stop=8)
keggresids

# Define plotting function 
plot_pathway = function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa", new.signature=FALSE)
# plot multiple pathways (the plots will save in the disk and returns a throwaway list object)
detach("package:dplyr", unload=TRUE)
tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa"))

## GO analysis
# data(go.sets.hs)
# data(go.subs.hs)
# gobpsets = go.sets.hs[go.subs.hs$BP]
# gobpres = gage(foldchanges, gsets=gobpsets, same.dir=FALSE)
# lapply(gobpres, head)
