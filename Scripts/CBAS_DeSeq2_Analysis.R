#http://www-huber.embl.de/users/klaus/Teaching/DESeq2Predoc2014.html#preparing-count-matrices consulted 07.03.2016

library("EDASeq")
library("DESeq2")
library("pheatmap")
library("ggplot2")
library("gplots")
library("geneplotter")
library("RColorBrewer")
library("genefilter")

library("BiocParallel")
register(MulticoreParam(4))


input_counts<-"~/Repos/CBAS/Count_Matrix/CBAS_Bleaching_RSEM_Expression_Matrix.counts"

input_sample_information<-"~/Repos/CBAS/Count_Matrix/CBAS_Bleaching_RSEM_Expression_Matrix.info"


CBAS_DE<-read.csv(input_counts, head=T, sep="\t")

CBAS_DE_INFO<-read.csv(input_sample_information, head=T)

# Remove genes that have no counts over all samples
CBAS_DE <- CBAS_DE[(rowSums(CBAS_DE) > 0),]

CBAS_DeSeq<-DESeqDataSetFromMatrix(countData=CBAS_DE, colData=CBAS_DE_INFO, design=~condition)

#relevel the factors so that control is the reference
CBAS_DeSeq$condition<-relevel(CBAS_DeSeq$condition,ref="control")

#estimate size factors
#Thus, if all size factors are roughly equal to one, the libraries have been sequenced equally deeply.
CBAS_DeSeq<-estimateSizeFactors(CBAS_DeSeq)

#get rows with all non-zero counts
non_zero_rows<-apply(counts(CBAS_DeSeq), 1, function(x){all(x>0)})

#get all rows
all_rows<-apply(counts(CBAS_DeSeq), 1, function(x){all(x>=0)})

#number of non-zero rows:
sum(non_zero_rows)

#number of rows:
sum(all_rows)

#cummulative distribution of normalized counts for non-zero rows
multiecdf(counts(CBAS_DeSeq, normalized=T)[non_zero_rows, ], xlab="mean counts", xlim=c(0,1000))

#density of normalized counts
multidensity(counts(CBAS_DeSeq, normalized=T)[non_zero_rows, ], xlab="mean counts", xlim=c(0,1000))

#sample comparisons
MDPlot(counts(CBAS_DeSeq, normalized=T)[non_zero_rows, ],c(1,2), main=paste( colnames(CBAS_DeSeq)[1], "vs.", colnames(CBAS_DeSeq)[2]))
MDPlot(counts(CBAS_DeSeq, normalized=T)[non_zero_rows, ],c(1,3), main=paste( colnames(CBAS_DeSeq)[1], "vs.", colnames(CBAS_DeSeq)[3]))
MDPlot(counts(CBAS_DeSeq, normalized=T)[non_zero_rows, ],c(1,4), main=paste( colnames(CBAS_DeSeq)[1], "vs.", colnames(CBAS_DeSeq)[4]))
MDPlot(counts(CBAS_DeSeq, normalized=T)[non_zero_rows, ],c(1,5), main=paste( colnames(CBAS_DeSeq)[1], "vs.", colnames(CBAS_DeSeq)[5]))
MDPlot(counts(CBAS_DeSeq, normalized=T)[non_zero_rows, ],c(1,6), main=paste( colnames(CBAS_DeSeq)[1], "vs.", colnames(CBAS_DeSeq)[6]))
MDPlot(counts(CBAS_DeSeq, normalized=T)[non_zero_rows, ],c(1,7), main=paste( colnames(CBAS_DeSeq)[1], "vs.", colnames(CBAS_DeSeq)[7]))
MDPlot(counts(CBAS_DeSeq, normalized=T)[non_zero_rows, ],c(1,8), main=paste( colnames(CBAS_DeSeq)[1], "vs.", colnames(CBAS_DeSeq)[8]))
MDPlot(counts(CBAS_DeSeq, normalized=T)[non_zero_rows, ],c(2,3), main=paste( colnames(CBAS_DeSeq)[2], "vs.", colnames(CBAS_DeSeq)[3]))
MDPlot(counts(CBAS_DeSeq, normalized=T)[non_zero_rows, ],c(2,4), main=paste( colnames(CBAS_DeSeq)[2], "vs.", colnames(CBAS_DeSeq)[4]))
MDPlot(counts(CBAS_DeSeq, normalized=T)[non_zero_rows, ],c(2,5), main=paste( colnames(CBAS_DeSeq)[2], "vs.", colnames(CBAS_DeSeq)[5]))
MDPlot(counts(CBAS_DeSeq, normalized=T)[non_zero_rows, ],c(2,6), main=paste( colnames(CBAS_DeSeq)[2], "vs.", colnames(CBAS_DeSeq)[6]))
MDPlot(counts(CBAS_DeSeq, normalized=T)[non_zero_rows, ],c(2,7), main=paste( colnames(CBAS_DeSeq)[2], "vs.", colnames(CBAS_DeSeq)[7]))
MDPlot(counts(CBAS_DeSeq, normalized=T)[non_zero_rows, ],c(2,8), main=paste( colnames(CBAS_DeSeq)[2], "vs.", colnames(CBAS_DeSeq)[8]))
MDPlot(counts(CBAS_DeSeq, normalized=T)[non_zero_rows, ],c(2,9), main=paste( colnames(CBAS_DeSeq)[2], "vs.", colnames(CBAS_DeSeq)[9]))
MDPlot(counts(CBAS_DeSeq, normalized=T)[non_zero_rows, ],c(3,4), main=paste( colnames(CBAS_DeSeq)[3], "vs.", colnames(CBAS_DeSeq)[4]))
MDPlot(counts(CBAS_DeSeq, normalized=T)[non_zero_rows, ],c(3,5), main=paste( colnames(CBAS_DeSeq)[3], "vs.", colnames(CBAS_DeSeq)[5]))
MDPlot(counts(CBAS_DeSeq, normalized=T)[non_zero_rows, ],c(3,6), main=paste( colnames(CBAS_DeSeq)[3], "vs.", colnames(CBAS_DeSeq)[6]))
MDPlot(counts(CBAS_DeSeq, normalized=T)[non_zero_rows, ],c(3,7), main=paste( colnames(CBAS_DeSeq)[3], "vs.", colnames(CBAS_DeSeq)[7]))
MDPlot(counts(CBAS_DeSeq, normalized=T)[non_zero_rows, ],c(3,8), main=paste( colnames(CBAS_DeSeq)[3], "vs.", colnames(CBAS_DeSeq)[8]))
MDPlot(counts(CBAS_DeSeq, normalized=T)[non_zero_rows, ],c(3,9), main=paste( colnames(CBAS_DeSeq)[3], "vs.", colnames(CBAS_DeSeq)[9]))
MDPlot(counts(CBAS_DeSeq, normalized=T)[non_zero_rows, ],c(4,5), main=paste( colnames(CBAS_DeSeq)[4], "vs.", colnames(CBAS_DeSeq)[5]))
MDPlot(counts(CBAS_DeSeq, normalized=T)[non_zero_rows, ],c(4,6), main=paste( colnames(CBAS_DeSeq)[4], "vs.", colnames(CBAS_DeSeq)[6]))
MDPlot(counts(CBAS_DeSeq, normalized=T)[non_zero_rows, ],c(4,7), main=paste( colnames(CBAS_DeSeq)[4], "vs.", colnames(CBAS_DeSeq)[7]))
MDPlot(counts(CBAS_DeSeq, normalized=T)[non_zero_rows, ],c(4,8), main=paste( colnames(CBAS_DeSeq)[4], "vs.", colnames(CBAS_DeSeq)[8]))
MDPlot(counts(CBAS_DeSeq, normalized=T)[non_zero_rows, ],c(4,9), main=paste( colnames(CBAS_DeSeq)[4], "vs.", colnames(CBAS_DeSeq)[9]))
MDPlot(counts(CBAS_DeSeq, normalized=T)[non_zero_rows, ],c(5,6), main=paste( colnames(CBAS_DeSeq)[5], "vs.", colnames(CBAS_DeSeq)[6]))
MDPlot(counts(CBAS_DeSeq, normalized=T)[non_zero_rows, ],c(5,7), main=paste( colnames(CBAS_DeSeq)[5], "vs.", colnames(CBAS_DeSeq)[7]))
MDPlot(counts(CBAS_DeSeq, normalized=T)[non_zero_rows, ],c(5,8), main=paste( colnames(CBAS_DeSeq)[5], "vs.", colnames(CBAS_DeSeq)[8]))
MDPlot(counts(CBAS_DeSeq, normalized=T)[non_zero_rows, ],c(5,9), main=paste( colnames(CBAS_DeSeq)[5], "vs.", colnames(CBAS_DeSeq)[9]))
MDPlot(counts(CBAS_DeSeq, normalized=T)[non_zero_rows, ],c(6,7), main=paste( colnames(CBAS_DeSeq)[6], "vs.", colnames(CBAS_DeSeq)[7]))
MDPlot(counts(CBAS_DeSeq, normalized=T)[non_zero_rows, ],c(6,8), main=paste( colnames(CBAS_DeSeq)[6], "vs.", colnames(CBAS_DeSeq)[8]))
MDPlot(counts(CBAS_DeSeq, normalized=T)[non_zero_rows, ],c(6,9), main=paste( colnames(CBAS_DeSeq)[6], "vs.", colnames(CBAS_DeSeq)[9]))
MDPlot(counts(CBAS_DeSeq, normalized=T)[non_zero_rows, ],c(7,8), main=paste( colnames(CBAS_DeSeq)[7], "vs.", colnames(CBAS_DeSeq)[8]))
MDPlot(counts(CBAS_DeSeq, normalized=T)[non_zero_rows, ],c(7,9), main=paste( colnames(CBAS_DeSeq)[7], "vs.", colnames(CBAS_DeSeq)[9]))
MDPlot(counts(CBAS_DeSeq, normalized=T)[non_zero_rows, ],c(8,9), main=paste( colnames(CBAS_DeSeq)[8], "vs.", colnames(CBAS_DeSeq)[9]))


#transform the data to rlog
CBAS_DeSeq_rlog<-rlogTransformation(CBAS_DeSeq, blind = T)
#produce a heat plot using the transformed distances
CBAS_distances<-dist(t(assay(CBAS_DeSeq_rlog)))

heatmap.2(as.matrix(CBAS_distances), trace="none", col=rev(colorRampPalette(brewer.pal(9, "GnBu"))(100)))

#PCA of samples
plotPCA(CBAS_DeSeq_rlog, intgroup=c("condition"))

#remove outlier samples: Cauris_141021_6T
outliers=c("Cauris_20141021_6_T")
CBAS_DeSeq_NoOutliers<-CBAS_DeSeq[,!(colnames(CBAS_DeSeq) %in% outliers)]

#repeat above analyses
CBAS_DeSeq_NoOutliers_rlog<-rlogTransformation(CBAS_DeSeq_NoOutliers, blind = T)

#produce a heat plot using the transformed distances
CBAS_distances<-dist(t(assay(CBAS_DeSeq_NoOutliers_rlog)))

heatmap.2(as.matrix(CBAS_distances), trace="none", col=rev(colorRampPalette(brewer.pal(9, "GnBu"))(100)))

#PCA of samples
plotPCA(CBAS_DeSeq_NoOutliers_rlog, intgroup=c("condition"))

#heat map of most variable genes
CBAS_most_variable_genes<-head(order(rowVars(assay(CBAS_DeSeq_NoOutliers_rlog)), decreasing = T), 750)

heatmap.2(assay(CBAS_DeSeq_NoOutliers_rlog)[CBAS_most_variable_genes, ], scale="row", trace = "none", dendrogram = "column", col= colorRampPalette(rev(brewer.pal(9, "RdBu")))(255))


###
#
#Differential expression analysis
#
#####
CBAS_DeSeq_NoOutliers<-estimateSizeFactors(CBAS_DeSeq_NoOutliers)
CBAS_DeSeq_NoOutliers<-estimateDispersions(CBAS_DeSeq_NoOutliers)

plotDispEsts(CBAS_DeSeq_NoOutliers)

#wald test
CBAS_DeSeq_NoOutliers<-nbinomWaldTest(CBAS_DeSeq_NoOutliers)
CBAS_DeSeq_NoOutliers_Results<-results(CBAS_DeSeq_NoOutliers, pAdjustMethod = "BH")

#number of DE genes at 0.01 significance level
table(CBAS_DeSeq_NoOutliers_Results$padj < 0.01)

#number of DE genes at 0.01 significance level and log fold change > 2
table(CBAS_DeSeq_NoOutliers_Results$padj < 0.01 & abs(CBAS_DeSeq_NoOutliers_Results$log2FoldChange) > 2)

#number of overexpressed genes at 0.01 significance level and log fold change > 2
table(CBAS_DeSeq_NoOutliers_Results$padj < 0.01 & CBAS_DeSeq_NoOutliers_Results$log2FoldChange > 2)

#number of underexpressed genes at 0.01 significance level and log fold change < -2
table(CBAS_DeSeq_NoOutliers_Results$padj < 0.01 & CBAS_DeSeq_NoOutliers_Results$log2FoldChange < -2)

#plot of p values
hist(CBAS_DeSeq_NoOutliers_Results$pvalue, main = "Bleached vs. Control", xlab="p-values")

plotMA(CBAS_DeSeq_NoOutliers_Results, alpha=0.01)

#get significant genes
degs_names<-rownames(subset(CBAS_DeSeq_NoOutliers_Results, padj<0.01 & abs(log2FoldChange) > 2))
write.csv(degs_names, "~/Desktop/CBAS/DeSeq2/Results/all_degs_p001_2LFC.csv")

#get significantly overexpressed in treatment:
over_deg_names<-rownames(subset(CBAS_DeSeq_NoOutliers_Results, padj<0.01 & log2FoldChange > 2))
write.csv(over_deg_names, "~/Desktop/CBAS/DeSeq2/Results/over_degs_p001_2LFC.csv")
heatmap.2(log2(counts(CBAS_DeSeq_NoOutliers, normalized=T)[rownames(CBAS_DeSeq_NoOutliers) %in% over_deg_names, ]+1), scale="row", trace = "none", dendrogram = "column", col= colorRampPalette(rev(brewer.pal(9, "RdBu"))))

#get significantly underexpressed in treatment:
under_deg_names<-rownames(subset(CBAS_DeSeq_NoOutliers_Results, padj<0.01 & log2FoldChange < -2))
write.csv(under_deg_names, "~/Desktop/CBAS/DeSeq2/Results/under_degs_p001_2LFC.csv")
heatmap.2(log2(counts(CBAS_DeSeq_NoOutliers, normalized=T)[rownames(CBAS_DeSeq_NoOutliers) %in% under_deg_names, ]+1), scale="row", trace = "none", dendrogram = "column", col= colorRampPalette(rev(brewer.pal(9, "RdBu"))))

###############################
#
#
#Prepare data for Go analysis.
#
#ontologies: MF=Molecular Function, BP=Biochemical Process, CC=Cellular Component
#
##########################################################

#get base mean to get the background genes for GO analysis
BaseMean<-as.matrix(CBAS_DeSeq_NoOutliers_Results[, "baseMean", drop=F])

#find 10 genes with similar expression profile to serve as background in go enrichment analysis
GO_Background<-genefinder(BaseMean, degs_names, 8, method="manhattan")
GO_Over_Background<-genefinder(BaseMean, over_deg_names, 8, method="manhattan")
GO_Under_Background<-genefinder(BaseMean, under_deg_names, 8, method="manhattan")

#get the names of the selected genes
GO_Background_Gene_names<-rownames(BaseMean)[as.vector(sapply(GO_Background, function(x)x$indices))]
GO_Over_Background_Gene_names<-rownames(BaseMean)[as.vector(sapply(GO_Over_Background, function(x)x$indices))]
GO_Under_Background_Gene_names<-rownames(BaseMean)[as.vector(sapply(GO_Under_Background, function(x)x$indices))]

#remove degs from background if necessary
GO_Background_Gene_names<-setdiff(GO_Background_Gene_names, degs_names)
GO_Over_Background_Gene_names<-setdiff(GO_Over_Background_Gene_names, over_deg_names)
GO_Under_Background_Gene_names<-setdiff(GO_Under_Background_Gene_names, under_deg_names)

#save the names:
write.csv(GO_Background_Gene_names, "~/Desktop/CBAS/DeSeq2/Results/Back_degs_p001_2LFC.csv")
write.csv(GO_Over_Background_Gene_names, "~/Desktop/CBAS/DeSeq2/Results/Back_over_degs_p001_2LFC.csv")
write.csv(GO_Under_Background_Gene_names, "~/Desktop/CBAS/DeSeq2/Results/Back_under_degs_p001_2LFC.csv")

#number of background genes
length(GO_Background_Gene_names)
length(GO_Over_Background_Gene_names)
length(GO_Under_Background_Gene_names)

#plot the back and forground
multidensity( list(fore=log2(CBAS_DeSeq_NoOutliers_Results[degs_names, "baseMean"]),  back=log2(CBAS_DeSeq_NoOutliers_Results[GO_Background_Gene_names, "baseMean"])),  xlab="log2 mean counts", main = "Matching for enrichment analysis")
multidensity( list(fore=log2(CBAS_DeSeq_NoOutliers_Results[over_deg_names, "baseMean"]),  back=log2(CBAS_DeSeq_NoOutliers_Results[GO_Over_Background_Gene_names, "baseMean"])),  xlab="log2 mean counts", main = "Matching for enrichment analysis")
multidensity( list(fore=log2(CBAS_DeSeq_NoOutliers_Results[under_deg_names, "baseMean"]),  back=log2(CBAS_DeSeq_NoOutliers_Results[GO_Under_Background_Gene_names, "baseMean"])),  xlab="log2 mean counts", main = "Matching for enrichment analysis")
