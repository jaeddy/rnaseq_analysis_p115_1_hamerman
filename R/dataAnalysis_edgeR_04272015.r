# read in data
rm(list=ls());gc()
library(limma)
library(biomaRt)
library(edgeR)

## whalen functions
######## E. Whalen's code
# computes normalization from libraries and then filters genes by having % samples with at least N counts and 
setUpDGEList <- function(countData, filterCount=1, filterPercentage=0.1)
{
	d <- DGEList(counts=countData)
	d <- calcNormFactors(d)

	keepRows <- rowSums(round(cpm(d$counts)) >= filterCount) >= filterPercentage*ncol(countData)
	print(table(keepRows))

	curDGE<-d[keepRows,]
	return(curDGE)
}

setUpSamples<-function(countData, designData, subsetSampleIndex)
{
curDesign<-designData[subsetSampleIndex,]
# remove any individuals that don't have a stim and non-stim sample
if (any(table(as.character(curDesign$donorID))==1))
{
remIndiv<-names(which(table(as.character(curDesign$donorID))==1))
curDesign<-curDesign[-which(curDesign$donorID %in% remIndiv),]
}
countData<-countData[,which(colnames(countData) %in% curDesign$Library.ID)]

return(list(countData=countData, curDesign=curDesign))
}
## end whalen function


cnts = read.csv(file="~/Desktop/BRI/RNAseq/P99/counts/combined_counts.csv",row.names=1); names(cnts) = gsub("\\.[0-9]+$","",names(cnts))
mets = read.csv(file="~/Desktop/BRI/RNAseq/P99/metrics/combined_metrics.csv",row.names=1); rownames(mets) = gsub("-[0-9]+$","",rownames(mets))

pdf(file="~/Desktop/BRI/RNAseq/P99/results/P99QCanalysisplots.pdf")
#par(mfrow=c(2,2))
x = mets$MEDIAN_CV_COVERAGE
y = mets$fastq_total_reads
z = mets$UNPAIRED_READS_EXAMINED / y
plot(x,y, xlim=c(0,max(max(x),1)),ylim= c(0,max(y)), pch=16,col="slategrey",cex=.5, xlab="Med CV Coverage", ylab="total reads")
plot(x,z, xlim=c(0,max(max(x),1)), pch=16,col="slategrey",cex=.5, xlab="Med CV Coverage", ylab="% aligned")
plot(y,z, pch=16,col="slategrey",cex=.5, xlab="total reads", ylab="% aligned")
plot(y,z, pch=16,col="slategrey",cex=.5, xlab="total reads", ylab="% aligned",type="n", main="Median CV Coverage")
text(y,z,labels=signif(x,2),cex=.65)


retDGE = setUpDGEList(countData=cnts, filterCount=50, filterPercentage=0.05)


designData     = read.csv(file="~/Desktop/BRI/RNAseq/P99/P99-1SampleSummary.csv")
designSub      = designData[match(colnames(retDGE$counts), designData$lib),]

plotMDS(retDGE$counts,cex=.5, main="MDS, lib labeled",labels=designSub$donorId) # may have to kick out one individual
dev.off()

des          = data.frame(donorID = as.factor(designSub$donorId))
des$cellType = designSub$studyGroup

#### method 1: voom with quality weights and use donor as a random effect
oneMM      <- model.matrix(~cellType, data=des)
DMvNSvoom1 <- voomWithQualityWeights(retDGE, design=oneMM, plot=TRUE)

# just looking at the voom object
names(DMvNSvoom1)
dim(DMvNSvoom1$E)
dim(DMvNSvoom1$weights)
DMvNSvoom1$E[1:5,1:5]
DMvNSvoom1$weights[1:5,1:5]
DMvNSvoom1$sample.weights[1:5]

M  = apply(DMvNSvoom1$E,1,mean)
SD = apply(DMvNSvoom1$E,1,sd) 
CV = SD / M

temp = DMvNSvoom1$E
temp = temp[CV > .5 & M > 4,]

readM = signif(y /1000000,2)

pdf(file="~/Desktop/BRI/RNAseq/P99/results/P99analysisplots.pdf")
plot(hclust(as.dist(1-cor(temp,method="spearman"))),labels = paste(des$donorID, des$cellType,readM))
plot(hclust(as.dist(1-cor(temp,method="kendall"))),labels = paste(des$donorID, des$cellType,readM))
#plot(hclust(dist(temp)),labels = paste(des$donorID, des$cellType,readM))

plot(density(DMvNSvoom1$E[,1]))
mus = c()
for(i in 1:ncol(DMvNSvoom1$E))
{
	dens = density(DMvNSvoom1$E[,i])
	lines(dens,col=i)
	mus = c(mus,mean(dens$y[dens$x > 6 & dens$x < 8]))
}
labs = designSub$lib[order(mus)]
cols = (1:10)[order(mus)]
text(11,seq(.1,.3,by =.2/10),labs,col=cols,cex=.65)

cells = as.numeric(gsub(".*_","",designSub[,1]))
totCnts =  apply(cnts,2,sum)
plot(cells, totCnts,pch=16,col="slategrey")

plot(cells, totCnts,type="n")
text(cells,totCnts,designSub$lib,cex=.6)

out = as.data.frame(DMvNSvoom1$E)

gc()
mart     = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
ens2Gene = getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values=rownames(out), mart=mart)
ens2Gene = ens2Gene[match(rownames(out), ens2Gene$ensembl_gene_id),]
out$geneSymbol = ens2Gene[,2]

allGene = getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values=rownames(cnts), mart=mart)
allGene = allGene[match(rownames(cnts), allGene$ensembl_gene_id),]


aGenes  = read.csv("~/Desktop/BRI/RNAseq/P99/amedeeGenesOfInterest.csv",header=F)


temp = out[ out$geneSymbol %in% aGenes[,1],1:10]; rownames(temp) = out[ out$geneSymbol %in% aGenes[,1],11,]
rowLabs = out[ out$geneSymbol %in% aGenes[,1],11,]
colLabs = paste(des$donorID, des$cellType)
heatmap(as.matrix(temp),labRow=rowLabs,labCol=colLabs,margins=c(12,5))
#heatmap(as.matrix(temp),distfun=function(x){as.dist(1-cor(t(x)))},labCol=colLabs,margins=c(12,5) )
# IL13RA2 is not in the filtered cnts with PALX etc

#plot(hclust(as.dist(1-cor(t(temp)))))
plot(hclust(dist(temp)))

dev.off()

write.csv(out,"~/Desktop/BRI/RNAseq/P99/expressedGenes.csv")























