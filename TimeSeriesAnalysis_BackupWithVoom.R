################################################################################
# File: TimeSeriesAnalysis.R												   #
# Purpose: Identify gene clusters that correspond to temporal patterns.		   #
# Created: May 1, 2019 														   #
# Author: Adam Faranda														   #
################################################################################

############################ Setup Environment #################################
setwd('/Users/afaranda/Documents/CISC683_Project')
# wd<-getwd()
source('BuildDataMatrix.R')
source('PreprocessingFunctions.R')
source('PrincipalComponents.R')
source('ClusteringFunctions.R')

############## Apply Preprocessing Steps and Evaluate Results ##################

ft<-raw[[1]]
ft$Class<-as.factor(paste("Hour",ft$Hours_PCS,sep=""))
counts<-(raw[[2]][6:nrow(raw[[2]]),])
png("RawCounts.png", width=240, height=150)
	plotPrinComp(counts, raw[[1]], groupCol=4)
dev.off()

# Normalize using edgeR's TMM method
ecpm<-edgeRcpm(counts)
png("TMM_Normalized.png", width=240, height=150)
	plotPrinComp(ecpm, raw[[1]], groupCol=4)
dev.off()

# Filter by variance
ecpm.filter<-varianceFilter(ecpm, threshold=500)
png("TMM_Normalized_variance20.png", width=240, height=150)
	plotPrinComp(ecpm.filter, ft, groupCol=8)
dev.off()

# Combat batch effects -- shows some improvement
ecmb<-wrapCombat(ecpm, ft)
ecmb.filter<-varianceFilter(ecmb, threshold=500)
png("ComBat_variance500.png", width=240, height=150)
	plotPrinComp(ecmb.filter, ft, groupCol=8)
dev.off()

# Combat batch effects use intercept only 
ecmb<-wrapCombat_intOnly(ecpm, ft)
ecmb.filter<-varianceFilter(ecmb, threshold=500)
png("ComBat_IntOnly_variance500.png", width=240, height=150)
	plotPrinComp(ecmb.filter, ft, groupCol=8, legendTitle="Interval")
dev.off()

# Add 1 and log transform to facilitate heat map
for(i in 2:ncol(ecmb.filter)){
	ecmb.filter[ecmb.filter[,i] < 0, i] <-1
}
pheatmap(log(ecmb.filter[2:ncol(ecmb.filter)], 2))

# Try several different methods -- get a range of values
# Tabulate Cluster statistics for each method attempted. 
h<-wrapHclust(ecmb.filter, ft,d.meth='minkowski', h.method = "median", d.p=0.2)
plotHclust(h, ft)
tabulate_H_Clusters(h, ft)


# Try some Bi-Clustering -- not looking good
mat<-as.matrix(log(ecmb.filter[,2:ncol(ecmb.filter)]))

mat.dbi<-mat[,1:9]
mat
res<-biclust(t(mat.dbi), method=BCCC(), delta=2.2, alpha=1.8)

res<-biclust(mat.dbi, method=BCPlaid())
res<-biclust(t(mat), method=BCSpectral())


# Try to get differentially expressed genes and use that to guid analysis
#row.names(counts)<-counts$Ensembl -- the following is from RNASeq 123
#counts<-counts[, 2:ncol(counts)]
y<-DGEList(
	counts = counts[,1:9], 
	genes = row.names(counts)
)
keep = rowSums(cpm(y) > 1) >= 2
y.filter<-y[keep,]
y.filter$samples$lib.size<-colSums(y.filter$counts)
y.filter<-calcNormFactors(y.filter)
y.filter<-estimateCommonDisp(y.filter)
y.filter<-estimateTagwiseDisp(y.filter)
y.filter$samples$group<-droplevels(ft$Class[1:9])
design<-model.matrix(~0+group, data=y.filter$samples)
colnames(design)<-gsub('group', '',colnames(design))
cont.matrix<-makeContrasts(
	Hour0vsHour48 = Hour48 - Hour0,
	Hour24vsHour48 = Hour48 - Hour24,
	Hour0vsHour24 = Hour24 - Hour0, 
	levels = colnames(design)
)

v<-voom(y.filter, design)
vfit<-lmFit(v)
vfit.cont<-contrasts.fit(vfit, cont.matrix)
efit<-eBayes(vfit)





