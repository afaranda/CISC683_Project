################################################################################
# File: TimeSeriesAnalysis.R												   #
# Purpose: Identify gene clusters that correspond to temporal patterns.		   #
# Created: May 1, 2019 														   #
# Author: Adam Faranda														   #
################################################################################

############################ Setup Environment #################################
setwd('/Users/afaranda/Documents/CISC683_Project')
library(dplyr)
# wd<-getwd()
source('BuildDataMatrix.R')
source('PreprocessingFunctions.R')
source('PrincipalComponents.R')
source('ClusteringFunctions.R')

############## Apply Preprocessing Steps and Evaluate Results ##################

ft<-raw[[1]]
ft$Class<-as.factor(paste("Hour",ft$Hours_PCS,sep=""))
counts<-(raw[[2]][6:nrow(raw[[2]]),])
# png("RawCounts.png", width=240, height=150)
	# plotPrinComp(counts, raw[[1]], groupCol=4)
# dev.off()

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
# for(i in 2:ncol(ecmb.filter)){
	# ecmb.filter[ecmb.filter[,i] < 0, i] <-1
# }
# pheatmap(log(ecmb.filter[2:ncol(ecmb.filter)], 2))

# Try several different methods -- get a range of values
# Tabulate Cluster statistics for each method attempted. 
h<-wrapHclust(ecmb.filter, ft,d.meth='minkowski', h.method = "median", d.p=0.2)
plotHclust(h, ft)
tabulate_H_Clusters(h, ft)


# Try some Bi-Clustering -- not looking good with variance filtered gene
mat<-as.matrix(log(ecmb.filter[,2:ncol(ecmb.filter)]))

mat.dbi<-mat[,1:9]
mat
res<-biclust(t(mat.dbi), method=BCCC(), delta=2.2, alpha=1.8)

res<-biclust(mat.dbi, method=BCPlaid())
res<-biclust(t(mat), method=BCSpectral())


# Try to get differentially expressed genes and use that to guide analysis
# deg_master<-data.frame(
	# Lab = character(),
	# DownGroup = character(),
	# UpGroup = character(),
	# genes = character(),
	# logFC = numeric(),
	# logCPM = numeric(),
	# PValue = numeric(),
	# FDR = numeric()
# )
# deg_master[1,]<-c("DROP", NA, NA, NA, NA, NA, NA, NA)
dbi.contrasts <-list(c("Hour0", "Hour24"), c("Hour0", "Hour48"), c("Hour24", "Hour48"))
dna.contrasts <-list(c("Hour0", "Hour6"), c("Hour0", "Hour24"), c("Hour6", "Hour24"))
for( i in unique(ft$Seq_Lab)){
	if(i == 'DNA'){
		for(j in dna.contrasts){
			g1<-ft[ft$Seq_Lab == i & ft$Class == j[1],'Sample_Number']
			g2<-ft[ft$Seq_Lab == i & ft$Class == j[2],'Sample_Number']
			deg<-edgeRPairwise(raw[[2]][6:nrow(raw[[2]]),], ft, group.dn=g1, group.up=g2)[[2]]
			deg$Lab <- i
			deg$DownGroup = j[1]
			deg$UpGroup = j[2]
			if(!exists('deg_master')){ rm(deg_master); deg_master<-deg}
			else{deg_master<-rbind(deg_master, deg)}
		}
	}
	else{
		for(j in dbi.contrasts){
			g1<-ft[ft$Seq_Lab == i & ft$Class == j[1],'Sample_Number']
			g2<-ft[ft$Seq_Lab == i & ft$Class == j[2],'Sample_Number']
			deg<-edgeRPairwise(raw[[2]][6:nrow(raw[[2]]),], ft, group.dn=g1, group.up=g2)[[2]]
			deg$Lab <- i
			deg$DownGroup = j[1]
			deg$UpGroup = j[2]
			if(!exists('deg_master')){rm(deg_master); deg_master<-deg}
			else{deg_master<-rbind(deg_master, deg)}
		
		}
	}
}


# Get DBI Degs -- Try Clustering Biclustering -- not much better
deg_master %>%
	filter(FDR < 0.05, abs(logFC) > 0.5) %>%
	group_by(Lab, DownGroup, UpGroup) %>%
	summarize(n(), Upregulated=sum(logFC > 1), Downregulated=sum(logFC < 1))

# Get Genes that are significant in DBI based on pairwise contrasts
dbi.gl<-inner_join(
	deg_master %>% 
	filter(FDR < 0.05, abs(logFC) > 0.6) %>%
	filter(DownGroup == 'Hour0', UpGroup == 'Hour24', Lab=='DBI') %>%
	select(genes, logFC) %>%
	rename( H0vsH24_logFC= 'logFC'),
	
	deg_master %>% 
	filter(FDR < 0.05, abs(logFC) > 0.6) %>%
	filter(DownGroup == 'Hour24', UpGroup == 'Hour48', Lab=='DBI') %>%
	select(genes, logFC) %>%
	rename(H24vsH48_logFC = 'logFC'),
	by='genes'
)

ecpm.filter<-ecpm[ecpm$ID %in% dgl$genes, 1:10]
row.names(ecpm.filter)<-ecpm.filter$ID
ecpm.filter<-ecpm.filter[,2:10]
mat<-as.matrix(ecpm.filter)

d<-dist(mat)

getClusterStat(mat, tabulate_k_means(mat, idCol=0, transpose=F, ks=1:5), k=1)



