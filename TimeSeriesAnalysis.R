################################################################################
# File: TimeSeriesAnalysis.R												       #
# Purpose: Identify gene clusters that correspond to temporal patterns.		   #
# Created: May 1, 2019 														   #
# Author: Adam Faranda														   #
################################################################################

############################ Setup Environment #################################
setwd('/Users/afaranda/Documents/CISC683_Project')
library(dplyr)
library(cluster)
library(reshape2)
# wd<-getwd()
source('BuildDataMatrix.R')
source('PreprocessingFunctions.R')
source('PrincipalComponents.R')
source('ClusteringFunctions.R')
######################### Apply Preprocessing Steps ############################

# Add Class attribute to feature definition table
ft<-raw[[1]]
ft$Class<-as.factor(paste("Hour",ft$Hours_PCS,sep=""))
ft$Class<-factor(
	ft$Class, 
	levels=unique(
		as.character(ft$Class[order(ft$Hours_PCS)])
	)
)


counts<-(raw[[2]][6:nrow(raw[[2]]),])      # Extract Raw Counts
ecpm<-edgeRcpm(counts)                     # Normalize using edgeR's TMM method
ecmb<-wrapCombat_intOnly(ecpm, ft)         # Correct for batch effects
ecmb<-fixCombatNegatives(ecmb)             # replace negative values with min +ve
ecpm.log<-log(ecmb[,2:19])				   # Apply log transformation Combat

############## Apply Variance Filters; Plot Principal Components ###############
varRanks <-c(10, 50, 100, 200)                # Try different variance filters
for( v in varRanks){
	print(v)
	ecpm.filter <-varianceFilter(ecpm, threshold=v)
	ecmb.filter <-varianceFilter(ecmb, threshold=v)
	
	# Plot results for TMM Normalized Data
	f1<-paste('ECPM_Samples_Top_', v,'_Ranked_Class.png')
	f2<-paste('ECPM_Samples_Top_', v,'_Ranked_Lab.png')
	png(f1, width=240, height=150)
		print(plotPrinComp(ecpm.filter, ft, groupCol=8, idCol=1))
	dev.off()
	png(f2, width=240, height=150)
		print(plotPrinComp(ecpm.filter, ft, groupCol=4, idCol=1))
	dev.off()

	# Plot results for batch adjusted TMM data
	f1<-paste('ECMB_Samples_Top_', v,'_Ranked_Class.png')
	f2<-paste('ECMB_Samples_Top_', v,'_Ranked_Lab.png')
	png(f1, width=240, height=150)
		print(plotPrinComp(ecmb.filter, ft, groupCol=8, idCol=1))
	dev.off()
	png(f2, width=240, height=150)
		print(plotPrinComp(ecmb.filter, ft, groupCol=4, idCol=1))
	dev.off()
}

############ Analyze Sample Clusters at desired Variance Threshold #############
distm <-c('euclidean', 'manhattan')			  # Try different distance methods
linkm <-c('complete', 'average', 'single')    # Try different linkage methods
trees <-c(1,2,3,4,5,6)                        # Different levels k
v =200
ecpm.filter<-varianceFilter(ecpm, threshold=v)
ecmb.filter<-varianceFilter(ecmb, threshold=v)

clustStats<-rbind(
	summarizeSampleClusters(
		data=ecpm, distm=distm, linkm=linkm, v=v, label='ecpm'
		
	),
	summarizeSampleClusters(
		data=ecmb, distm=distm, linkm=linkm, v=v, label='ecmb'
	)
)

write.csv(clustStats, 'Sample_Cluster_Statistics.csv')

################# Analyze Gene Clusters at Variance Threshold ####################
v = 200
ecmb.filter<-varianceFilter(ecmb, threshold=v)
row.names(ecmb.filter)<-ecmb.filter$ID
mat<-ecmb.filter[,2:ncol(ecpm.filter)]

h<-wrapHclust(log(mat), idCol=0, transpose=F, d.meth='manhattan', h.method='complete')
ktable<-tabulate_H_Clusters(h, ks=1:20)
ct<-reshapeClusterTable(log(mat), ktable, ft, k=20)

s<-as.data.frame(
	silhouette(cutree(h, 20), dist(log(mat), method='manhattan'))[,1:3]
)
s.means<-s %>% 
	group_by(cluster) %>%
	summarize(Mean_Silhouette=mean(sil_width))

cstat<-inner_join(
	ct %>% 
		group_by(Cluster) %>%
		summarize(Count = n()/18, Mean_logCPM=mean(value), STDev_logCPM = sd(value)),		
	s %>%
		group_by(cluster) %>%
		dplyr::rename(Cluster = cluster) %>%
		summarize(Mean_Silhouette=mean(sil_width)),
		by='Cluster'
)
write.csv(cstat, 'Gene_Cluster_Stats.csv')
			
cl<- cstat %>%
	filter(Count > 3, STDev_logCPM < 1)
cl<-cl$Cluster

bp<-ggplot(data=ct[ct$Cluster %in% cl,], mapping=aes(x=Class, y=value, color=Class)) + 
	geom_boxplot() + facet_grid( . ~Cluster)

png('Gene_Cluster_Profiles.png', width=960, height=280)
	print(bp)
dev.off()

