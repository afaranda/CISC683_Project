################################################################################
# File: ClusteringFunctions.R												   #
# Purpose: Wrappers for various clustering algorithms, and plotting functions. #
# Created: May 1, 2019														   #
# Author: Adam Faranda														   #
################################################################################
library(dplyr)
library(ggfortify)

wrapHclust<-function(df, idCol=1,
	transpose=T, d.meth="euclidean", d.p = 2,
	h.method = "average"
 ){
 	if(!is.null(idCol) & idCol !=0){
		idCol<-ifelse(is.character(idCol), colNum(x, idCol), idCol)
		row.names(df)<-df[,idCol]
		df<-df[,setdiff(1:ncol(df), idCol)]
	}
	if(transpose){df<-t(df)}
	d<-dist(df, method=d.meth, p=d.p)
	h<-hclust(d, method=h.method)
}

# Plot dendrograms for heirarchical clusters -- only use when clustering samples
plotHclust<-function(h, ft, sampleCol=7, labelCol=8, colorCol=4, main=''){
	labColor <-ifelse(ft[h$order,colorCol] == 'DNA', "red", "black")
	x<-1:nrow(ft)
	lab <- ft[h$order, labelCol]
	print(lab)
	plot(
		h, labels = F,
		main = main,
		xlab = '',
		sub = '', 
		hang = -1
	)
	#print(labColor)
	text(x=x, labels = lab, col=labColor, srt = 90, xpd=NA, adj=c(1.2,0.5))
}


tabulate_H_Clusters<-function(h, ks=c(1:5) ){
	hc_table <- data.frame(
		Sample = h$labels
	)
	for (k in ks){
		if(k <= nrow(hc_table)){
			c=data.frame(c=cutree(h, k))
			c$Sample = row.names(c)
			hc_table<-merge(
				hc_table, c,
				by = 'Sample'
			)
			name<-paste("k_eq_", k, sep='')
			names(hc_table)[grep('c', names(hc_table))]<-name
		}
	}
	row.names(hc_table)<-hc_table$Sample
	hc_table<-hc_table[,setdiff(1:ncol(hc_table), grep('Sample', names(hc_table)))]
	hc_table
}


tabulate_k_means<-function(df, idCol=1, ks=c(1:5), transpose=T, method='Hartigan-Wong'){
	if(!is.null(idCol) & idCol != 0){
		idCol<-ifelse(is.character(idCol), colNum(x, idCol), idCol)
		row.names(df)<-df[,idCol]
		df<-df[,setdiff(1:ncol(df), idCol)]
	}
	if(transpose){df<-t(df)}
	hc_table <- data.frame(
		Sample = row.names(df)
	)
	for (k in ks){
		if(k <= nrow(hc_table)){
			c=data.frame(c=kmeans(df, nstart=nrow(hc_table), centers=k, algorithm=method )[[1]])
			c$Sample = row.names(c)
			hc_table<-merge(
				hc_table, c,
				by = 'Sample'
			)
			name<-paste("k_eq_", k, sep='')
			names(hc_table)[grep('c', names(hc_table))]<-name
		}
	}
	row.names(hc_table)<-hc_table$Sample
	hc_table<-hc_table[,setdiff(1:ncol(hc_table), grep('Sample', names(hc_table)))]
	hc_table
}
# Requires dplyr -- cant have variable number of columns
# getClusterMeans<-function(mat, ktable, k){
	# mat<-as.data.frame(mat)
	# ktable<-as.data.frame(ktable[k])
	# mat$ID<-row.names(mat)
	# ktable$ID<-row.names(ktable)
	
	# print(names(ktable)[1])
	# print(head(ktable))
	# x<-inner_join(mat, ktable, by='ID') %>%
		# group_by_(names(ktable)[1])	%>%
		# summarize( n(),
			# mean(Sample_1), mean(Sample_2), mean(Sample_3), 
			# mean(Sample_4), mean(Sample_5), mean(Sample_6), 
			# mean(Sample_7), mean(Sample_8), mean(Sample_9)
		# )
	# as.data.frame(x)
# }

apfunc<-function(x, applyFun='mean'){
	if(applyFun == 'mean'){
		return(mean(x, na.rm=T))
	} else if(applyFun == 'sd'){
		return(sd(x, na.rm=T))
	} else if(applyFun == 'max'){
		return(max(x, na.rm=T))
	} else if(applyFun == 'min'){
		return(min(x, na.rm=T))
	} else if(applyFun == 'length'){
		return(length(x))
	}
}

getClusterStat<-function(mat, ktable, k, applyFun='mean', transpose=F){
	if(transpose){mat<-t(mat)}
	mat<-as.data.frame(mat)
	ktable<-as.data.frame(ktable[k])
	mat$ID<-row.names(mat)
	kname<-as.character(names(ktable))
	ktable$ID<-row.names(ktable)
	
	x<-merge(mat, ktable, by='ID')
	row.names(x)<-x$ID
	x<-x[2:ncol(x)]
	
	for (i in unique(x[,kname])){
		y<-x[x[,kname] == i,setdiff(names(x), kname)]
		m<-as.list(apply(y, 2, apfunc, applyFun=applyFun))
		m[['count']]<-nrow(y)
		m[['cluster']]<-i
		#nc<-grep('kname', names(x))
		#x<-x[, c(nc, setdiff(1:ncol(x), nc))]
		if(exists('clust_table')){clust_table<-rbind(clust_table, m)}
		else{clust_table<-as.data.frame(m)}
	}
	clust_table<-as.data.frame(clust_table)
	row.names(clust_table)<-clust_table$cluster
	clust_table
}


plotCluster<-function(mat, ktable, k, c=1, ft){
	# Get cluster statistics 
	clustMeans<-getClusterStat(mat, ktable, k, applyFun='mean')
	clustSDs<-getClusterStat(mat, ktable, k, applyFun='sd')
	clustSize<-getClusterStat(mat, ktable, k, applyFun='length')
	
	# Drop Count, Cluster columns from statistics tables
	clustMeans<-clustMeans[,setdiff(names(clustMeans), c('cluster', 'count'))]
	clustSDs<-clustSDs[,setdiff(names(clustSDs), c('cluster', 'count'))]
	clustSize<-clustSize[,setdiff(names(clustSize), c('cluster', 'count'))]

	# Transpose 
	Mean<-as.data.frame(t(clustMeans))[c]
	names(Mean)<-'Mean'
	
	StDev<-as.data.frame(t(clustSDs))[c]
	names(StDev)<-'StDev'
	
	Size<-as.data.frame(t(clustSize))[c]
	names(Size)<-'Size'
	
	df<-merge(Mean, StDev, by=0)
	row.names(df)<-df$Row.names
	df<-df[,setdiff(names(df), 'Row.names')]
	
	df<-merge(df, Size, by=0 )
	row.names(df)<-df$Row.names
	df<-df[,setdiff(names(df), 'Row.names')]
	
	row.names(ft) <-ft$Sample_Number
	df$Group <- ft[row.names(df),'Sample_Number' ]
	df
	#boxplot(mean ~ )
	
	
	

}









