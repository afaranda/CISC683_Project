################################################################################
# File: NormalizationFunctions.R											   #
# Purpose: Implement various normalization strategies as functions. Each       #
#          function takes a count matrix as input and returns a matrix of      #
#          normalized counts                                                   #
# Created: April 30, 2019 													   #
# Author: Adam Faranda														   #
################################################################################
library(edgeR)

# Get the position(s) (column number) of (a) column(s) by name from a data frame
colNum<-function(df, n){
	positions<-c()
	for(i in n){
		positions<-c(
			positions,
			grep(i, names(df))
		)
	}
	positions
}

# Basic Counts Per Million -- for each sample, divide each count by the sample's
# total library size
basicCPM<-function(x, idCol=1){
	idCol<-ifelse(is.character(idCol), colNum(x, idCol), idCol)
	cs<-colSums(x[,setdiff(1:ncol(x), idCol)])
	for (i in setdiff(1:ncol(x), idCol)){
		x[,i]<-(x[,i] / cs[i-1]) * 1000000
	}
	x
} 

# Wrapper for EdgeR's cpm function
edgeRcpm<-function(x, idCol=1, prior=2){
	idCol<-ifelse(is.character(idCol), colNum(x, idCol), idCol )
	row.names(x)<-x[,idCol]
	
	y<-DGEList(
		counts = x[,setdiff(1:ncol(x), idCol)],
		genes = row.names(x)
	)
	y<-calcNormFactors(y)
	x<-as.data.frame(2^cpm(y, prior.count=prior, log=T))
	print(head(x))
	x['ID']<-row.names(x)
	row.names(x)<-1:nrow(x)
	n<-ncol(x)
	x[,c(n, 1:(n-1))]
}



x<-data.frame(
	ID=c('A1','A2','A3','A4','A5','A6','A7','A8','A9','A10'),
	Samp1=c(100,100,100,100,100,100,100,100,100,100),
	Samp2=c(1000,800,200,400,600,500,500,100,900,0),
	Samp3=c(1,1,1,1,1,1,1,1,1,1),
	Samp4=c(2400,1600,2400,800,3200,400,400,800,3200,800),
	stringsAsFactors=F
)




