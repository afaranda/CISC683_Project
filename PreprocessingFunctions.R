################################################################################
# File: NormalizationFunctions.R											   #
# Purpose: Implement various normalization strategies as functions. Each       #
#          function takes a count matrix as input and returns a matrix of      #
#          normalized counts                                                   #
# Created: April 30, 2019 													   #
# Author: Adam Faranda														   #
################################################################################
library(edgeR)
library(sva)

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

# Wrapper for EdgeR's cpm function, transforms to log2 to add prior count
# then back to cpm (prior only works when log=T).
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

# Simple variance filter -- take the top n or % rows of a matrix based on
# overall variance.  Assumes Column 1 is an ID column.
varianceFilter<-function(
df, idCol=1, threshold=10,
mode='absolute', decreasing=T

){
    idCol<-ifelse(is.character(idCol), colNum(df, idCol), idCol)
    df$TotalVar <-apply(df[,setdiff(1:ncol(x), idCol)], 1, var, na.rm=T)
    tvCol<-grep('TotalVar', names(df))
    if(!mode %in% c('absolute', 'percent')){
        print("threshold mode must be 'absolute' or 'percent'")
        return(NULL)
    }
    else if(mode == 'absolute'){
        threshold <- ifelse(threshold <=nrow(df), threshold, nrow(df))
        return(df[order(df$TotalVar, decreasing=decreasing),][1:threshold, -tvCol])
    }
    else if(mode == 'percent'){
        threshold<-ifelse(0 < threshold,
        ifelse(threshold <= 100, threshold, 100),
        100
        )
        n<-round((threshold * nrow(df))/100)
        return(df[order(df$TotalVar, decreasing=decreasing),][1:n, -tvCol])
    }
}




# Simple standard dev filter -- take the top n or % rows of a matrix based on
# overall standard deviation  Assumes Column 1 is an ID column.
stdevFilter<-function(
df, idCol=1, threshold=10,
mode='absolute', decreasing=T

){
    idCol<-ifelse(is.character(idCol), colNum(x, idCol), idCol)
    df$TotalVar <-apply(df[,setdiff(1:ncol(x), idCol)], 1, sd, na.rm=T)
    tvCol<-grep('TotalVar', names(df))
    if(!mode %in% c('absolute', 'percent')){
        print("threshold mode must be 'absolute' or 'percent'")
        return(NULL)
    }
    else if(mode == 'absolute'){
        threshold <- ifelse(threshold <=nrow(df), threshold, nrow(df))
        return(df[order(df$TotalVar, decreasing=decreasing),][1:threshold, -tvCol])
    }
    else if(mode == 'percent'){
        threshold<-ifelse(0 < threshold,
        ifelse(threshold <= 100, threshold, 100),
        100
        )
        n<-round((threshold * nrow(df))/100)
        return(df[order(df$TotalVar, decreasing=decreasing),][1:n, -tvCol])
    }
}


# Combat correct (sva) -- assumes all, or all but one columns are expression
# measurements -- apply combat batch correction to a data matrix
wrapCombat<-function(df, ft, groupCol=3, batchCol=4, idCol=1){
    if(idCol !=0 | !is.null(idCol)){
        idCol<-ifelse(is.character(idCol), colNum(x, idCol), idCol)
        idtype<-names(df)[idCol]
        samples<-setdiff(1:ncol(df), idCol)
        row.names(df)<-df[,idCol]
        df<-df[,samples]
    }
    
    batchCol<-ifelse(is.character(groupCol),
        batchCol,
        names(ft)[batchCol]
    )

    groupCol<-ifelse(is.character(groupCol),
        groupCol,
        names(ft)[groupCol]
    )
    
    
    ft[,groupCol]<-as.factor(ft[, groupCol])
    
    
    model<-model.matrix(formula(paste('~', groupCol)), ft)
    df<-ComBat(df, batch=ft[,batchCol], mod=model)
    df[idtype]<-row.names(df)
    row.names(df)<-1:nrow(df)
    n<-ncol(df)
    df<-df[,c(n, 1:n-1)]
    df
}

wrapCombat_intOnly<-function(df, ft, groupCol=3, batchCol=4, idCol=1){
    if(idCol !=0 | !is.null(idCol)){
        idCol<-ifelse(is.character(idCol), colNum(df, idCol), idCol)
        idtype<-names(df)[idCol]
        samples<-setdiff(1:ncol(df), idCol)
        row.names(df)<-df[,idCol]
        df<-df[,samples]
    }
    
    batchCol<-ifelse(is.character(groupCol),
        batchCol,
        names(ft)[batchCol]
    )

    groupCol<-ifelse(is.character(groupCol),
        groupCol,
        names(ft)[groupCol]
    )
    
    
    ft[,groupCol]<-as.factor(ft[, groupCol])
    
    
    model<-model.matrix(~1, data = ft)
    df<-ComBat(df, batch=ft[,batchCol], mod=model)
    df[idtype]<-row.names(df)
    row.names(df)<-1:nrow(df)
    n<-ncol(df)
    df<-df[,c(n, 1:n-1)]
    df
}

# ComBat introduces negative expression values, re-assign these to the
# minimum observed value.  This is required for log transformation
fixCombatNegatives<-function(df, idCol = 1){
	mpv<-as.numeric(
		as.matrix(
			df[,setdiff(1:ncol(df), idCol)]
		)
	)
	mpv<-min(mpv[mpv >0])
	for(i in setdiff(1:ncol(df), idCol)){
		df[df[,i]<0,i]<-mpv
	}
	df
}

# EdgeR Pairwise Contrast
edgeRPairwise<-function(
	df, ft, idCol=1, sampCol=7, group.dn=c(1,2,3), 
	group.up=c(4,5,6), groupCol=8
){
	print(group.dn)
	idCol<-ifelse(is.character(idCol), colNum(df, idCol), idCol)
	if(is.character(group.dn)){
		group.dn<-as.numeric(row.names((ft[ft[,sampCol]%in% group.dn,])))
		group.up<-as.numeric(row.names((ft[ft[,sampCol]%in% group.up,])))
	}
	
	if(is.numeric(group.dn)){
		print(group.dn)
		gr.dn<-ft[group.dn, sampCol]
		print(group.up)
		gr.up<-ft[group.up, sampCol]
	}
	
	row.names(df)<-df[,idCol]
	counts<-df[, setdiff(1:ncol(df), idCol)]
	samples<-c(gr.dn, gr.up)
	counts<-counts[,samples]
	
	y<-DGEList(
		counts=counts, 
		group=droplevels(ft[c(group.dn, group.up), groupCol]), 
		genes=row.names(df)	
	)
	# Filter DGEList: Remove genes where fewer than two samples have a cpm > 1
	keep <- rowSums(cpm(y) > 1) >= 2
	y.filter <-y[keep, ]
	y.filter$samples$lib.size <- colSums(y.filter$counts) # Fix library size after filtering
	
	# Calculate Normalization factors and dispersion estimates
	y.filter <-calcNormFactors(y.filter)
	y.filter <-estimateCommonDisp(y.filter, verbose=T)
	y.filter <-estimateTagwiseDisp(y.filter)
	
	# Calculate Differential Expression
	et<-exactTest(y.filter, pair=levels(y.filter$samples$group))
	degSet<-topTags(et, n=25000)@.Data[[1]]
	
	return(list(y.filter, degSet))
}


scaleCenterByRow<-function(mat){
	a<-apply(mat, 1, mean, na.rm=T)
	s<-apply(mat, 1, sd, na.rm=T)
	return((mat - a)/s)
}

x<-data.frame(
	ID=c('A1','A2','A3','A4','A5','A6','A7','A8','A9','A10'),
	Samp1=c(100,100,100,100,100,100,100,100,100,100),
	Samp2=c(1000,800,200,400,600,500,500,100,900,0),
	Samp3=c(1,1,1,1,1,1,1,1,1,1),
	Samp4=c(2400,1600,2400,800,3200,400,400,800,3200,800),
	stringsAsFactors=F
)




