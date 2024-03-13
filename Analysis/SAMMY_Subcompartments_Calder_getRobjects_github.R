###  SAMMY seq analyses paper - October 2022 ###
## Call subcompartment: re-implementation of CALDER (https://github.com/CSOgroup/CALDER)
rm(list=ls())

## --- Libraries --- ##
require(data.table)
require(doParallel)
require(GenomicRanges)
require(R.utils)
require(factoextra)
require(maptools)
require(CALDER)
require(patchwork)
require(dendextend)
require(scales)




## -- Project -- ##
prefix.path<-"/storage/data/FF/esalviat/SAMMY/Fibroblast/Paper_202301_January/Github/SAMMY-seq/"
prefix.name<-"Fibroblast_paper_202301_"



## -- Load utility functions -- ##
source(paste0(prefix.path,"Rscripts/SAMMY_Subcompartments_Calder_UtilityFunctions_github.R"))



### --- Input --- ###

## -- Main directories -- ##
#dir.plot<-paste0(prefix.path,"Plot/")
dir.track<- paste0(prefix.path,"Robjects/Track/")
dir.extra<- paste0(prefix.path,"Robjects/TrackExtra/")
dir.hic<-paste0(prefix.path,"Data/4DN/")

dir.out<-paste0(prefix.path,"Robjects/Subcompartment/")
if(!dir.exists(dir.out)) dir.create(dir.out)


## -- Main variables -- ##
binsize<-"50000"
dist.method<-"euclidean"

select.cor<-c("4f.untreated.C002-rep1_techX",        "4f.untreated.C004-rep1",       "4f.untreated.C004-rep2_techX")

## Palette
col.ab.8<-c(B.2.2="#ee4d5a",B.2.1="#f97b57",B.1.2="#f3ad6a",B.1.1="#ecda9a",A.2.2="#c4e6c3",A.2.1="#6dbc90",A.1.2="#36877a",A.1.1="#1d4f60")


### --- START script: Calculate sub-compartments for each chromosomes --- ###
CHR<-sapply(strsplit(dir(dir.track)[grep(paste0("_",binsize,"_"),dir(dir.track))],split="[_\\.]"),function(x) x[6])
ALL.file.track<-dir(dir.track)[grep(paste0("_",binsize,"_"),dir(dir.track))]
ALL.file.hic<-dir(dir.hic)[grep(paste0("_",binsize,"_"),dir(dir.hic))]
ALL.file.extra<-dir(dir.extra)[grep(paste0("_",binsize,"_"),dir(dir.extra))]



## --- Convert cooler output HiC matrix in Calder format -- ##
## Required bins and counts table (cooler)
## Run on cooler CLI (https://cooler.readthedocs.io/en/latest/cli.html#cooler-dump)
#cooler dump --table bins --header --out HFF_Finn_2019_50000_bins.txt.gz 4DNFIMDOXUT8.mcool::/resolutions/50000
#cooler dump --range chr18:0-80373285 --balanced --header --out HFF_Finn_2019_50000_chr18_counts.txt.gz 4DNFIMDOXUT8.mcool::/resolutions/50000

# Load bin tables (all chromosomes)
key<-c("HFF",paste0("_",binsize,"_"),"bins")
ind.key<-as.numeric(names(which(table(unlist(sapply(key,grep,dir(dir.hic))))==length(key))))
file.bin<-dir(dir.hic)[ind.key]
B<-data.table::fread(file=paste0(dir.hic,file.bin))


# HiC matrix in cooler format
for(chr in CHR){

	cat(">>>> ",chr,": ")
	name.save<-paste0(dir.hic,"HFF_CalderFormat_",binsize,"_",chr,".txt.gz")
	if(file.exists(name.save)){ cat("HiC matrix in Calder format already exists! Skip!\n"); next }

	# Read HiC matrix	
	key<-c("HFF_Finn",paste0("_",binsize,"_"),paste0("_",chr,"_"))
	ind.key<-as.numeric(names(which(table(unlist(sapply(key,grep,dir(dir.hic))))==length(key))))
	file.hic<-dir(dir.hic)[ind.key]

	cat("read HiC.. ")
	D<-data.table::fread(file=paste0(dir.hic,file.hic))
	D$bin1_start<-B$start[D$bin1_id+1]
	D$bin2_start<-B$start[D$bin2_id+1]

	take<-c("bin1_start","bin2_start","balanced")
	D<-D[,..take]

	cat("write new format.. ")
	# 0-based starts
	data.table::fwrite(D,file=name.save,col.names=FALSE, quote=FALSE, sep="\t")

	cat("\nend!\n\n")
}





## -- Calculate Subcompartments -- ##

for(i in 1:length(CHR)){
	
	chr<-CHR[i]
	cat("---> ",chr,": ")

	## check if file exist already ##
	name.chr<-paste0(dir.out,prefix.name,"subcalder_","ownblocks_",binsize,"_",chr,".RData")
	check0<-file.exists(name.chr)
	if(check0){ cat("\n\t",name.chr,"\n\tFile already exits.. next!\n\n"); next }


	### ---- Read SAMMY tracks --- ###
	# It reads tracks (already trasformed in Robjects) and trasform them in a distance matrix
	
	## key tracks
	key.track<-c(paste0("_",binsize),paste0("_",chr,"\\."))
	ind.track<-table(unlist(lapply(key.track,grep,ALL.file.track)))
	ind.track<-as.numeric(names(ind.track))[ind.track==length(key.track)]
	file.track<-ALL.file.track[ind.track]


	## key tracks extra (histon marks)
	key.extra<-c(paste0("_",binsize),paste0("_",chr,"\\."))
	ind.extra<-table(unlist(lapply(key.extra,grep,ALL.file.extra)))
	ind.extra<-as.numeric(names(ind.extra))[ind.extra==length(key.extra)]
	file.extra<-ALL.file.extra[ind.extra]


	## key hic
	key.hic<-c("_CalderFormat",paste0("_",binsize),paste0("_",chr,"\\."))
	ind.hic<-table(unlist(lapply(key.hic,grep,ALL.file.hic)))
	ind.hic<-as.numeric(names(ind.hic))[ind.hic==length(key.hic)]
	file.hic<-ALL.file.hic[ind.hic]


	# Return both distance matrix and genomic coordinates of bins
	cat("\n\t(1) Read tracks & hic:\n")
	
	# Load HiC data in calder format (already converted)
	A.hic<-read.hic.calder(dir.file=dir.hic,name.file=file.hic,bin.size=binsize,trans.01=FALSE,trans.log2=FALSE)

	AA<-lapply(select.cor,function(x) read.SAMMY.calder(dir.file=dir.track,name.file=file.track,bin.size=binsize,select=x,metric=dist.method,frac.select=NULL,bin.select=rownames(A.hic) ) )
	names(AA)<-select.cor

	AA<-append(AA,list(HiC=list(A=A.hic)))
	select.all<-append(select.cor,"HiC")


	### --- 1.1. Crossproduct correlation --- ###
	cat("\n\t(2) Cross-Correlation:\n")
	# Inverse hyperbolic tangent transformation (amplify signal >0.3); avoid Inf with (1+1E-7)
	CC<-lapply(select.all,function(x){ 
		cat("\t-- ",x,": ")
		compute.fastcor.calder(AA[[x]]$A,cor.cor=TRUE,trans.atanh=TRUE)
	})
	names(CC)<-select.all


	### --- 1.2. Identify domains/blocks boundaries --- ##
	cat("\n\t(3) Block boundaries.. ")

	BB<-lapply(select.all,function(x){
		cat("\t",x,"..\n")
		get.blocks.calder(A=CC[[x]],binsize,chr=chr)
	})
	names(BB)<-select.all


	### --- 3/4. Correlation of Trend matrix by block --- ###
	# Summarized (mean) contact/distance matrix  by block
	# Built the trend matrix (lags=4)
	# Correlation of trend matrix (atanh and not.scaled)

	# same blocks from HiC
	cat("\n\t(4) Trend-block matrices:\n")
	
	TT<-lapply(select.all,function(x){
		cat("\t-- ",x,": ")
		cor.trend.blocks(AA[[x]]$A,blocks=BB[[x]],lag=4,trans.atanh=TRUE,scale=TRUE,metric="mean")
	})
	names(TT)<-select.all


	### --- Get subcompartments --- ###
	# i.Principal component analysis on trend matrix
	# ii. k-means clustering
	# iii. sort with projected PC1/PC2

	#dir.out.track<-"/storage/data/FF/esalviat/SAMMY/Fibroblast/Compartments/Robjects/Track/"
	#files.track<-dir(dir.out.track)[grep(paste0("_",binsize,"_"),dir(dir.out.track))]
	TR.marks<-read.marks(dir.file=dir.extra,name.file=file.extra,bin.size=binsize,select="E055",bin.select=rownames(AA[[1]]$A))
	ind.active<-which(colnames(TR.marks) %in% c("Id","H3K36me3"))


	# $Bin: Bin-level subcompartments and reference block (i.e., domain)
	# $Block: Block-level subcompartments with pc1, pc2 and projected pca for sorting
	# $Dendro: iterative knn-results in dendrogram format 
	cat("\n\t(5) Subcompartments.. ")
	AB<-lapply(select.all,function(x){
		cat(x,".. ")
		get.subcompartment.calder(T=TT[[x]],blocks=BB[[x]],chr=chr,active="H3K36me3",tracks=TR.marks,n.comp=10,const.comp=5)
	})
	names(AB)<-select.all
	## Add some checks for the size and the match of columns (all(A.6f.JQ1$Bin==A.6f.DMSO$Bin))


	## check plots
	#plot.check<-plot.pca12(AB.obj=AB,select=c(select.cor,"HiC"),colour=col.ab.8 )
	#ggsave(plot.check,file=paste0(dir.plot,"PcaSub/",prefix.name,"pcasubcompartments_",chr,"_",binsize,".pdf"),height=10,width=0.6*(length(select.cor)+1))

	cat("\n\t(6) Save results..\n")
	Subcomp<-append(AB,AA[[1]]$Bin)
	names(Subcomp)[length(Subcomp)]<-"Bin"
	save(Subcomp,file=name.chr)

	cat("\n\n\n")
}

