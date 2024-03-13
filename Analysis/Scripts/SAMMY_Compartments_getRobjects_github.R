###  SAMMY-seq analyses paper - January 2023 ###
# Create Robjects: 

### Experimental setup:
# 4f: many cells, DNAseI: S2S, S2L, (S3), (S4)
# 10kh: few cells, DNAseI: (S2), S3, S4

### Binsizes:
# 250 kb: A/B compartments, eigenvectors


rm(list=ls())
### --- libraries --- ###
require(data.table)
require(GenomicRanges)
#require(rtracklayer) # read tracks
require(pbapply) # parallelize apply functions
require(Matrix) 


## -- Project -- ##
prefix.path<-"/storage/data/FF/esalviat/SAMMY/Fibroblast/Paper_202301_January/Github/SAMMY-seq/"
prefix.name<-"Fibroblast_paper_202301_"

## -- Load utility functions -- ##
source(paste0(prefix.path,"Rscripts/SAMMY_Compartments_utilityFunctions_github.R"))

### --- Input --- ###

## -- Main directories -- ##
dir.hic<-paste0(prefix.path,"Data/4DN/")
dir.hmm<-paste0(prefix.path,"Data/ChromHMM/")



## HiC data (contact matrix Human Foreskin Fibroblast, 2019 - 4DNFIMDOXUT8)
dir.out.hic<-paste0(prefix.path,"Robjects/HiC/") 
## Bin genomic region with assciated chromHMM states (E055_15_coreMarks_hg38lift_dense.bed)
dir.out.bin<-paste0(prefix.path,"Robjects/Bin/")
## Bin RPKM for each protocols and fraction 
dir.out.track<-paste0(prefix.path,"Robjects/Track/")

dir.out.cor<-paste0(prefix.path,"Robjects/Correlation/")
dir.out.pca<-paste0(prefix.path,"Robjects/Eigenvector/")

## -- Main variables -- ##
binsize<-c("250000")
CHR<-c("chr18")
select.pca<-"PC1"

select.cor<-c(
	## untreated
	"3f.untreated.C002-rep1",     "3f.untreated.C004-rep1", 
	"4f10K.untreated.C002-rep1",     "4f10K.untreated.C004-rep1",     "4f10K.untreated.C004-rep2",
	"4f.untreated.C002-rep1_techX","4f.untreated.C004-rep1","4f.untreated.C004-rep2_techX"
)


## ChromHMM state only for chr18
# To repeat the analysis on all chromosomes download the complete files E055_15_coreMarks_hg38lift_dense.bed and E056_15_coreMarks_hg38lift_dense.bed
# Available at: https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/
name.hmm<-list(
	state15_E055="E055_15_coreMarks_hg38lift_dense_small.bed",
	state15_E056="E056_15_coreMarks_hg38lift_dense_small.bed")


####### ------ Start script ------ #######


### --- Read 2D Object (i.e., HiC) --- ###
# Required: a matrix of bin and a matrix of counts (cooler dump output, balanced matrix)

## Run on cooler CLI (https://cooler.readthedocs.io/en/latest/cli.html#cooler-dump)
#cooler dump --table bins --header --out HFF_Finn_2019_250000_bins.txt.gz 4DNFIMDOXUT8.mcool::/resolutions/250000
#cooler dump --range chr18:0-80373285 --balanced --header --out HFF_Finn_2019_250000_chr18_counts.txt.gz 4DNFIMDOXUT8.mcool::/resolutions/250000

ind.hic.bin<-sapply(binsize,function(bs)  Reduce("intersect",lapply(c("bin",paste0("_",bs)),function(x) grep(x,dir(dir.hic))))  )
file.bin<- dir(dir.hic)[ind.hic.bin]
names(file.bin)<-binsize


for(i in 1:length(CHR)){
	
	chr<-CHR[i]
	cat("--- ",chr," ---\n ")
	
	for(j in 1:length(binsize)){

		bin.j<-binsize[j]
		cat("-> ",bin.j,": ")

		## check if file exist already ##
		name.chr<-paste0(dir.out.hic,prefix.name,"hic_",bin.j,"_",chr,".RData")
		check0<-file.exists(name.chr)
		if(check0){ cat("\n\t",name.chr,"\n\tFile already exits.. next!\n\n"); next }


		file.count.chr<- dir(dir.hic)[Reduce("intersect",lapply(c("count",paste0("_",chr,"_"),paste0("_",bin.j)),function(x) grep(x,dir(dir.hic))))]
		if(length(file.count.chr)==0){ cat("\n\t","\n\tCounts file does not exist.. next!\n\n"); next }
		if(length(file.count.chr)>1){ cat("\n\t","\n\tAmbiguous name for counts file.. next!\n\n"); next }

		cat("read.. ")
		HiC<-read.hic(directory=dir.hic, file.count=file.count.chr, file.bin=file.bin[bin.j], chr=chr)

		cat("save.. ")
		save(HiC,file=name.chr)

		cat("\n")
	}
	cat("\n\n")

}




### --- Compute correlation matrices (and save bin in GenomicRanges coordinates)--- ####
if( !dir.exists(dir.out.cor) ) dir.create(dir.out.cor)

for(i in 1:length(CHR)){

	chr<-CHR[i]
	cat("--- ",chr," ---\n ")

	for(j in 1:length(binsize)){

		bin.j<-binsize[j]
		cat("-> ",bin.j,": ")

		## check if file exist already ##
		name.chr<-paste0(dir.out.cor,prefix.name,"corr_",bin.j,"_",chr,".RData")

		check0<-file.exists(name.chr)
		if(check0){ cat("\n\t",name.chr,"\n\tFile already exits.. next!\n\n"); next }

		cat("read.. ")
		## -- Inizialize list -- ##
		M.cor<-vector("list",length(select.cor)+1)
		names(M.cor)<-c((select.cor),"HiC")
		
		## -- Load Tracks -- ##
		file.chr.track<-dir(dir.out.track)[Reduce("intersect",lapply(c(paste0(chr,"\\.RData"),paste0("_",bin.j) ),function(x) grep(x,dir(dir.out.track))))]
		
		TR.split.chr<-get(load(paste0(dir.out.track,file.chr.track)))
		## -- Load HiC -- ##
		file.chr.hic<-dir(dir.out.hic)[Reduce("intersect",lapply(c(paste0(chr,"\\.RData"),paste0("_",bin.j)),function(x) grep(x,dir(dir.out.hic))))]
		HiC.chr<-get(load(paste0(dir.out.hic,file.chr.hic)))


		## check original sizes of the matrices
		check1<-all(HiC.chr$original.dim==sapply(unlist(TR.split.chr),length))
		if(!check1){ cat("Error: in",chr,"tracks and HiC matrix are not compatible"); next  }

		## Remove the same Bins both in tracks and HiC (Telomere or Centromere and other bins with problems, low mappability)
		ind.rem.1D<- lapply(TR.split.chr[select.cor],function(x){ data<-sapply(x,function(y) y$score); which(apply(data,1,sum)==0 ) })
		ind.rem.1D<-unique(unlist(ind.rem.1D))

		ind.rem.chr<-unique(c(HiC.chr$ind.remove,ind.rem.1D))
		

		cat("compute cor.. ")
		## 1D correlation
		M.cor[select.cor]<-lapply(TR.split.chr[select.cor],get.correlation.1D,ind.remove=ind.rem.chr)

		M.cor[["HiC"]]<-get.correlation.2D(Count=HiC.chr$Count,ind.remove=ind.rem.chr,orig.dim=HiC$original.dim)$Correlation
		
		## check size of matrices and bin ids
		check2<-all(table(unlist(lapply(M.cor,colnames)))==length(M.cor))
		if(!check2){ cat("Error: in",chr,"not all correlation matrices have the same bin ids"); next  }

		## -- Save bins coordinates in GenomicRanges format -- ##
		GR.bin<-granges(HiC$Bins)
		GR.bin$Id<-1:length(GR.bin)
		GR.bin<-GR.bin[-ind.rem.chr]

		check3<-all(colnames(M.cor[[1]])==GR.bin$Id)
		if(!check3){ cat("Error: in",chr,"M.cor and GR.bin Object not compatible"); next  }

		cat("bin.. ")
		name.chr.bin<-paste0(dir.out.bin,prefix.name,"bininfo_",bin.j,"_",chr,".RData")
		save(GR.bin,ind.rem.chr,file=name.chr.bin)


		cat("save.. ")
		save(M.cor,file=name.chr)
		cat("\n")
	}

	cat("\n\n")
}

### --- Compartments Using Standard PCA --- ###
if( !dir.exists(dir.out.pca) ) dir.create(dir.out.pca)

for(i in 1:length(CHR)){

	chr<-CHR[i]
	cat("--- ",chr," ---\n ")

	for(j in 1:length(binsize)){

		bin.j<-binsize[j]
		cat("-> ",bin.j,": ")

		## check if file exist already ##
		name.chr<-paste0(dir.out.pca,prefix.name,"pca_",bin.j,"_",chr,".RData")
		check0<-file.exists(name.chr)
		if(check0){ cat("\n\t",name.chr,"\n\tFile already exits.. next!\n\n"); next }

		cat("read.. ")
		## Load correlation matrices
		file.chr.cor<-dir(dir.out.cor)[Reduce("intersect",lapply(c(paste0(chr,"\\.RData"),paste0("_",bin.j) ),function(x) grep(x,dir(dir.out.cor))))]
		M.cor<-get(load(paste0(dir.out.cor,file.chr.cor)))

		cat("pca:\n")
		PCA<-lapply(1:length(M.cor),function(k){
			cat("\t-",names(M.cor)[k])
			
			x<-M.cor[[k]]
			pca<-prcomp(x,center = FALSE, scale =TRUE)#,na.action = na.omit)
			eig.vec1<-pca$rotation[,select.pca]

			return(eig.vec1)
			cat("\n")
		})
		cat("\n")
		names(PCA)<-names(M.cor)

		cat("save.. ")
		save(PCA,file=name.chr)
		cat("\n")
	}

	cat("\n\n")

}

#### --- Chromatin states in bins --- ###
## ChromHMM: download chromHMM states
GR.states<-lapply(name.hmm,function(x) rtracklayer::import(con=paste0(dir.hmm,x),format="bed"))

for(k in 1:length(GR.states)){

	cat("----- ",names(GR.states)[k]," -----\n")
	gr.state.k<-GR.states[[k]]

	for(i in 1:length(CHR)){

		chr<-CHR[i]
		cat("--- ",chr," ---\n ")

		for(j in 1:length(binsize)){

			bin.j<-binsize[j]
			cat("-> ",bin.j,": ")

			## check if file (with states) exist already ##
			name.chr<-paste0(dir.out.bin,prefix.name,"bininfo",names(GR.states)[k],"_",bin.j,"_",chr,".RData")
			check0<-file.exists(name.chr)
			if(check0){ cat("\n\t",name.chr,"\n\tFile already exits.. next!\n\n"); next }

			cat("read.. ")	

			## -- Load Bin matrices -- ##
			file.chr.bin<-dir(dir.out.bin)[Reduce("intersect",lapply(c("_bininfo_",paste0(chr,"\\.RData"),paste0("_",bin.j)),function(x) grep(x,dir(dir.out.bin))))]
			GR.bin<-get(load(paste0(dir.out.bin,file.chr.bin)))
			GR.bin<-get.states.by.bin(Bin=GR.bin,State=gr.state.k)
			

			cat("save.. ")
			save(GR.bin,file=name.chr)
			cat("\n")
		}

		cat("\n\n")

	}

}


















