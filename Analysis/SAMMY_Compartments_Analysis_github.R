### AB compartment analysis ###

rm(list=ls())
## -- libraries -- ##
require(data.table)
require(GenomicRanges)
#require(rtracklayer) # read tracks

require(reshape2)

require(ggplot2) # plot
require(patchwork) # merge multiple plots
require(RColorBrewer) # color density

require(ggforce)


## -- Project -- ##
prefix.path<-"/storage/data/FF/esalviat/SAMMY/Fibroblast/Paper_202301_January/Github/SAMMY-seq/"
prefix.name<-"Fibroblast_paper_202301_"

## -- Load utility functions -- ##
source(paste0(prefix.path,"Rscripts/SAMMY_Compartments_utilityFunctions_github.R"))

## -- Inputs -- ##

# Output directories (i.e., directories that contains Robjects)
dir.out.trackextra<-paste0(prefix.path,"Robjects/TrackExtra/")
dir.out.cor<-paste0(prefix.path,"Robjects/Correlation/")
dir.out.pca<-paste0(prefix.path,"Robjects/Eigenvector/")
dir.out.bin<-paste0(prefix.path,"Robjects/Bin/")


dir.plot<- paste0(prefix.path,"Plot/")
if(!dir.exists(dir.plot) ) dir.create(dir.plot)

dir.out.result<-paste0(prefix.path,"Results/")
if(!dir.exists(dir.out.result) ) dir.create(dir.out.result)


## -- Main variables -- ##
binsize<-"250000"
CHR<-"chr18"

scale.egv<-TRUE   # Divide by the abs(max) (visual purpose)
ref.hic<-"HiC"    # Reference data, to compare with

### Code of ChromHMM states (Roadmap)
code<-c("E055","E056")


select.cor.untreated<-c(
	"3f.untreated.C002-rep1",     "3f.untreated.C004-rep1", 
	"4f10K.untreated.C002-rep1",     "4f10K.untreated.C004-rep1",     "4f10K.untreated.C004-rep2",
	"4f.untreated.C002-rep1_techX",    "4f.untreated.C004-rep1",        "4f.untreated.C004-rep2_techX")

name.new.hic.untreated<-c(
	"3f.untreated.C002-rep1"="3f.C002-rep1",
	"3f.untreated.C004-rep1"="3f.C004-rep1",
	"4f10K.untreated.C002-rep1"="4f10K.C002.1",     
	"4f10K.untreated.C004-rep1"="4f10K.C004.1",
	"4f10K.untreated.C004-rep2"="4f10K.C004.2",
	"4f.untreated.C002-rep1_techX"="4f.C002.1tX", 
	"4f.untreated.C004-rep1"="4f.C004.1",
	"4f.untreated.C004-rep2_techX"="4f.C004.2tX",
	"HiC"="HiC")


## -- Color Palette -- ##
col.ab<-c("A-A"="#00A19D","B-B"="#FFB344","A->B"="#FFE194","B->A"="#B8DFD8")
col.sub<-c("A.1.1"="#FF0000","A.1.2"="#FF4848","A.2.1"="#FF9191","A.2.2"="#FFDADA",
		"B.1.1"="#DADAFF","B.1.2"="#9191FF","B.2.1"="#4848FF","B.2.2"="#0000FF","gap"="gray")

# The order of col.state control the order of the bars in the chromatin states plot (top-to-bottom)
col.state=c("16_gap"="black","9_Het"="#8A91D0","14_ReprPCWk"="#C0C0C0","13_ReprPC"="#808080", 
	"8_ZNF/Rpts"="#66CDAA","15_Quies"="#FFFFFF",
	"5_TxWk"="#006400", "4_Tx"="#008000","2_TssAFlnk"="#FF4500","1_TssA"="#FF0000",
	"7_Enh"="#FFFF00","6_EnhG"="#C2E105","12_EnhBiv"="#BDB76B",
	"11_BivFlnk"="#E9967A","10_TssBiv"="#CD5C5C", "3_TxFlnk"="#32CD32")

col.replica<-c(HiC="#A7D2CB",'C002-rep1'="#004D80",'C004-rep1'="#006C65",'C004-rep2'="#7FC87F","C001-rep1"="#970E53")

## -- Layout -- ##
# chromosome size to adjust the width of eigevencotrs plots
chr.size<-read.table(file=paste0(prefix.path,"Data/4DN/HFF_Finn_2019_4DNFIMDOXUT8_chromosomes.txt") )
colnames(chr.size)<-c("chr","size")




####### ------ Start script ------ #######


### --- Prepare data --- ########
## Load eigenvectors
files.pca<-dir(dir.out.pca)[grep(paste0("_",binsize,"_"),dir(dir.out.pca))]
PCA.all<-lapply(files.pca,function(x,dir){ get(load( paste0(dir,x)  )) },dir=dir.out.pca)
names(PCA.all)<-sapply(strsplit(files.pca,split="[_\\.]"),function(x) x[[6]])

## Change the sign (based on H3K36me3: active marker)
files.track<-dir(dir.out.trackextra)[grep(paste0("_",binsize,"_"),dir(dir.out.trackextra))]

Track.active<-lapply(files.track,function(x,dir){ 
	y<-get(load( paste0(dir,x)))
	return(y[["E055"]][["H3K36me3"]])
},dir=dir.out.trackextra)
names(Track.active)<-sapply(strsplit(files.track,split="[_\\.]"),function(x) x[[6]])
	
PCA.all<-lapply(names(PCA.all),function(x,DD,mm){    
	change.sign.activemarker(E.list=DD[[x]],marker=mm[[x]])
},DD=PCA.all,mm=Track.active)
names(PCA.all)<-sapply(strsplit(files.pca,split="[_\\.]"),function(x) x[[6]])


## Original Eigenvectors Distributions
PCA.all.notscaled<-PCA.all

if(scale.egv){

	# Divide each eigenvector by the maximum of the absolute values
	PCA.all<-lapply(names(PCA.all),function(x,DD){
		egv.chr<-DD[[x]]
		egv.chr.norm<-lapply(egv.chr,function(y) y/max(abs(y)) )
		return(egv.chr.norm)	
	},DD=PCA.all)
	names(PCA.all)<-sapply(strsplit(files.pca,split="[_\\.]"),function(x) x[[6]])

}



## Load bin genomic coordinate with ChromHMM info
files.bin<-lapply(code,function(x) dir(dir.out.bin)[grep(paste0("bininfostate15_",x,"_",binsize,"_"),dir(dir.out.bin))] )
names(files.bin)<-code

Bin.all.15<-lapply(code,function(e){
	bin.all.e<-lapply(files.bin[[e]],function(x,dir){ get(load( paste0(dir,x)  )) },dir=dir.out.bin)
	names(bin.all.e)<-sapply(strsplit(files.bin[[e]],split="[_\\.]"),function(x) x[[7]])
	return(bin.all.e)
})
names(Bin.all.15)<-code


HMM<-lapply(code,function(e){
	hmm.e<-lapply(CHR,function(x) get.datatable.binstates.chr(B=Bin.all.15[[e]][[x]],chr=x))
	names(hmm.e)<-CHR
	return(hmm.e)
})
names(HMM)<-code

## Chromosomes analysed
CHR<-intersect(names(PCA.all),names(Bin.all.15[[1]]))


## Get A and B compartments, using HiC reference 
# $res: all the information for each bin
# $n: number of A-A, B-B, A->B, B->A compartments (HiC->other)
AB.hic<-lapply(CHR,function(x,DD) get.datatable.eigenvector.chr(E=DD[[x]],chr=x,reference=ref.hic),DD=PCA.all)
names(AB.hic)<-CHR



##### ----- Figure 1c: Correlation matrices plots ----- ####
# chromosome: chr18
# protocols: HiC, 4f (C002-rep1)
# eigenvectors: 4 colours
# layout: vertical
# ChromHMM: 15 states

fig.chr<-"chr18"
fig.sel<-c("HiC","4f.untreated.C002-rep1_techX")
fig.state<-"state15"
fig.code<-"E055"
fig.col.state15<-col.state

for(i in 1:length(CHR)){
	
	fig.chr<-CHR[i]
	cat("---",fig.chr,"---")

	## Read only necessary objects
	files.cor<-dir(dir.out.cor)[grep(paste0("_",binsize,"_"),dir(dir.out.cor))]
	files.cor<-files.cor[grep(paste0("_",fig.chr,"."),files.cor)]

	M.cor<-lapply(files.cor,function(x,dir){ get(load( paste0(dir,x)  )) },dir=dir.out.cor)
	names(M.cor)<-sapply(strsplit(files.cor,split="[_\\.]"),function(x) x[[6]])
	M.cor<-M.cor[[fig.chr]][fig.sel]


	## Correlation matrix (data.table)
	dt.M.cor<-lapply(1:length(M.cor),function(i){ 
		x<-M.cor[[i]]

		# Row and columns contains Ids
		tab<-data.table::data.table(reshape2::melt(x))

		colnames(tab)<-c("row","col","corr")
		tab$fraction<-names(M.cor)[i]
		return(tab)
	})
	dt.M.cor<-rbindlist(dt.M.cor)

	## Eigenvectors
	dt.pca.chr<-AB.hic[[fig.chr]]$res[,c("Id","value","fraction","type"),with=FALSE]
	colnames(dt.pca.chr)<-c("col","eigen","fraction","type")
	dt.pca.chr<-dt.pca.chr[dt.pca.chr$fraction %in% fig.sel,]


	## Chromatin states (data.table)
	# Nb. HMM contains only selected bins
	hmm.chr<-Bin.all.15[[fig.code]][[fig.chr]]
	dt.hmm<-data.table(Id=hmm.chr$Id,hmm.chr$ChromHMM)
	dt.hmm<-reshape2::melt(dt.hmm,id.vars="Id",value.name="bp",variable.name="state")
	dt.hmm$state<-factor(paste0(dt.hmm$state),levels=(names(fig.col.state15)))


	## Layout parameters
	br.ax<-seq(0,max(c(dt.M.cor$col,dt.pca.chr$col) ),by=20)
	lab.ax<-br.ax*as.numeric(binsize)/1000000
    col.matrix<-colorRampPalette(brewer.pal(9,"RdBu")[c(-1,-9)])(201)


    ## Object with all plots in ggformat

    # Correlation matrices
	PLOT.cor<-vector("list",length(fig.sel))
	names(PLOT.cor)<-fig.sel

	for(fr in fig.sel){

		PLOT.cor[[fr]]<-ggplot(data=dt.M.cor[dt.M.cor$fraction==fr, ],aes(x=row,y=col,fill=corr ))+
			geom_tile(na.rm=TRUE)+
			theme_classic()+
			coord_equal()+
			facet_grid(~fraction)+
			scale_x_continuous(expand=c(0,5,0,5),breaks=br.ax,labels=lab.ax)+
			scale_y_continuous(position="right",trans="reverse",expand=c(0,5,0,5),breaks=br.ax,labels=lab.ax)+
			scale_fill_gradientn(colours=col.matrix, na.value = "white",limits = c(-1.01, 1.01))+

			labs(fill="Correlation\n",x="",y="Coordinates (Mb)")+
			theme(panel.spacing = unit(1.2, "lines"),
				axis.text.x=element_text(size=6),
				axis.text.y=element_text(size=6,angle=90),
				axis.title=element_text(size=10), 
				legend.key.width= unit(0.8, "lines"))
			
	}

    # Eigenvectors
	PLOT.eig<-vector("list",length(fig.sel))
	names(PLOT.eig)<-fig.sel

	for(fr in fig.sel){

		PLOT.eig[[fr]]<-ggplot(data=dt.pca.chr[dt.pca.chr$fraction==fr,],aes(x=col,y=eigen))+
			geom_col(aes(fill=type),colour=NA,width=1)+
			theme_classic()+
			facet_grid(~fraction)+
			scale_x_continuous(expand=c(0,5,0,5),breaks=br.ax,labels=lab.ax)+
			scale_fill_manual(values=col.ab)+
			scale_colour_manual(values=col.ab)+
			labs(fill="Correlation\n",x="",y="eigen")+
			theme(panel.spacing = unit(1.2, "lines"),axis.text=element_text(size=6),axis.title=element_text(size=10),legend.position="none")

	}

	## Plot states for each bin
	plot.state<-ggplot(data=dt.hmm)+
		geom_bar(aes(x=Id, y=bp,fill=state), stat="identity",size=0.1,colour=NA,width=1)+
		theme_classic()+
		scale_fill_manual(name="",values=alpha(fig.col.state15,0.6))+
		scale_colour_manual(name="",values=alpha(fig.col.state15,0.6))+
		scale_x_continuous(expand=c(0,5,0,5),breaks=br.ax,labels=lab.ax)+
		theme(legend.position="none")+
		#theme(legend.text = element_text(size=6.5),
		#		legend.key.height= unit(0.3, 'cm'),
		#        legend.key.width= unit(0.2, 'cm'),
		#        legend.position=c(1, 1),
		#        legend.justification=c(1, 0))+ # (, alto/basso=0,1)
		#labs(subtitle=paste0("Chromosome ",gsub("chr","",fig.chr)," ( ", format(as.numeric(binsize),big.mark=",",scientific = FALSE), " bp )"))+
		xlab("")+ylab("")+
		guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))

	cat("\n")
	ggsave(
		PLOT.cor[[fig.sel[1]]]+PLOT.eig[[fig.sel[1]]]+
		plot.state+
		PLOT.eig[[fig.sel[2]]]+PLOT.cor[[fig.sel[2]]]+
		plot_layout(ncol =1, heights=c(7,0.5,0.5,0.5,7) ),
		file=paste0(dir.plot,"Figure1c_correlationmatrices_",binsize,"_",fig.code,"_",fig.chr,".pdf"),width=8,height=17)


}


##### ----- Suppl.Figure 4c: Eigenvectors and ChromHMM ----- ####
# chromosome: genome-wide
# protocols: HiC, 4f (all replicas)
fig.chr<-"chr18"
fig.sel<-c("HiC","4f.untreated.C002-rep1_techX","4f.untreated.C004-rep1","4f.untreated.C004-rep2_techX")
fig.state<-"state15"
fig.code<-"E055"
fig.ref<-c("HiC")
fig.col.state15<-col.state


## Eigenvectors (data.table)
dt.pca.chr<-AB.hic[[fig.chr]]$res[,c("Id","value","fraction","type"),with=FALSE]
colnames(dt.pca.chr)<-c("col","eigen","fraction","type")
dt.pca.chr<-dt.pca.chr[dt.pca.chr$fraction %in% fig.sel,]


## Chromatin states (data.table)
# Nb. HMM contains only selected bins
hmm.chr<-Bin.all.15[[fig.code]][[fig.chr]]
dt.hmm<-data.table(Id=hmm.chr$Id,hmm.chr$ChromHMM)
dt.hmm<-reshape2::melt(dt.hmm,id.vars="Id",value.name="bp",variable.name="state")
dt.hmm$state<-factor(paste0(dt.hmm$state),levels=(names(fig.col.state15)))

## Cor Plots
br.ax<-seq(0,max(c(dt.pca.chr$col) ),by=20)
lab.ax<-br.ax*as.numeric(binsize)/1000000

PLOT.eig<-vector("list",length(fig.sel))
names(PLOT.eig)<-fig.sel

for(fr in fig.sel){

	PLOT.eig[[fr]]<-ggplot(data=dt.pca.chr[dt.pca.chr$fraction==fr,],aes(x=col,y=eigen))+
		geom_col(aes(fill=type),colour=NA,width=1)+
		theme_classic()+
		facet_grid(fraction~.)+
		scale_x_continuous(expand=c(0,5,0,5),breaks=br.ax,labels=lab.ax)+
		scale_y_continuous(expand=c(0,0.4,0,0.4))+
		scale_fill_manual(values=col.ab)+
		scale_colour_manual(values=col.ab)+
		labs(fill="Correlation\n",x="",y="eigen")+
		theme(panel.spacing = unit(1.2, "lines"),axis.text=element_text(size=6),axis.title=element_text(size=10),
			#legend.key.width= unit(0.8, "lines")
			legend.position="none")

}


## Plot states for each bin
plot.state<-ggplot(data=dt.hmm)+
	geom_bar(aes(x=Id, y=bp,fill=state), stat="identity",size=0.1,colour=NA,width=1)+
	theme_classic()+
	scale_fill_manual(name="",values=alpha(fig.col.state15,0.6))+
	scale_colour_manual(name="",values=alpha(fig.col.state15,0.6))+
	scale_x_continuous(expand=c(0,5,0,5),breaks=br.ax,labels=lab.ax)+
	theme(legend.position="none")+
	theme(legend.text = element_text(size=6.5),
			legend.key.height= unit(0.3, 'cm'),
	        legend.key.width= unit(0.2, 'cm'),
	        legend.position=c(1, 1),
	        legend.justification=c(1, 0))+ # (, alto/basso=0,1)
	labs(subtitle=paste0("Chromosome ",gsub("chr","",fig.chr)," ( ", format(as.numeric(binsize),big.mark=",",scientific = FALSE), " bp )"))+
	xlab("")+ylab("")+
	guides(
		colour = guide_legend(nrow = 1),
		fill = guide_legend(nrow = 1))


cat("\n")
ggsave(
	PLOT.eig[[ fig.sel[1] ]]+
	plot.state+
	PLOT.eig[[ fig.sel[2] ]]+PLOT.eig[[ fig.sel[3] ]]+PLOT.eig[[ fig.sel[4] ]]+
	plot_layout(ncol =1, heights=c(1,1.3,1,1,1) ),
	file=paste0(dir.plot,"SupplFigure4_eigenvectors_",binsize,"_",fig.code,"_",fig.chr,".pdf"),width=15,height=8)




