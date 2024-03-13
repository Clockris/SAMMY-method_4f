## --- Utility functions --- ##


### ---- get Robjects --- ###
sparse.cor3 <- function(x){
    n <- nrow(x)

    cMeans <- colMeans(x)
    cSums <- colSums(x)

    # Calculate the population covariance matrix.
    # There's no need to divide by (n-1) as the std. dev is also calculated the same way.
    # The code is optimized to minize use of memory and expensive operations
    covmat <- tcrossprod(cMeans, (-2*cSums+n*cMeans))
    crossp <- as.matrix(Matrix::crossprod(x))
    covmat <- covmat+crossp

    sdvec <- sqrt(diag(covmat)) # standard deviations of columns
    covmat/crossprod(t(sdvec)) # correlation matrix
}


get.correlation.1D<-function(Tracks,ind.remove=NULL,variable="score",method="Pearson",scale=FALSE,center=FALSE){
	
	Data.cor<-lapply(Tracks,function(x) (x$score) )
	Data.cor<-data.table::data.table(Reduce("cbind",Data.cor))
	colnames(Data.cor)<-names(Tracks)
	bin<-1:nrow(Data.cor)

	## Remove some bins (e.g., centromere/telomere)
	if(length(ind.remove)>0){
		Data.cor<-Data.cor[-ind.remove,]
		bin<-bin[-ind.remove]
	}

	## Scale or center each track
	if( scale==TRUE | center ==TRUE ){
		Data.cor<-apply(Data.cor,2,function(x) scale(x,scale=scale,center=center) )
	}

	## Compute pairwise correlation between all pair of bins
	Mat.cor<-cor(t(Data.cor))
	colnames(Mat.cor)<-rownames(Mat.cor)<-bin

	return(Mat.cor)
}

read.hic<-function(directory, file.count, file.bin, chr){

	## Read counts
	count<-data.table::fread(file=paste0(directory,file.count))
	# Conver: 1-based
	count$bin1_id<-count$bin1_id+1
	count$bin2_id<-count$bin2_id+1

	## Read bins
	bin<-data.table::fread(file=paste0(directory,file.bin))
	bin$Id<-1:nrow(bin)

	# mantain only chromosome of interest
	ind.take<-which(bin$chr==chr)
	bin.chr<-bin[ind.take,]

	GR.bin.chr<-GRanges(
		seqnames=Rle(bin.chr$chrom),
		ranges=IRanges(start=as.numeric(bin.chr$start)+1,end=as.numeric(bin.chr$end) ),
		strand=Rle(rep("*",nrow(bin.chr))),
		bin.chr[,-c(1:3)])


	## Built matrix
	M<-Matrix::sparseMatrix(i=count$bin1_id,j=count$bin2_id,dims=c(max(ind.take),max(ind.take)),x=count$balanced,symmetric = TRUE)
	M.chr<-M[ind.take,ind.take]

	## Remove not balanced rows (i.e., weight == NA)
	ind.rem<-which(is.na(bin.chr$weight))
	#M.chr<-M.chr[-ind.rem,-ind.rem]
	#GR.bin.chr<-GR.bin.chr[-ind.rem,]

	Res<-list(Count=M.chr,Bins=GR.bin.chr,original.dim=length(ind.take),ind.remove=ind.rem,balanced=TRUE)
	return(Res)
}

get.correlation.2D<-function(Count,ind.remove=NULL,orig.dim){

	bin.ids<-1:orig.dim
	if(length(ind.remove)>0){
		Count<-Count[-ind.remove,-ind.remove]
		bin.ids<-bin.ids[-ind.remove]
	}

	n<-ncol(Count)
	df.count<-data.table::data.table(summary(Count))
	df.count$distance<-df.count$j-df.count$i

	# Could be improved....
	cat("Calculate expected matrix..\n")
	exp<-pbapply::pbsapply(0:(n-1),function(dd,data,n){
		## All elements at fixed distance dd
		vec<-data$x[data$distance==dd]

		# Fill the vector with 0s (they are removed in the sparse representation, but shuold be considered in the mean computation)
		n.zeros<- (n-dd)-length(vec)
		all.vec<-c(vec,rep(0,n.zeros))

		return(mean(all.vec,na.rm=TRUE))
	},data=df.count,n=n)
	df.exp<-data.table::data.table(distance=0:(n-1),exp=exp)
	df.count<-merge(x=df.count,y=df.exp,by="distance",all.x=TRUE)
	df.count$norm<- df.count$x/df.count$exp

	Mat.norm<-Matrix::sparseMatrix(i=df.count$i,j=df.count$j,dims=c(n,n),x=df.count$norm,symmetric = TRUE)	
	
	cat("Calculate correlation matrix..\n")
	Mat.cor<- sparse.cor3(Mat.norm)

	colnames(Mat.cor)<-rownames(Mat.cor)<-bin.ids
	colnames(Mat.norm)<-rownames(Mat.norm)<-bin.ids

	return(list(Correlation=Mat.cor,Normalized=Mat.norm))

}


get.states.by.bin<-function(Bin,State){

	all.state<-unique(mcols(State)[,"name"])

	name.gap<-paste0(length(all.state)+1,"_gap")
	all.state<-c(all.state,name.gap)
	all.state<-all.state[order(as.numeric(sapply(strsplit(all.state,split="_"),function(x) x[1])),decreasing=FALSE )]

	mcols(State)<-State$name
	names(mcols(State))<-"name"


	ind.state<-findOverlaps(Bin,State)

	info.state.overlap<-data.table(
		ind.bin=ind.state@from,
		ind.state=ind.state@to,
		state= match(paste0(State$name[ind.state@to]),all.state ),
		overlap=pintersect(granges(Bin[ind.state@from]),granges(State[ind.state@to]))@ranges@width )


	M.hmm<-sparseMatrix(
		i=info.state.overlap$ind.bin,
		j=info.state.overlap$state,
		x=info.state.overlap$overlap,
		dims=c(length(Bin),length(all.state)) )
	
	Bin$ChromHMM<-data.table(as.matrix(M.hmm))
	colnames(Bin$ChromHMM)<-all.state

	Bin$ChromHMM[[name.gap]]<-Bin@ranges@width-apply(Bin$ChromHMM,1,sum)

	# Sanity check: sum of all state must be lowe than bin size
	check0<-sum(Bin$ChromHMM[[name.gap]]<0)
	if(check0!=0) warning("The sum of bp in all the ChromHMM states is greater than the bin size!")

	return(Bin)
}

get.subcompartment.by.bin<-function(Bin,SC){

	all.type<-c("A.1.1","A.1.2","A.2.1","A.2.2","B.1.1","B.1.2","B.2.1","B.2.2","gap")
	ind.over<-findOverlaps(Bin,SC)

	SC$type<-paste0(SC$type)
	
	# Slow but it works
	Tab.subcomp<-pbapply::pblapply(1:length(Bin),function(i,all){
		#cat(i," ")
		bb<-Bin[i,]

		## Subset only states that overlap with the bin
		sub<-subsetByOverlaps(SC, bb)
		sub<-pintersect(sub,rep(bb,length(sub)))
		sub<-sort(c(sub,intersect(gaps(sub),bb)))

		sub$type[is.na(sub$type)]<-"gap"

		## check if there are some overlaps between subcompartments
		check<-all(table(findOverlaps(sub,sub)@from)==1)
		if(!check){ cat("\noverlapping subcompartments in bin ",i,"\n") }

		
		db<-data.table::data.table(bp=sub@ranges@width,type=sub$type)
		tab<-db[, sum(bp),by=type]
		
		res<-tab$V1
		names(res)<-tab$type
		res<-res[all]
		res[is.na(res)]<-0
		names(res)<-all
	
		return(res)

	},all=all.type)

	DT.subcomp<-data.table::data.table(Reduce("rbind",Tab.subcomp))

	mcols(Bin)<-cbind(mcols(Bin),DT.subcomp)
	return(Bin)


}


### ---- Plots ---- ###

change.sign.activemarker<-function(E.list,marker){

	# bin number
	marker$bin<-1:length(marker)
	bin<-as.numeric(names(E.list[[1]]))

	marker<-marker[match(bin,marker$bin),]

	E.list<-lapply(E.list,function(x,vec){
		# should be positively associated: A is high value of active marker
		return( x*sign(cor(x,vec)) )
	},vec=marker$score)

	return(E.list)
}


get.datatable.eigenvector.chr<-function(E,chr,others=NULL,reference){

	if(is.null(others)){
		select<-names(E)
	} else {
		select<-others[others %in% names(E)]
	}


	dt.ab.dif<-data.table::data.table(
		value=unlist(E[select]),
		fraction=rep(names(E[select]),times=sapply(E[select],length)),
		bin=unlist(lapply(sapply(E[select],length),seq_len)),
		value.ref=rep(E[[reference]],length(select)),
		chr=chr,
		Id=as.numeric(unlist(lapply(E[select],names)))
	)


	dt.ab.dif$sign<-sign(dt.ab.dif$value)
	dt.ab.dif$sign.ref<-sign(dt.ab.dif$value.ref)

	dt.ab.dif$type<-paste(dt.ab.dif$sign.ref,dt.ab.dif$sign,sep=";") # ref is HiC
	dt.ab.dif$type<-factor(paste0(dt.ab.dif$type),levels=c("-1;-1","1;-1","-1;1","1;1"),labels=c("B-B","A->B","B->A","A-A"))


	## Add absolute number of discordant/concordant
	dt.ab.num<-reshape2::melt( table(value=dt.ab.dif$sign,value.ref=dt.ab.dif$sign.ref,fraction=dt.ab.dif$fraction),value.name="number")
	dt.ab.num$type<-paste(dt.ab.num$value.ref,dt.ab.num$value,sep=";") # ref is HiC
	
	dt.ab.num$type<-factor(paste0(dt.ab.num$type),levels=c("-1;-1","1;-1","-1;1","1;1"),labels=c("B-B","A->B","B->A","A-A"))
	dt.ab.num$chr<-chr

	return(list(res=dt.ab.dif,n=dt.ab.num))

}

get.datatable.binstates.chr<-function(B,chr){

	# B$id refers to the original size of the matrices
	# Add bin number (after removal)
	B$bin<-1:length(B)
	B<-B[,-which(colnames(mcols(B))=="Id")]

	# Merge info on states and classification of states
	data.bin<-B$ChromHMM
	data.bin$bin<-B$bin
	
	dt.bin.states<-reshape2::melt(data.bin,id.vars="bin",value.name = "bp",variable.name="state")
	dt.bin.states$chr<-chr

	return(dt.bin.states)

}

get.datatable.subcompartment.chr<-function(B,chr){

	B$bin<-1:length(B)
	
	# Merge info on states and classification of states
	data.bin<-as.data.table(mcols(B))
	colnames(data.bin)<-colnames(mcols(B))
	dt.bin.subc<-reshape2::melt(data.bin,id.vars="bin",value.name = "bp",variable.name="type")
	dt.bin.subc$chr<-chr

	# add the sign to the B compartments
	dt.bin.subc$bp.sign<-dt.bin.subc$bp
	ind.b<-grep("B",dt.bin.subc$type)
	dt.bin.subc$bp.sign[ind.b]<- -1* dt.bin.subc$bp.sign[ind.b]


	return(dt.bin.subc)
}


get.datatable.subcompartments.chr<-function(B,chr){

	# B$id refers to the original size of the matrices
	# Add bin number (after removal)
	B$bin<-1:length(B)
	B<-B[,-which(colnames(mcols(B))=="Id")]

	# Merge info on states and classification of states
	data.bin<-as.data.table(mcols(B))
	colnames(data.bin)<-colnames(mcols(B))
	dt.bin.states<-reshape2::melt(data.bin,id.vars="bin",value.name = "bp",variable.name="state")
	dt.bin.states$chr<-chr

	return(dt.bin.states)

}


plot.eigenvector.state<-function(E,B,reference,col.ab,col.state,name.sorted,name.new=NULL,labs.size=5){


	## Mantain only the selected fractions + reference
	E$res<-E$res[which(paste0(E$res$fraction) %in% c(name.sorted,reference)),]


	## Plot same bin 
	select.bin<-unique(E$res$bin)
	select.bin<-sort(select.bin)
	B.sel<-B[which(B$bin %in% select.bin),]
	
	ind.ref<-which(paste0(E$res$fraction)==reference)
	chr<-unique(B$chr)

	## Order chromatin states (based on the order of col.state variable)
	B.sel$state<-factor(paste0(B.sel$state),levels=names(col.state))


	## Order fractions (based on name.sorted variable)
	if(!is.null(name.new)){
		E$res$fraction<-factor(paste0(E$res$fraction),levels=c(name.sorted,reference),labels=name.new[c(name.sorted,reference)] )
	} else { 
		E$res$fraction<-factor(paste0(E$res$fraction),levels=c(name.sorted,reference))
	}


	## Limits for y-axis
	const<-100
	ymin<-min(c(E$res$value,E$res$value.ref))
	ymax<-max(c(E$res$value,E$res$value.ref))
	ylim<- c(-1,1)* ceiling( max(abs(c(ymin,ymax))) * const )/const

	## Limits for x-axis
	xlim<-c(min(select.bin)-1,max(select.bin)+1) 

	## Barplot for the eigenvectors of the reference
	plot.ref<-ggplot(data=E$res[ind.ref,])+
		geom_col(aes(x=bin, y=ifelse(value<0, value, 0), fill=type,colour=type),size=0.1)+
		geom_col(aes(x=bin, y=ifelse(value>0, value, 0), fill=type,colour=type),size=0.1)+
		theme_classic()+
		scale_fill_manual(values=alpha(col.ab,0.5))+
		scale_colour_manual(values=alpha(col.ab,0.5))+
		scale_x_continuous(limits=xlim,expand=c(0,1,0,1))+
		scale_y_continuous(breaks=c(ylim[1],0,ylim[2]),limits=ylim)+
		facet_grid(fraction~.)+labs(title= paste0("Chromosome ",gsub("chr","",chr)), subtitle="First eigenvectors" )+
		xlab("")+ylab("")

	#ggsave(plot.ref,file=paste0(dir.plot,"Eigenvectors_States_ref_",chr,"_",binsize,".pdf"),width=15,height=2.5)


	## Barplot for the eigenvectors to compare with the reference 
	plot.AB<-ggplot(data=E$res[-ind.ref,])+
		geom_col(aes(x=bin, y=ifelse(value<0, value, 0), fill=type,colour=type),size=0.1)+
		geom_col(aes(x=bin, y=ifelse(value>0, value, 0), fill=type,colour=type),size=0.1)+
		theme_classic()+
		theme(strip.text.y=element_text(size=labs.size))+
		scale_fill_manual(values=alpha(col.ab,0.5))+
		scale_colour_manual(values=alpha(col.ab,0.5))+
		scale_x_continuous(limits=xlim,expand=c(0,1,0,1))+
		scale_y_continuous(breaks=c(ylim[1],0,ylim[2]),limits=ylim)+
		facet_grid(fraction~.)+#labs(title= paste0("Chromosome ",gsub("chr","",chr)), subtitle="First eigenvectors" )+
		xlab("")+ylab("")
	#ggsave(plot.AB,file=paste0(dir.plot,"Eigenvectors_States_others_",chr,"_",binsize,".pdf"),width=15,height=6)


	## Plot states for each bin
	plot.state<-ggplot(data=B.sel)+
		geom_bar(aes(x=bin, y=bp,fill=state,colour=state), stat="identity",size=0.1)+
		theme_classic()+
		scale_fill_manual(name="",values=alpha(col.state,0.6))+
		scale_colour_manual(name="",values=alpha(col.state,0.6))+
		scale_x_continuous(limits=xlim,expand=c(0,1,0,1))+
		theme(legend.text = element_text(size=6.5),
			legend.key.height= unit(0.3, 'cm'),
	        legend.key.width= unit(0.2, 'cm'),
	        legend.position=c(1, 1),
	        legend.justification=c(1, 0))+ # (, alto/basso=0,1)
		labs(subtitle="ChromHMM states for each bin")+xlab("")+ylab("")+
		guides(
			colour = guide_legend(nrow = 1),
			fill = guide_legend(nrow = 1))
	
	#ggsave(plot.state,file=paste0(dir.plot,"Eigenvectors_state_",chr,"_",binsize,".pdf"),width=15,height=2.5)
	#ggsave( (plot.ref / plot_spacer() / plot.AB ) / plot.state  + plot_layout(height=c(1.5,-0.5,5.5,2.5)),file=paste0(dir.plot,"Eigenvectors_Allinone_",chr,"_",binsize,".pdf"),width=18,height=7)

	return(list(ref=plot.ref,others=plot.AB,state=plot.state))
}


plot.eigenvector.subcompartment<-function(E,B,reference,col.ab,col.sub,binsize){

	binsize<-as.numeric(binsize)

	## Plot same bin 
	select.bin<-unique(E$res$bin)
	select.bin<-sort(select.bin)
	B.sel<-B[which(B$bin %in% select.bin),]
	
	ind.ref<-which(E$res$fraction==reference)
	chr<-unique(B$chr)

	## Limits for y-axis
	const<-100
	ymin<-min(c(E$res$value,E$res$value.ref))
	ymax<-max(c(E$res$value,E$res$value.ref))
	ylim<- c(-1,1)* ceiling( max(abs(c(ymin,ymax))) * const )/const

	## Limits for x-axis
	xlim<-c(min(select.bin)-1,max(select.bin)+1) 


	## Barplot for subcompartments (eight classes)
	# wrong size & gap
	B.sel<-B.sel[which(B.sel$type!="gap"),]
	ind.id.take<- which( !(paste0(B.sel$bin) %in% names(tab[which(tab>250000)]) ) )

	plot.sub.ab<-ggplot(data=B.sel[ind.id.take,])+
		geom_col(aes(x=bin, y=ifelse(bp.sign<0, bp.sign, 0), fill=type,colour=type),size=0.1)+
		geom_col(aes(x=bin, y=ifelse(bp.sign>0, bp.sign, 0), fill=type,colour=type),size=0.1)+
		theme_classic()+
		scale_fill_manual(values=col.sub)+
		scale_colour_manual(values=col.sub)+
		scale_x_continuous(limits=xlim,expand=c(0,1,0,1))+
		scale_y_continuous(breaks=c(-binsize,0,binsize),limits=c(-1,1)*binsize)+
		#facet_grid(fraction~.)+#labs(title= paste0("Chromosome ",gsub("chr","",chr)), subtitle="First eigenvectors" )+
		xlab("")+ylab("")+
		theme( legend.key.height= unit(0.4, 'cm'),legend.key.width= unit(0.7, 'cm') )

	#ggsave(plot.sub.ab,file=paste0(dir.plot,"Eigenvectors_Subcomp_",chr,"_",binsize,".pdf"),width=15,height=2.5)



	## Barplot for the eigenvectors of the reference
	plot.ref<-ggplot(data=E$res[ind.ref,])+
		geom_col(aes(x=bin, y=ifelse(value<0, value, 0), fill=type,colour=type),size=0.1)+
		geom_col(aes(x=bin, y=ifelse(value>0, value, 0), fill=type,colour=type),size=0.1)+
		theme_classic()+
		scale_fill_manual(values=alpha(col.ab,0.5))+
		scale_colour_manual(values=alpha(col.ab,0.5))+
		scale_x_continuous(limits=xlim,expand=c(0,1,0,1))+
		scale_y_continuous(breaks=c(ylim[1],0,ylim[2]),limits=ylim)+
		facet_grid(fraction~.)+labs(title= paste0("Chromosome ",gsub("chr","",chr)), subtitle="First eigenvectors" )+
		xlab("")+ylab("")

	#ggsave(plot.ref,file=paste0(dir.plot,"Eigenvectors_States_",chr,"_",binsize,".pdf"),width=15,height=2.5)


	## Barplot for the eigenvectors to compare with the reference 
	plot.AB<-ggplot(data=E$res[-ind.ref,])+
		geom_col(aes(x=bin, y=ifelse(value<0, value, 0), fill=type,colour=type),size=0.1)+
		geom_col(aes(x=bin, y=ifelse(value>0, value, 0), fill=type,colour=type),size=0.1)+
		theme_classic()+
		scale_fill_manual(values=alpha(col.ab,0.5))+
		scale_colour_manual(values=alpha(col.ab,0.5))+
		scale_x_continuous(limits=xlim,expand=c(0,1,0,1))+
		scale_y_continuous(breaks=c(ylim[1],0,ylim[2]),limits=ylim)+
		facet_grid(fraction~.)+#labs(title= paste0("Chromosome ",gsub("chr","",chr)), subtitle="First eigenvectors" )+
		xlab("")+ylab("")
	#ggsave(plot.AB,file=paste0(dir.plot,"Eigenvectors_States_",chr,"_",binsize,".pdf"),width=15,height=6)



	#ggsave( (plot.ref / plot_spacer() / plot.AB ) / plot.state  + plot_layout(height=c(1.5,-0.5,5.5,2.5)),file=paste0(dir.plot,"Eigenvectors_Allinone_",chr,"_",binsize,".pdf"),width=18,height=7)

	return(list(ref=plot.ref,others=plot.AB,subcompartment=plot.sub.ab))
}




plot.compartment.density<-function(E.all,others,reference,col.ab,chr=NULL,n.contour=15,add.point=FALSE,names.new=NULL){

	th.text<-0.4
	th.legend<-0.2


	if(length(chr)>0){		
		
		## Coordinates
		E.xy<-E.all[[chr]]$res
		E.xy<-E.xy[E.xy$fraction %in% c(others,reference),]
		
	
		## Text (number of compartments)
		E.n.form<-E.all[[chr]]$n
		E.n.form<-E.n.form[E.n.form$fraction %in% c(others,reference),]

		E.n.form$value<-E.n.form$value*th.text
		E.n.form$value.ref<-E.n.form$value.ref*th.text


	} else {

		E.xy<-rbindlist(lapply(E.all,function(x){ return(x$res) }))
		# select fractions
		E.xy<-E.xy[E.xy$fraction %in% c(others,reference),]
		
		## Test format
		E.n<-rbindlist(lapply(E.all,function(x) x$n))
		E.n<-E.n[E.n$fraction %in% c(others,reference),]

		# get value and value.ref
		tab<-unique(E.n[,c("value","value.ref","type")])
		E.n.form<-reshape2::melt(tapply(E.n$number,list(fraction=E.n$fraction,type=E.n$type),sum))
		colnames(E.n.form)<-c("fraction","type","number")
		E.n.form<-merge(E.n.form,tab,all.x=TRUE)
		E.n.form<-E.n.form[!is.na(E.n.form$number),]

		E.n.form$number<-paste0(format(E.n.form$number,big.mark=","))
		E.n.form$value<-E.n.form$value*th.text
		E.n.form$value.ref<-E.n.form$value.ref*th.text

		
	}

	if(!is.null(names.new)){
		E.xy$fraction<-factor(paste0(E.xy$fraction),levels=c(reference,others),labels=names.new[c(reference,others)])
		E.n$fraction<-factor(paste0(E.n$fraction),levels=c(reference,others),labels=names.new[c(reference,others)])
		E.n.form$fraction<-factor(paste0(E.n.form$fraction),levels=c(reference,others),labels=names.new[c(reference,others)])

		reference<-names.new[reference]

	} else {
		E.xy$fraction<-factor(paste0(E.xy$fraction),levels=c(reference,others))
		E.n$fraction<-factor(paste0(E.n$fraction),levels=c(reference,others))
		E.n.form$fraction<-factor(paste0(E.n.form$fraction),levels=c(reference,others))

	}

	# Limits and breaks
	lims<- c(-1,1)
	br.axes<- seq(from=lims[1],to=lims[2],by=0.25)
	
	## Palette
	get.palette.density<-colorRampPalette(brewer.pal(9, "YlGnBu"))
	col.density<-c("white",rev(get.palette.density(n.contour-1)))

	frac<-setdiff(levels(E.n.form$fraction),reference)
	PLOT<-vector("list",length(frac))
	names(PLOT)<-frac


	for(i in 1:length(frac)){
		
		# Remove reference
		ind.take<-which(E.xy$fraction == frac[i])
		ind.take.txt<-which(E.n.form$fraction == frac[i])
		

		PLOT[[i]]<-ggplot()+
			geom_density_2d_filled(data=E.xy[ind.take,],aes(x=value,y=value.ref,fill=..level..,alpha=..level..),bins=n.contour,colour="white",size=0.1)+
			geom_text(data=E.n.form[ind.take.txt,],aes(x=value,y=value.ref,label=number),size=2.8)+
			
			geom_label(data=NULL,aes(x=0.2,y=0.2),colour="white",label="A-A",fill=col.ab["A-A"],fontface = "bold",size=2)+
			geom_label(data=NULL,aes(x=-0.2,y=-0.2),colour="white",label="B-B",fill=col.ab["B-B"],fontface = "bold",size=2)+
			geom_label(data=NULL,aes(x=-0.2,y=0.2),colour="white",label="A->B",fill=col.ab["A->B"],fontface = "bold",size=2)+
			geom_label(data=NULL,aes(x=0.2,y=-0.2),colour="white",label="B->A ",fill=col.ab["B->A"],fontface = "bold",size=2)+

			scale_fill_manual(values =col.density )+
			scale_colour_manual(values=col.ab)+
			geom_hline(yintercept=0,size=0.1)+
			geom_vline(xintercept=0,size=0.1)+
			theme_classic()+coord_fixed()+
			scale_x_continuous(expand=c(0,0.02,0,0.02),limits=lims)+#,,limits=lims)+
			scale_y_continuous(expand=c(0,0.02,0,0.02),limits=lims)+#,limits=lims)+
			facet_grid(~fraction)+ylab(paste0(reference," first eigenvector"))+xlab(paste0(frac[i]," first eigenvector"))+
			theme(legend.position="none",panel.spacing = unit(1.2, "lines"))

		if(add.point) PLOT[[i]]<-PLOT[[i]]+geom_point(data=E.xy[ind.take,],aes(x=value,y=value.ref),shape=21,fill="darkgray",alpha=0.1,size=0.2)

		#ggsave(PLOT[[i]],file=paste0(dir.plot,"ABdifference_",frac[i],"_",chr,"_",binsize,".pdf"),width=3,height=3)
	}


	
	return(PLOT)

}


plot.compartment.summary.chr<-function(E.all,type.plot="Number",col.ab,reference,others,names.new=NULL,span.thick=250){


	E.n<-rbindlist(lapply(E.all,function(x){ return(x$n) }))
	#E.n<-E.n[-which(E.n$fraction==reference),]
	E.n<-E.n[which(E.n$fraction %in% others),]

	## Add column with the summary for all chromosomes
	E.total.n<-reshape2::melt(tapply(E.n$number,list(value=E.n$value,value.ref=E.n$value.ref,fraction=E.n$fraction,type=E.n$type),sum),value.name="number")
	E.total.n<-data.table(E.total.n[!is.na(E.total.n$number),])
	E.total.n$chr<-"Genome"
	E.total.n<-E.total.n[,colnames(E.n),with=FALSE]

	E.n<-rbind(E.n,E.total.n)

	tab.chr<-reshape2::melt(tapply(E.n$number,list(E.n$chr,paste0(E.n$fraction)),sum))
	colnames(tab.chr)<-c("chr","fraction","total")

	E.n.perc<-merge(E.n,tab.chr,by=c("chr","fraction"),all.x=TRUE)
	E.n.perc$perc<-E.n.perc$number/E.n.perc$total*100

	# Bars order
	E.n.perc$type<-factor(paste0(E.n.perc$type),levels=c("A-A","B->A","A->B","B-B"))
	E.n.perc$chr<-factor(paste0(E.n.perc$chr),levels=c("Genome",paste0("chr",1:22),"chrX","chrY"))

	plot<-NULL

	if(!is.null(names.new)){
		E.n.perc$fraction<-factor(paste0(E.n.perc$fraction),levels=others,labels=names.new[others])
		tab.chr$fraction<-factor(paste0(tab.chr$fraction),levels=others,labels=names.new[others])
		reference<-names.new[reference]

	} else {
		E.n.perc$fraction<-factor(paste0(E.n.perc$fraction),levels=c(reference,others))
		tab.chr$fractions<-factor(paste0(tab.chr$fraction),levels=c(reference,others))
		
	}


	if(type.plot=="Number"){

		n.max<-max(tapply(E.n$number,list(E.n$fraction,E.n$chr),sum),na.rm=TRUE)
		br.n<-seq(0,(n.max %/% span.thick)*span.thick,by=span.thick)

		plot<-ggplot(data=E.n.perc,aes(x=chr,y=number))+
			geom_bar(aes(fill=type),stat="identity",position = "stack",width=0.8)+
			theme_classic()+facet_grid(~fraction)+
			scale_fill_manual(values=col.ab)+
			theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=7))+
			scale_y_continuous(breaks=br.n,expand=c(0,0,0.1,0))+
			xlab("")+ylab("Number of bins")+labs(subtitle=paste0("Number of bins for each type of compartment, compared to the reference (",reference,")") )
	}

	if(type.plot=="Percentage"){

		plot<-ggplot(data=E.n.perc,aes(x=chr,y=perc))+
			geom_bar(aes(fill=type),stat="identity",position = "stack",width=0.8)+
			geom_text(data=tab.chr,aes(x=chr,y=+Inf,label=total),angle=90,size=2,hjust=1.2)+
			theme_classic()+facet_grid(~fraction)+
			scale_fill_manual(values=col.ab)+
			theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=7))+
			scale_y_continuous(breaks=(0:10)*10,expand=c(0,0,0.15,0))+
			xlab("")+ylab("Percentage of bins")+labs(subtitle=paste0("Percentage of bins for each type of compartment, compared to the reference (",reference,")") )

	} 

	if(length(plot)==0) { cat("Only 'Percentage' or 'Number' are valid 'type.plot' options") }

	#ggsave(plot,file=paste0(dir.plot,"ABdifference_resume_",binsize,".pdf"),width=13,height=3)

	return(plot)

}


## Inputs
#ab<-AB.hic
#bin.hmm<-Bin.all.15[["E055"]]
#sel.fractions<-fig.sel
#sel.states<-names(fig.col.state15)
#sel.chromosomes<-CHR

plot.compartment.change.perc<-function(ab,bin.hmm,sel.fractions,sel.states,sel.chromosomes,layout="vertical"){

	# Merge info on ab and states
	SS<-vector("list",length(sel.chromosomes))
	names(SS)<-sel.chromosomes

	
	for(i in 1:length(sel.chromosomes)){

		chr.i<-sel.chromosomes[i]

		ab.chr<-ab[[chr.i]]$res[,c("chr","fraction","bin","Id","type"),with=FALSE]
		ab.chr<-ab.chr[fraction %in% sel.fractions,]

		hmm.chr<-as.data.table( bin.hmm[[chr.i]]$ChromHMM[,sel.states,with=FALSE] )
		hmm.chr$Id<-bin.hmm[[chr.i]]$Id

		states.chr<-merge(ab.chr,hmm.chr,by="Id",all.x=TRUE)
		SS[[chr.i]]<-states.chr
	}
	SS.all<-rbindlist(SS)
	

	# Total bps occupy by each state
	tot.bp.state<-apply( SS.all[,names(fig.col.state15),with=FALSE] ,2,sum)/length(fig.sel)
	tot.bp.state<-data.table(state=names(tot.bp.state),tot=tot.bp.state)


	SS.all<-reshape2::melt(SS.all,id.vars=c("Id","chr","fraction","bin","type"),value.name="bp",variable.name="state")


	## -- Data barplots: % of states -- ##
	# Count number of bp per state and fraction and type
	tab.all.state<-data.table(reshape2::melt(tapply(SS.all$bp,list(SS.all$fraction,SS.all$type,SS.all$state),sum)))
	colnames(tab.all.state)<-c("fraction","type","state","bp")
	tab.all.state<-tab.all.state[!is.na(tab.all.state$bp),]

	# Get percentages
	tab.all.state<-merge(tab.all.state,tot.bp.state,by="state",all.x=TRUE)
	tab.all.state$perc<-tab.all.state$bp/tab.all.state$tot*100
	tab.all.state$ab<-sapply(strsplit(paste0(tab.all.state$type),split="[->]"),function(x) x[length(x)])



	## -- Data Heatmap: differences between HiC and fractions -- ##
	tab.diff.state<-data.table(reshape2::melt( tapply(tab.all.state$bp,list(tab.all.state$fraction,tab.all.state$ab,tab.all.state$state),sum ) ) )
	colnames(tab.diff.state)<-c("fraction","type","state","bp")
	tab.diff.state<-merge(tab.diff.state,tot.bp.state,by="state",all.x=TRUE)
	tab.diff.state$perc<-tab.diff.state$bp/tab.diff.state$tot*100

	## differences respect to "A" (it's complementary to B)
	tab.diff.state.A<-tab.diff.state[tab.diff.state$type=="A",]

	ind.hic<-which(tab.diff.state.A$fraction=="HiC")
	tab.diff.hic<-tab.diff.state.A[ind.hic,]
	tab.diff.state.A<-tab.diff.state.A[-ind.hic,]

	tab.diff.state.A<-merge(tab.diff.state.A,tab.diff.hic,by=c("state","type"))
	tab.diff.state.A$diff<- tab.diff.state.A$perc.x-tab.diff.state.A$perc.y


	## sort differences by mean of selected fractions

	#ind.ref<-which(tab.diff.state.A$fraction.x==sort.by)
	name.state.sorted<-names(sort(tapply(tab.diff.state.A$diff,tab.diff.state.A$state,mean),decreasing=TRUE))
	#name.state.sorted<- tab.diff.state.A$state[ind.ref][order(tab.diff.state.A$diff[ind.ref],decreasing=TRUE)]

	tab.diff.sub<-tab.diff.state.A[,c("state","fraction.x","diff")]
	col.max<-round(max(abs(tab.diff.sub$diff),na.rm=TRUE)+1)

	tab.diff.sub$state<-droplevels(factor(paste0(tab.diff.sub$state),levels=(name.state.sorted)))
	tab.diff.sub$fraction.x<-factor(paste0(tab.diff.sub$fraction.x),levels=rev(setdiff(sel.fractions,"HiC")),labels=gsub("4f.untreated.","" ,rev(setdiff(sel.fractions,"HiC")) )  )


	# figures options
	tab.all.state$type<-factor(paste0(tab.all.state$type),levels=rev(c("A-A","B->A","A->B","B-B")))
	tab.all.state$fraction<-factor(paste0(tab.all.state$fraction),levels=sel.fractions,labels=gsub("4f.untreated.","" ,sel.fractions)  )
	tab.all.state$state<-droplevels(factor(paste0(tab.all.state$state),levels=name.state.sorted))


	#### --- Table to save --- ###
	tab.save<-split(tab.all.state[,c("state","type","bp","perc"),with=FALSE],tab.all.state$fraction)
	TAB<-lapply(names(tab.save),function(x){
		y<-tab.save[[x]]

		y.st<-split(y[,c("state","bp","perc"),with=FALSE],y$type)
		y.st<-y.st[ which(sapply(y.st,nrow)>0)]
		for(k in names(y.st)){  colnames(y.st[[k]])[2:3]<-paste0(x,".",k,".",colnames(y.st[[k]])[2:3])   }

		Y<-Reduce(function(...) merge(..., all = TRUE), y.st)
		return(Y)

	})
	TAB.ALL<-Reduce(function(...) merge(..., all = TRUE), TAB)

	diff.save<-split(tab.diff.sub[,c("state","diff"),with=FALSE],tab.diff.sub$fraction)
	for(k in names(diff.save)){  colnames(diff.save[[k]])[2]<-paste0(k,".",colnames( diff.save[[k]] )[2] )   }
	DIFF.ALL<-Reduce(function(...) merge(..., all = TRUE), diff.save)

	TAB.ALL<-merge(TAB.ALL,DIFF.ALL)


	### --- PLOTS --- ###
	if(layout=="vertical"){



		plot.diff.state.bar<-ggplot(data=tab.all.state, aes(y=perc,x=fraction,fill=type) )+
			geom_bar(stat="identity")+
			scale_fill_manual(values=col.ab)+
			facet_grid(~state)+
			#scale_x_continuous(breaks=(0:10)/10*100 )+
			#scale_y_discrete(position = "right")+
			theme_classic()+
			theme(strip.background = element_blank(),
				axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1,size=4),
		  		strip.text.x = element_blank(),
		  		panel.spacing = unit(0.1, "lines"),
		  		axis.text.y=element_text(size=5))+
			xlab("")+ylab("")
		#ggsave(plot.diff.state.bar,file=paste0(dir.plot,"Bars.pdf"),width=7,height=1.8)


		plot.diff.state.heat<-ggplot(data=tab.diff.sub,aes(y=fraction.x,x=state,fill=diff))+
			#geom_tile()+
			geom_point(aes(size=abs(diff)),shape=21)+
			#geom_text(aes(label=round(diff,1)),size=2.5)+
			scale_fill_gradient2(name="SAMMY(A%) - HiC(A%)", low="#B73E3E", mid="white", high="#227C70",midpoint = 0)+
			facet_grid(~state,scales="free_x")+
			scale_size(guide = 'none')+
			theme_classic()+
			theme(axis.text.x=element_text(angle=45,size=7,hjust=1),
				strip.background = element_blank(),
		  		strip.text.x = element_blank(),
		  		panel.spacing = unit(0.1, "lines"),
		  		legend.key.height=unit(0.4,"cm"),
		  		legend.key.width=unit(0.2,"cm"),
		  		legend.title = element_text(size=5),legend.text = element_text(size=5))+
			xlab("")+ylab("")

		#ggsave(plot.diff.state.heat,file=paste0(dir.plot,"HeatCircle.pdf"),width=7,height=1.8)
	

		res.plot<-list(barplot=plot.diff.state.bar,heatmap=plot.diff.state.heat,table=TAB.ALL)
		return( res.plot )
		
	} else { 
		return(NULL) 
	}



}




# Should be improved
plot.compartment.summary.chr.replica<-function(E.all,type.plot="Number",col.ab,reference,others,names.new=NULL,span.thick=250){


	E.n<-rbindlist(lapply(E.all,function(x){ return(x$n) }))
	#E.n<-E.n[-which(E.n$fraction==reference),]
	E.n<-E.n[which(E.n$fraction %in% others),]

	tab.chr<-reshape2::melt(tapply(E.n$number,list(E.n$chr,paste0(E.n$fraction)),sum))
	colnames(tab.chr)<-c("chr","fraction","total")

	E.n.perc<-merge(E.n,tab.chr,by=c("chr","fraction"),all.x=TRUE)
	E.n.perc$perc<-E.n.perc$number/E.n.perc$total*100

	# Bars order
	E.n.perc$type<-factor(paste0(E.n.perc$type),levels=c("A-A","B->A","A->B","B-B"))
	E.n.perc$chr<-factor(paste0(E.n.perc$chr),levels=c(paste0("chr",1:22),"chrX","chrY"))

	plot<-NULL

	# Find replicas
	E.n.perc$replica<-sapply(strsplit(paste0(E.n.perc$fraction),split="\\."),function(x) x[3])
	E.n.perc$fraction.group<-sapply(strsplit(paste0(E.n.perc$fraction),split="\\."),function(x) paste0(x[1],".",x[2],collapse="")  )


	if(!is.null(names.new)){
		E.n.perc$fraction<-factor(paste0(E.n.perc$fraction),levels=others,labels=names.new[others])
		tab.chr$fraction<-factor(paste0(tab.chr$fraction),levels=others,labels=names.new[others])
		reference<-names.new[reference]

	} else {
		E.n.perc$fraction<-factor(paste0(E.n.perc$fraction),levels=c(reference,others))
		tab.chr$fractions<-factor(paste0(tab.chr$fraction),levels=c(reference,others))
		
	}




	if(type.plot=="Number"){

		n.max<-max(tapply(E.n$number,list(E.n$fraction,E.n$chr),sum),na.rm=TRUE)
		br.n<-seq(0,(n.max %/% span.thick)*span.thick,by=span.thick)

		plot<-ggplot(data=E.n.perc,aes(x=fraction,y=number))+
			geom_bar(aes(fill=type),stat="identity",position = "stack",width=0.8)+
			theme_classic()+facet_grid(~chr)+
			scale_fill_manual(values=col.ab)+
			theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=7))+
			scale_y_continuous(breaks=br.n,expand=c(0,0,0.1,0))+
			xlab("")+ylab("Number of bins")+labs(subtitle=paste0("Number of bins for each type of compartment, compared to the reference (",reference,")") )
	}

	if(type.plot=="Percentage"){

		plot<-ggplot(data=E.n.perc,aes(x=fraction,y=perc))+
			geom_bar(aes(fill=type),stat="identity",position = "stack",width=0.8)+
			#geom_text(data=unique(tab.chr[,c(1,3)]),aes(x=chr,y=+Inf,label=total))+
			theme_classic()+facet_grid(~chr,scales="free")+
			scale_fill_manual(values=col.ab)+
			theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=7))+
			scale_y_continuous(breaks=(0:10)*10,expand=c(0,0,0,1))+
			xlab("")+ylab("Percentage of bins")+labs(subtitle=paste0("Percentage of bins for each type of compartment, compared to the reference (",reference,")") )

	} 

	if(length(plot)==0) { cat("Only 'Percentage' or 'Number' are valid 'type.plot' options") }

	#ggsave(plot,file=paste0(dir.plot,"ABdifference_resume_",binsize,".pdf"),width=20,height=3)

	return(plot)

}



plot.compartment.summary.genome<-function(E.all,type.plot="Number",col.ab,reference,others,names.new=NULL,span.thick=250){


	E.n<-rbindlist(lapply(E.all,function(x){ return(x$n) }))
	#E.n<-E.n[-which(E.n$fraction==reference),]
	E.n<-E.n[which(E.n$fraction %in% others ),]

	## Add column with the summary for all chromosomes
	E.total.n<-reshape2::melt(tapply(E.n$number,list(value=E.n$value,value.ref=E.n$value.ref,fraction=E.n$fraction,type=E.n$type),sum),value.name="number")
	E.total.n<-data.table(E.total.n[!is.na(E.total.n$number),])
	E.total.n$chr<-"Genome"
	E.total.n<-E.total.n[,colnames(E.n),with=FALSE]

	#E.n<-rbind(E.n,E.total.n)

	tab.tot<-reshape2::melt(tapply(E.total.n$number,list(E.total.n$chr,paste0(E.total.n$fraction)),sum))
	colnames(tab.tot)<-c("chr","fraction","total")

	E.n.perc<-merge(E.total.n,tab.tot,by=c("chr","fraction"),all.x=TRUE)
	E.n.perc$perc<-E.n.perc$number/E.n.perc$total*100

	# Bars order
	E.n.perc$type<-factor(paste0(E.n.perc$type),levels=c("A-A","B->A","A->B","B-B"))
	
	E.n.perc$chr<-droplevels(factor(paste0(E.n.perc$chr),levels=c("Genome",paste0("chr",1:22),"chrX","chrY")))
	tab.tot$chr<-droplevels(factor(paste0(tab.tot$chr),levels=c("Genome",paste0("chr",1:22),"chrX","chrY")))



	plot<-NULL

	if(!is.null(names.new)){
		E.n.perc$fraction<-factor(paste0(E.n.perc$fraction),levels=others,labels=names.new[others])
		tab.tot$fraction<-factor(paste0(tab.tot$fraction),levels=others,labels=names.new[others])
		reference<-names.new[reference]

	} else {
		E.n.perc$fraction<-factor(paste0(E.n.perc$fraction),levels=c(reference,others))
		tab.chr$fractions<-factor(paste0(tab.chr$fraction),levels=c(reference,others))
		
	}


	if(type.plot=="Number"){

		n.max<-max(tapply(E.n$number,list(E.n$fraction,E.n$chr),sum),na.rm=TRUE)
		br.n<-seq(0,(n.max %/% span.thick)*span.thick,by=span.thick)

		plot<-ggplot(data=E.n.perc,aes(x=fraction,y=number))+
			geom_bar(aes(fill=type),stat="identity",position = "stack",width=0.8)+
			theme_classic()+
			#facet_grid(~fraction)+
			scale_fill_manual(values=col.ab)+
			theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=7))+
			scale_y_continuous(breaks=br.n,expand=c(0,0,0.1,0))+
			xlab("")+ylab("Number of bins")+labs(subtitle=paste0("Number of bins for each type of compartment, compared to the reference (",reference,")") )
	}

	if(type.plot=="Percentage"){

		plot<-ggplot(data=E.n.perc,aes(y=fraction,x=perc))+
			geom_bar(aes(fill=type),stat="identity",position = "stack",width=0.8)+
		
			theme_classic()+
		
			scale_fill_manual(values=col.ab)+
			theme(plot.title=element_text(size=6.5),
				plot.subtitle=element_text(size=5.5),
				#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=7)
				axis.text.x = element_text(size=6)
				)+
			scale_x_continuous(breaks=(0:10)*10,expand=c(0,0,0.05,0))+
			ylab("")+xlab("Percentage of bins")+
			labs(
				#title=paste0("Percentage of bins for each type of compartment, compared to the reference (",reference,")"),
				subtitle=paste0("binsize: ",format(as.numeric(binsize),big.mark=",")," bp; number of bins ",format(max(tab.tot$total),big.mark="," ) ) 
			)

	} 

	if(length(plot)==0) { cat("Only 'Percentage' or 'Number' are valid 'type.plot' options") }

	#ggsave(plot,file=paste0(dir.plot,"ABdifference_resume_",binsize,".pdf"),width=5,height=3)

	return(plot)

}





plot.distribution.state.barplot<-function(E,B,state,chr=NULL,col.ab,trans="identity",select.show=NULL,names.new=NULL,sort.name=NULL){

	
	if(length(chr)>0){
		# fix for multiple chr
		B.all<-B[[chr]]
		E.all<-E[[chr]]$res

	} else {

		B.all<-rbindlist(B)
		E.all<-rbindlist(lapply(E,function(x) x$res))
		
	}

	if(!is.null(select.show)){
		E.all<-E.all[E.all$fraction %in% select.show,]
	}

	if(!is.null(names.new)){
		E.all$fraction<-factor(paste0(E.all$fraction),levels=names(names.new), labels=names.new)
		E.all$fraction<-droplevels(E.all$fraction)
	}

	E.all<-E.all[,c("fraction","bin","chr","sign")]

	# Select and aggregate states
	ind.take<-which(paste0(B.all$state) %in% state)
	B.all.sel<-B.all[ind.take,]


	B.all.sel.split<-split(B.all.sel,by="state")
	BB<-lapply(names(B.all.sel.split),function(x){
		dt.x<-B.all.sel.split[[x]]
		dt.x<-dt.x[,c("chr","bin","bp")]
		colnames(dt.x)[3]<-x
		return(dt.x)
	})

	BB.all<-Reduce(function(...) merge(...,by=c("bin","chr"), all = TRUE,sort=FALSE), BB)
	tot.bp.state<-apply(BB.all[,..state],2,sum)

	BB.relfreq<-BB.all
	for(i in 1:length(state)) BB.relfreq[[state[i]]]<-BB.relfreq[[state[i]]]/tot.bp.state[state[i]]*100

	S<-merge(E.all,BB.relfreq,all.x=TRUE,by=c("bin","chr"), all = TRUE,sort=FALSE)
	S$sign<-factor(S$sign,levels=c(-1,1),labels=c("B","A"))

	# Sum up the relative frequency for each state
	Freq.state.all<-rbindlist(lapply(state,function(x){
		dt.x<-S[,c("bin","chr","fraction","sign",x),with=FALSE]
		tab.x<-reshape2::melt(tapply(dt.x[[x]],list(fraction=dt.x$fraction,type=dt.x$sign),sum),value.name = "rel.freq")
		tab.x$state<-x
		return(data.table(tab.x))
	}))



	if(!is.null(sort.name)){
		## Sort by prevalence of state in A name[1]-name[2]
		f1<-Freq.state.all[Freq.state.all$fraction==sort.name[1] & Freq.state.all$type=="A",]
		f2<-Freq.state.all[Freq.state.all$fraction==sort.name[2] & Freq.state.all$type=="A",]

		f1<-f1[match(state,f1$state),]
		f2<-f2[match(state,f2$state),]

		diff<-f1$rel.freq-f2$rel.freq
		names(diff)<-state
		diff<-sort(diff,decreasing=TRUE)

	}

	Freq.state.all$state<-factor(paste0(Freq.state.all$state),levels=names(diff))
	Freq.state.all$fraction<-factor(paste0(Freq.state.all$fraction),levels=rev(select.show))

	plot<-ggplot(data=Freq.state.all,aes(x=rel.freq,y=fraction))+
		geom_bar(aes(fill=type),stat="identity")+
		facet_grid(state~.)+
		scale_x_continuous(expand=c(0,0,0.05,0))+
		scale_fill_manual(values=c(A="#00A19D",B="#FFB344"))+
		theme_classic()+xlab("")+ylab("")
		#scale_y_continuous(expand=c(0.1,0,0.1,0),trans=trans)+
		#
		#facet_wrap(~fraction,scale="free_y",nrow=1)+xlab("")+
		#labs(
		#	title=paste0(nn," - Chromosome ",gsub("chr","",chr)," ( ", format(as.numeric(binsize),big.mark=",",scientific = FALSE), " bp )"),
		#	subtitle=paste0("Compartments: (",reference,") -> (",paste0(levels(S$fraction),collapse=","),")") )+
		#theme( panel.grid.major.y = element_line(colour = "grey",size=0.2),legend.position="none",
		#	plot.title = element_text(size = 10),plot.subtitle = element_text(size = 8))


	ggsave(plot,file=paste0(dir.plot,"State_resume_barsorted_",binsize,".pdf"),width=4,height=8)

	return(plot)

}






plot.distribution.state<-function(E,B,state,chr=NULL,col.ab,reference,trans="identity",select.show=NULL,names.new=NULL){

	
	if(length(chr)>0){
		# fix for multiple chr
		B.all<-B[[chr]]
		E.all<-E[[chr]]$res

	} else {

		B.all<-rbindlist(B)
		E.all<-rbindlist(lapply(E,function(x) x$res))
		
	}

	if(!is.null(select.show)){
		E.all<-E.all[E.all$fraction %in% select.show,]
	}

	if(!is.null(names.new)){
		E.all$fraction<-factor(paste0(E.all$fraction),levels=names(names.new), labels=names.new)
		E.all$fraction<-droplevels(E.all$fraction)
		reference<-names.new[reference]
	}

	E.all<-E.all[,c("fraction","bin","chr","type")]

	# Select and aggregate states
	ind.take<-which(paste0(B.all$state) %in% state)
	B.all.sel<-B.all[ind.take,][, sum(bp),by=list(bin,chr)]
	colnames(B.all.sel)<-c("bin","chr","bp")


	S<-merge(E.all,B.all.sel,by=c("bin","chr"),all.x=TRUE,)
	S.text<-reshape2::melt(table(S$fraction,S$type))
	colnames(S.text)<-c("fraction","type","number")
	S.text$bp<- +Inf

	nn<-paste(state,collapse=" ; ")

	if(!is.null(select.show)){
		S$fractions<-factor(paste0(S$fractions),levels=select.show,labels=select.show)
	}


	plot<-ggplot(data=S,aes(x=type,y=bp+1))+
		
		geom_text(data=S.text,aes(label=paste0("(n=",format(number,big.mark=","),")")),size=1.8,vjust=1.5)+
		
		geom_violin(aes(fill=type),scale="width",colour=NA,alpha=0.5)+
		geom_boxplot(outlier.size=0.2,outlier.color="darkgray",width=0.1,size=0.2,fill=NA)+
		theme_classic()+
		scale_y_continuous(expand=c(0.1,0,0.1,0),trans=trans)+
		scale_fill_manual(values=col.ab)+
		facet_wrap(~fraction,scale="free_y",nrow=1)+xlab("")+
		labs(
			title=paste0(nn," - Chromosome ",gsub("chr","",chr)," ( ", format(as.numeric(binsize),big.mark=",",scientific = FALSE), " bp )"),
			subtitle=paste0("Compartments: (",reference,") -> (",paste0(levels(S$fraction),collapse=","),")") )+
		theme( panel.grid.major.y = element_line(colour = "grey",size=0.2),legend.position="none",
			plot.title = element_text(size = 10),plot.subtitle = element_text(size = 8))


	#ggsave(plot,file=paste0(dir.plot,"State_resume_",binsize,".pdf"),width=13,height=3)

	return(plot)

}





plot.distribution.state.paper<-function(E,B,state,chr=NULL,col.ab,reference,select.show=NULL,names.new=NULL,binsize,breaks,return.plot=FALSE){

	
	if(length(chr)>0){
		# fix for multiple chr
		B.all<-B[[chr]]
		E.all<-E[[chr]]$res

	} else {

		B.all<-rbindlist(B)
		E.all<-rbindlist(lapply(E,function(x) x$res))
		
	}

	if(!is.null(select.show)){
		E.all<-E.all[E.all$fraction %in% select.show,]
	}

	if(!is.null(names.new)){
		E.all$fraction<-factor(paste0(E.all$fraction),levels=names(names.new), labels=names.new)
		E.all$fraction<-droplevels(E.all$fraction)
		reference<-names.new[reference]
	}

	E.all<-E.all[,c("fraction","bin","chr","type")]

	# Select and aggregate states
	ind.take<-which(paste0(B.all$state) %in% state)
	B.all.sel<-B.all[ind.take,][, sum(bp),by=list(bin,chr)]
	colnames(B.all.sel)<-c("bin","chr","bp")

	tot.state<-sum(B.all.sel$bp)


	S<-merge(E.all,B.all.sel,by=c("bin","chr"),all.x=TRUE,)
	S.text<-reshape2::melt(table(S$fraction,S$type))
	colnames(S.text)<-c("fraction","type","number")
	S.text$bp<- +Inf

	nn<-paste(state,collapse=" ; ")

	if(!is.null(select.show)){
		S$fraction<-factor(paste0(S$fraction),levels=select.show,labels=select.show)
	}


	## plot layout parameters
	th<-10000
	bs<-as.numeric(binsize)
	ylim<-c(0,ceiling(max(S$bp)/th)*th)

	breaks<-c(0,0.0001,0.001,0.01,0.1,1)
	#breaks<-c(0:10)/10

	S.text<-S.text[S.text$number>0,]


	RES<-list(data=S,text=S.text,parameters=list(th=th,bs=bs,ylim=ylim,breaks=breaks) )

	if(return.plot){
		
		plot<-ggplot(data=S,aes(x=type,y=bp))+
			
			geom_text(data=S.text,aes(label=paste0("(n=",format(number,big.mark=","),")")),size=1.8,vjust=1.5)+
			
			geom_violin(aes(fill=type),scale="width",colour=NA,alpha=0.5)+
			geom_boxplot(outlier.size=0.2,outlier.shape=NA,width=0.1,size=0.4,fill=NA)+
			theme_classic()+
			scale_y_continuous(expand=c(0.1,0,0.1,0),trans= scales::pseudo_log_trans(base=2),breaks=breaks*bs,labels=paste0(breaks*100,"%") )+
			scale_fill_manual(values=col.ab)+
			
			facet_grid(~fraction,scales="free_x",space = "free_x")+xlab("")+
			labs(
				title=paste0(nn," - Chromosome ",gsub("chr","",chr)," ( ", format(as.numeric(binsize),big.mark=",",scientific = FALSE), " bp )"),
				subtitle=paste0("Compartments: (",reference,") -> (",paste0(levels(S$fraction),collapse=","),")") )+
			theme( panel.grid.major.y = element_line(colour = "grey",size=0.2),legend.position="none",
				plot.title = element_text(size = 10),plot.subtitle = element_text(size = 8))

		RES<-c(RES,plot=list(plot))

	}
	#ggsave(plot,file=paste0(dir.plot,"State_resume_",binsize,".pdf"),width=13,height=3)

	return(RES)

}





plot.enhrichment.histone.chr<-function(E,TT,type="tile",chr=NULL,select.show=NULL,col.ab=NULL,marks,fun=mean,reference,names.new=NULL){


	if(!is.null(chr)){
		E<-E[chr]
		TT<-TT[chr]
	}

	TT.all<- unlist(GRangesList(TT))
	E.all<-rbindlist(lapply(E,function(x) x$res))

	#E.all$type<- c("B","A")[(E.all$sign==1)+1]

	if(!is.null(select.show)){
		E.all<-E.all[E.all$fraction %in% select.show,]
	}

	TT.all<-as.data.table(mcols(TT.all))
	take<-c("bin","Id","chr",rev(paste0(marks)))
	TT.all<-TT.all[,..take]

	## select bins
	bin.select<-unique(paste(E.all$chr,E.all$bin,sep=":"))
	TT.all<-TT.all[paste(TT.all$chr,TT.all$bin,sep=":") %in% bin.select,]

	## Negative values are not allowed (add a costant and scale)
	TT.all[,(marks) := lapply(.SD, function(x){  scale(x,center=TRUE,scale=TRUE)  } ), .SDcols = marks ]

	take<-c("fraction","bin","Id","chr","type")
	E.all<-E.all[,..take]


	TT.melt<-reshape2::melt(TT.all,c("bin","chr","Id"))
	
	exp.chr<-reshape2::melt(tapply(TT.melt$value,list(mark=TT.melt$variable),fun),value.name="exp")
	exp.chr<-data.table::data.table(exp.chr)

	## median/mean observed histone values by compartments/chr/histone mark
	E.all.merge<-merge(TT.melt,E.all,by=c("bin","chr","Id"),allow.cartesian=TRUE)

	## Show only median/mean enrichment
	if(type=="tile"){
		
		obs.chr<-reshape2::melt(tapply(E.all.merge$value,list(mark=E.all.merge$variable,type=E.all.merge$type,fraction=E.all.merge$fraction),fun),value.name="obs")
		obs.chr<-data.table::data.table(obs.chr)

		## data.table for plot
		OE<-merge(obs.chr,exp.chr,all.x=TRUE)
		OE$enrich<- OE$obs - OE$exp

		if(!is.null(names.new)){
			OE$fraction<-factor(paste0(OE$fraction),levels=rev(names(names.new)), labels=rev(names.new))
			OE$fraction<-droplevels(OE$fraction)
			reference<-names.new[reference]
		}


		#OE$type<-factor(paste0(OE$type),levels=c("B","A"))
		OE$mark<-factor(paste0(OE$mark),levels=marks,labels=names(marks))


		PLOT<-ggplot(data=OE,aes(x=type,y=fraction,fill=enrich))+
			geom_tile(colour="white",size=0.5)+
			geom_text(aes(label=round(enrich,1)),size=2.5,colour="white")+
			theme_classic()+
			facet_grid(~mark)+coord_equal()+
			scale_fill_gradient2(name="log2(Observed)-log2(Expected)",na.value = 'white')+
			theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),panel.spacing = unit(2, "lines"),
				legend.direction = 'horizontal',legend.title=element_text(size=9),legend.position=c(1,1.5),legend.justification=c(1,1),legend.text = element_text(size=6),legend.key.height= unit(0.3, 'cm'))+
			labs(title=paste0("Signal enrichment (",fun@generic[[1]],")"),
				subtitle=paste0("Reference: ",reference," - N.bins: ", format(nrow(TT.all),big.mark=",")," - Bin size: ",format(TT[[1]][1]@ranges@width,big.mark=","),"bp"),x="",y="")
	} 

	## Show all distributions
	if(type=="box"){

		colnames(E.all.merge)<-c("bin","chr","Id","mark","value","fraction","type")

		if(!is.null(names.new)){
			E.all.merge$fraction<-factor(paste0(E.all.merge$fraction),levels=names(names.new), labels=names.new)
			E.all.merge$fraction<-droplevels(E.all.merge$fraction)
			reference<-names.new[reference]
		}

		E.all.merge$mark<-factor(paste0(E.all.merge$mark),levels=marks,labels=names(marks))
		exp.chr$mark<-factor(paste0(exp.chr$mark),levels=marks,labels=names(marks))

		PLOT<-vector("list",length(marks))
		names(PLOT)<-names(marks)
		
		for(k in 1:length(marks)){

			#E.sub.k<-
			#E.sub.k$mark<-droplevels(E.sub.k$mark)

			PLOT[[k]]<-ggplot(data=E.all.merge[mark==names(marks)[k]],aes(x=type,y=value))+
				geom_hline(data=exp.chr[mark==names(marks)[k]],aes(yintercept=exp),col="darkgray")+
				geom_violin(aes(fill=type),colour=NA,size=0.2,alpha=0.6)+
				geom_boxplot(fill=NA,colour="black",size=0.3,width=0.1,outlier.size = 0.01,outlier.alpha=0.5)+
				theme_classic()+
				facet_grid(mark~fraction)+
				scale_fill_manual(values=col.ab)+
				labs(title=paste0("Scaled signal distribution for ",names(marks)[k]),x="",y="scaled(score)",
					subtitle=paste0("Reference: ",reference," - N.bins: ", format(nrow(TT.all),big.mark=",")," - Bin size: ",format(TT[[1]][1]@ranges@width,big.mark=","),"bp"))
		
		}


 
	}

	#ggsave(PLOT[[k]],file=paste0(dir.plot,"Prova3.pdf"),width=10,height=3.5)

	return(PLOT)

}




get.switch.agreement<-function(E,B,reference,others){


	select<-c(unique(unlist(others)),reference)
	
	M.all<-sign(do.call("cbind",E))
	M.all<-data.table::data.table(M.all)

	agree<-lapply(others,function(x,m){
		tab.x<-sapply(m[,..x],function(fr,ref){ as.numeric(fr!=ref) },ref=m[[reference]])
		vec.agr<-as.numeric(apply(tab.x,1,sum)==length(x))
	},m=M.all[,..select])

	
	agree<-data.table::data.table(do.call("cbind",agree))
	ind.take<-which(apply(agree,1,sum)!=0)
	
	agree$type<- c("B->A","A->B")[as.numeric(M.all[[reference]]==1)+1]
	mcols(B)<-cbind(M.all,agree,mcols(B))


	return(B[ind.take])
}

#https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E055_15_coreMarks_hg38lift_dense.bed.gz

## two or three elements: the central one in select variabe is the reference
# To be improved
plot.alluvional.AB.all<-function(E,select,col.ab){

	if(length(select)<2 | length(select)>3){ warning("Only three classes are supported"); return(NA) }

	reference<-select[2]
	#tab<-data.table::data.table(reshape2::melt(E))
	A<-rbindlist(lapply(1:length(E),function(i){
		x<-E[[i]]
		mat<-do.call("cbind",x)

		tab<-data.table(apply(mat,2,function(x) c("B","A")[as.numeric(x>0)+1] ))
		tab$chr<-names(E)[i]
		tab$bin<-rownames(mat)

		return(tab)
	}))

	Tab<-as.data.table(table(A[,select,with=FALSE]))
	Tab$Type.ref<-Tab[[reference]]	
	D<-to_lodes_form(Tab,axes=1:length(select))

	D$Freq<-D$N/ nrow(A)*100
	#D$Type<-paste0(D$Type.ref,"->",D$stratum)
	#D$Type[D$Type.ref==D$stratum]<-gsub("->","-",D$Type[D$Type.ref==D$stratum])

	


	col.all<-c(A=col.ab[["A-A"]],B=col.ab[["B-B"]])
	
	ww<-1/3
 	plot<-ggplot(data=D,aes(x = x, y=Freq,stratum = stratum, alluvium = alluvium)) +
  		geom_flow(aes(fill=Type.ref),alpha=0.3,width=ww)+
  		geom_stratum(aes(fill = stratum),colour=NA,width=ww)+
  		theme(legend.position = "bottom") +
  		#geom_text(stat = "stratum", aes(label = round(Freq) ),size=3)+
  		scale_fill_manual(values=col.all)+#scale_colour_manual(values=col.all)+
  		theme_classic()+
  		scale_x_discrete(limits = select, expand = c(.1, .1))+
  		scale_y_continuous(breaks=(0:10)*10,expand=c(0,0,0.05,0))+
  		theme(legend.position="none")+xlab("")+ylab("Percentage of bins")

	#ggsave(plot,file=paste0(prefix.path,"Plot/Alluvional_all8.pdf"),width=3,height=3)

	return(list(plot=plot,table=Tab))

}
