# Note: parameter names with "_" are from calder

# trans.01: if the matrix is normalize with ICE, multiply for a constant and allow the minimum number be ~ 1 
# trans.log2: transform in log2 (if trans.01 it's not raccomended)
read.hic.calder<-function(dir.file,name.file,bin.size,trans.01=TRUE,trans.log2=TRUE){

	#see on calder: CALDER::contact_mat_processing(contact_mat_file, bin_size=bin_size)

	cat("\n",name.file,": ")
	## Fixed parameters (by CALDER)
	sub_domains=TRUE
	save_intermediate_data=FALSE
	zero_ratio=0.1
	const<-10^(1:10)
	

	contact_mat_file<-paste0(dir.file,name.file)
	#chr<-gsub("chr","",chr)
	bin_size_initial<-as.numeric(bin.size)

	# Read matrix
	cat("read.. ")
	combined_xk_oe_raw <- data.table::fread(contact_mat_file) # only upper triangular matrix
	ind.ord<-order(combined_xk_oe_raw[,1],combined_xk_oe_raw[,2])
	combined_xk_oe_raw<-combined_xk_oe_raw[ind.ord,]

	# Retrive indeces using the start coordinates of bins and transform in 1-based
	combined_xk_oe_raw[,1]<-(combined_xk_oe_raw[,1]/bin_size_initial)+1
	combined_xk_oe_raw[,2]<-(combined_xk_oe_raw[,2]/bin_size_initial)+1
	
	combined_xk_oe_raw <- subset(combined_xk_oe_raw, !is.na(V3))
	combined_xk_oe <- combined_xk_oe_raw
	colnames(combined_xk_oe) = c('pos_1', 'pos_2', 'val')   

	# Create sparse matrix
	oe_size <- max(max(combined_xk_oe[[1]]), max(combined_xk_oe[[2]]))

	cat("sparse format.. ")
	mat_oe_sparse <- Matrix::Matrix(0, nrow=oe_size, ncol=oe_size)
	mat_oe_sparse[cbind(combined_xk_oe[[1]], combined_xk_oe[[2]])] <- combined_xk_oe[[3]]
	rownames(mat_oe_sparse) <- colnames(mat_oe_sparse) <- as.character( 1:nrow(mat_oe_sparse) )
	mat_oe_sparse <- Matrix::forceSymmetric(mat_oe_sparse, uplo='U')

	cat("remove blank.. ")
	A_oe <- CALDER::remove_blank_cols(mat_oe_sparse, sparse=TRUE, ratio=zero_ratio) ## has the same rows/cols as A

	A_oe_ret<-A_oe
	
	if(trans.01){	
		cat("transform.. ")
		mina<-min(summary(A_oe_ret)$x)
		check<-floor((mina * const ))
		sel.const<-const[min(which(check>0))-1]
		A_oe_ret<-A_oe_ret * sel.const
	}

	if(trans.log2){
		cat("log2.. ")
		A_oe_ret <- log2(A_oe_ret + 1) 
	}

	A_oe_ret<-as.matrix(A_oe_ret)

	cat("\n\n ")
	return(A_oe_ret)
}

## Note: possible riscale [0,1] as for hic 
read.SAMMY.calder<-function(dir.file,name.file,bin.size,select,bin.select=NULL,frac.select=NULL,metric="euclidean"){

	cat("\t",select,": read.. ")

	TR.chr<-get(load(paste0(dir.file,name.file)))
	TR.chr.sel<-data.table::data.table(do.call("cbind",(lapply(TR.chr[[select]],function(x) x$score))))

	if(!is.null(frac.select)){
		TR.chr.sel<-TR.chr.sel[,frac.select,with=FALSE]
	}

	GR.chr.sel<-granges(TR.chr[[select]][[1]])
	GR.chr.sel$Id<-1:length(GR.chr.sel)	

	# only for comparing with HiC matrix
	if(!is.null(bin.select)){
		## Ids compatible with calder (starts are 0-based)
		# all TR object have the same length
		Id<-( (start(TR.chr[[1]][[1]])-1  ) / as.numeric(bin.size) ) +1

		## Common bins between HiC and SAMMY/Chip-seq
		cat("select.. ")
		ind.take<-match(bin.select,Id)
		TR.chr.sel<-TR.chr.sel[ind.take,]
		rownames(TR.chr.sel)<-rownames(TR.chr.sel)<-bin.select

		GR.chr.sel<-GR.chr.sel[ind.take,]
		#names(GR.chr.sel)<-bin.select		

	} else { bin.select<-1:nrow(TR.chr.sel) }


	## Calculate eucledean distance between pairs of points (i.e., bins)
	# Each point is define in the n-dimensional space, where n is 3,4, or 6 based on the number of fractions or Chip-seq experiments
	cat("pairwise ",ncol(TR.chr.sel),"-dimensional distances.. ",sep="")	
	A<-as.matrix(dist(TR.chr.sel,method=metric))
	rownames(A)<-colnames(A)<-bin.select

	cat("\n")
	return(list(A=A,Bin=GR.chr.sel))
}



## Note: possible riscale [0,1] as for hic 
read.marks<-function(dir.file,name.file,bin.size,select,bin.select){

	#cat("\n",select,": read.. ")


	TR.chr<-get(load(paste0(dir.file,name.file)))
	TR.chr.sel<-data.table::data.table(do.call("cbind",(lapply(TR.chr[[select]],function(x) x$score))))
	
	
	## Ids compatible with calder (starts are 0-based)
	# all TR object have the same length
	Id<-( (start(TR.chr[[1]][[1]])-1  ) / as.numeric(bin.size) ) +1

	## Common bins between HiC and SAMMY/Chip-seq
	#cat("select.. ")
	ind.take<-match(bin.select,Id)
	TR.chr.sel<-TR.chr.sel[ind.take,]
	TR.chr.sel$Id<-bin.select

	#cat("\n")
	return(TR.chr.sel)
}




# cor.cor: if TRUE compute two times the correlation
# trans.atanh: trasnform the cor (or cor.cor) with the inverse hyperbolic tangent
# const: avoid +/- Inf in inverse hyperbolic tangent
compute.fastcor.calder<-function(A,cor.cor=TRUE,trans.atanh=TRUE,const=1+1E-7){

	cat("\t1.cor.. ")
	cA<-CALDER::fast_cor(A)

	if(cor.cor){
		cat("2.cor.. ")
		ccA<- CALDER::fast_cor(cA)
	} else { ccA<-cA }

	if(trans.atanh){
		cat("inv.hyper.tangent.. ")
		accA<- atanh( ccA / const)
	} else { accA<- ccA }

	cat("\n")
	return(accA)

}


# A: should be the ccaA_oe_compressed_log
get.blocks.calder<-function(A,bin.size,chr,window.sizes = 3){

	p.th <- ifelse(as.numeric(bin.size) < 40000, 0.05, 1)
	
	# change the indices to be sure of the consistences of boundary predictions
	info.index<-data.table(Id=rownames(A),index=paste0(1:nrow(A)))
	rownames(A)<-colnames(A)<-1:nrow(A)

	#chr_name = paste0("chr", chr)
	TD.out<- CALDER::generate_compartments_bed(input_mat=A,chr=gsub("chr","",chr),p_thresh=p.th, bin_size=as.numeric(bin.size), window.sizes = window.sizes, out_file_name=NULL, stat_window_size = NULL )



	blocks<-lapply(1:nrow(TD.out$domain),function(i,D,info) {
		x<-D[i,]
		cl.id<-info$Id[x$from.id:x$to.id]
		
		return(cl.id)
	},D=TD.out$domain,info=info.index)
	names(blocks)<-1:length(blocks)

	return(blocks)

}

# A: should be the contact/distance matrix
cor.trend.blocks<-function(A,blocks,lag=4,trans.atanh=TRUE,scale=TRUE,const=1+1E-7,metric="mean"){

	# Summarize by blocks (it takes the rownames of the blocks)
	cat("\tSummarize by blocks.. ")
	B <- CALDER::HighResolution2Low_k_rectangle(mat=A,row_split=blocks,col_split=blocks,sum_or_mean =metric)
	rownames(B)<-colnames(B)<-names(blocks)

	# Compute enrichment trends at different lags 
	cat("Trend.. ")
	n.block<-length(blocks)
	T.lags <- lapply( 1:lag, function(v,mat,n) { 
		n<-nrow(B)
		1 * (mat[, -(1:v)] > mat[, - n - 1 + (v:1)])  
	},mat=B,n=n.block)
	T <- do.call(cbind, T.lags)
	

	cat("Corr.. ")
	cT<- CALDER::fast_cor(t(T))

	if(trans.atanh){
		cat("Trans atanh.. ")
		acT<-atanh(cT/const)
	} else { acT<-cT }

	if(scale){
		cat("Scale..\n")
		acT.scaled<-scale(acT)
	} else {  acT.scaled<-acT }
	

	return(acT.scaled)

}


# T: trend matrix
# n.comp: number of principal compnent
# const.comp: devide the pca 2:ncomp by the constany (in this way the first remain predominant)
# error: standard deviation of random error to be added to duplicated rows in the data
get.subcompartment.calder<-function(T,blocks,chr,active,tracks,n.comp=10,const.comp=5,error=0.0001){


	## Take first 10 principal component
	PC.comp <- CALDER::get_PCs(T, which=1:n.comp)
	PC.comp[,2:n.comp] <- PC.comp[,2:n.comp]/const.comp

	## ----------- Note ----------- ##
	# Some rows have the exact same principal components
	# Add some noise to these rows (otherwise the bisection in k-means cluster is not feasible)
	# Error reported: more cluster centers than distinct data points
	dup<-setdiff(rownames(PC.comp),rownames(unique(PC.comp)))
	set.seed(1234)
	PC.comp[dup,]<-PC.comp[dup,]+rnorm(length(dup)*n.comp,mean=0,sd=error)
	## ---------------------------- ##

	## First PCA should have the same sign of active marker
	# active marker score by block
	bl.active<-sapply(1:length(blocks),function(i,bl,dt,ac) mean(dt[[active]][match(bl[[i]],dt$Id)]) ,bl=blocks, dt=tracks,ac=active )
	PC.comp[,1]<-PC.comp[,1]*sign(cor(PC.comp[,1],bl.active))


	## Complete k(=2)-iterative clustering, with eucledean distance return a hclust dendogram object)
	H.k2 <- CALDER::bisecting_kmeans(PC.comp)
	
	## Reorder blocks using the first projected linear component
	# Non-linear projection using the first two components
	new.pc1 <- CALDER::project_to_major_axis(PC.comp)
	ord.block<-CALDER::get_best_reorder(hc_hybrid_x_pro=H.k2, x_pro=new.pc1$x_pro) 

	H.k2.ord <- dendextend::rotate(x=H.k2, order=ord.block)

	## vector of 
	AB.sub<-CALDER::get_cluser_levels(H.k2.ord, k_clusters=Inf, balanced_4_clusters=FALSE)$cluster_labels
	
	AB.sub.dt<-data.table::data.table(
		chr=chr,
		bin=unlist(blocks),		
		block=rep(names(blocks),sapply(blocks, length)),
		sub=rep(AB.sub[names(blocks)],sapply(blocks, length)) )

	AB.sub.dt$sub.2<-substr(AB.sub.dt$sub,start=1,stop=1)
	AB.sub.dt$sub.4<-substr(AB.sub.dt$sub,start=1,stop=3)
	AB.sub.dt$sub.8<-substr(AB.sub.dt$sub,start=1,stop=5)

	info.pca<-data.table::data.table(
		block=names(blocks),
		sub=AB.sub[names(blocks)],
		sub.2=substr(AB.sub[names(blocks)],start=1,stop=1),
		sub.4=substr(AB.sub[names(blocks)],start=1,stop=3),
		sub.8=substr(AB.sub[names(blocks)],start=1,stop=5),
		pc1=PC.comp[,1],pc2=PC.comp[,2],
		new.pc1=new.pc1$x_pro
	)

	return(list(Bin=AB.sub.dt,Block=info.pca,Dendro=H.k2.ord))

}










### --- Plot --- ###


plot.pca12<-function(AB.obj,select,colour ){



	dt.pca<-rbindlist(lapply(select,function(x){ 
		y<-AB.obj[[x]]$Block
		y$protocol<-x
		return(y)
	}))
	dt.pca$protocol<-factor(paste0(dt.pca$protocol),levels=select )


	tab.nblock<-table(dt.pca$protocol)
	dt.nblock<-data.table(protocol=names(tab.nblock),nblock=as.numeric(tab.nblock),pc1=0,pc2=0)
	
	dt.nblock$protocol<-factor(paste0(dt.nblock$protocol),levels=select )

	br.x<- seq(from=-50,to=50,by=5)
	br.y<- seq(from=-10,to=10,by=2)



	plot.pca<-ggplot(data=dt.pca,aes(x=pc1,y=pc2))+
		geom_hline(yintercept=0,size=0.2)+geom_vline(xintercept=0,size=0.2)+
		geom_point(aes(colour=sub.8),alpha=0.5,size=1)+
		geom_label(data=dt.nblock,aes(label=nblock),size=3)+
		facet_wrap(~protocol,ncol=5,scales="free")+
		scale_colour_manual(values=colour)+
		scale_x_continuous(breaks=br.x,expand=c(0,0.1,0,0.1))+
		scale_y_continuous(breaks=br.y)+
		theme_classic()+theme(
			panel.grid.major = element_line(size = 0.4))

	#ggsave(plot.pca,file=paste0(dir.plot,prefix.name,"pca12.pdf"),width=6,height=10)

	return(plot.pca)

}




plot.eigenvector.state<-function(E,B,reference,col.ab,col.state,name.sorted,name.new=NULL){


	## Plot same bin 
	select.bin<-unique(E$res$bin)
	select.bin<-sort(select.bin)
	B.sel<-B[which(B$bin %in% select.bin),]
	
	ind.ref<-which(E$res$fraction==reference)
	chr<-unique(B$chr)

	## Order chromatin states (based on the order of col.state variable)
	B.sel$state<-factor(paste0(B.sel$state),levels=names(col.state))

	## Order fractions (based on name.sorted variable)
	if(!is.null(name.new)){
		E$res$fraction<-factor(paste0(E$res$fraction),levels=c(name.sorted,reference),labels=name.new[c(name.sorted,reference)])
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
	
	#ggsave(plot.state,file=paste0(dir.plot,"Eigenvectors_",chr,"_",binsize,".pdf"),width=15,height=2.5)

	#ggsave( (plot.ref / plot_spacer() / plot.AB ) / plot.state  + plot_layout(height=c(1.5,-0.5,5.5,2.5)),file=paste0(dir.plot,"Eigenvectors_Allinone_",chr,"_",binsize,".pdf"),width=18,height=7)

	return(list(ref=plot.ref,others=plot.AB,state=plot.state))
}



## ---- Utility functions: plots --- ###


cov.region.by.bin<-function(bin,peak,name=""){
	
	# find overalp 
	ind<-findOverlaps(bin,peak)

	# Fast!
	peak.int<-pintersect(bin[ind@from],peak[ind@to])
	cov<-tapply(peak.int@ranges@width,ind@from,sum)
	
	# create data.table with multiple information
	dt.cov<-data.table(bin=bin$Id,Flag=0,bp.tot=0)
	dt.cov$bp.tot[as.numeric(names(cov))]<-unlist(cov)
	dt.cov$bp.Mb<-(dt.cov$bp/bin@ranges@width)*1000
	dt.cov$Flag[dt.cov$bp.tot>0]<-1

	# rename columns
	colnames(dt.cov)[-1]<-paste0(name,".",colnames(dt.cov)[-1])

	return(dt.cov)
}


get.subcompartment.bin.chr<-function(chr,name.file,dir.file,type="sub.4",select){
	#MergedDT = Reduce(function(...) merge(..., all = TRUE), List_of_DTs)

	cat("-> ",chr,": ")
	take<-c("chr","bin",type)

	## -- Load subcompartments -- ##
	name.chr<-name.file[grep(paste0(chr,"\\."),name.file)]
	Sub.chr<-get(load(paste0(dir.file,name.chr)))
	
	Sub.chr<-Sub.chr[names(Sub.chr) %in% c(select,"Bin") ]
	Bin.chr<-Sub.chr$Bin
		
	
	## -- Info on DMSO and JQ1 -- ##
	A.sub<-lapply(select,function(x,SS){ 
		X<-SS[[x]]$Bin[,take,with=FALSE] 
		colnames(X)[colnames(X)==type]<-x
		return(X)
	},SS=Sub.chr)

	Info<-Reduce(function(...) merge(..., all.x = TRUE,sort=FALSE), A.sub)



	return(Info)
}




get.subcompartment.pairscore.bin.chr<-function(chr,name.file,dir.file,type="sub.4",select,new.names){
	#MergedDT = Reduce(function(...) merge(..., all = TRUE), List_of_DTs)

	cat("-> ",chr,": ")
	take.bin<-c("chr","bin","block",type)
	take.block<-c("block",type,"new.pc1")

	## -- Load subcompartments -- ##
	name.chr<-name.file[grep(paste0(chr,"\\."),name.file)]
	Sub.chr<-get(load(paste0(dir.file,name.chr)))
	
	Sub.chr<-Sub.chr[names(Sub.chr) %in% c(select,"Bin") ]
	Bin.chr<-Sub.chr$Bin
		
	
	## -- Info on DMSO and JQ1 -- ##
	A.sub<-lapply(select,function(x,SS){ 

		# Take blocks and bins
		Bl.x<-SS[[x]]$Block[,take.block,with=FALSE]
		Bi.x<-SS[[x]]$Bin[,take.bin,with=FALSE]

		# Rank blocks
		Bl.x$rank<-rank(Bl.x$new.pc1)
		BiBl.x<-merge(x=Bi.x,y=Bl.x,by=c("block",type),all.x=TRUE,sort=FALSE)
			
		return(BiBl.x)
	},SS=Sub.chr)
	names(A.sub)<-select


	Info<-merge(x=A.sub[[select[1]]],y=A.sub[[select[2]]],by=c("chr","bin"),sort=FALSE,suffixes=paste0(".",new.names) )

	#Info<-Reduce(function(...) merge(..., all.x = TRUE,sort=FALSE ), A.sub)
	
	return(Info)
}





plot.compartment.summary.chr<-function(E.all,type.plot="Number",col.ab,reference){

	#E.n<-#rbindlist(lapply(E.all,function(x){ return(x$n) }))
	#E.n<-E.n[-which(E.n$fraction==reference),]

	tab.chr<-as.data.table(table(E.all[,c("chr",reference),with=FALSE]))
	colnames(tab.chr)<-c("chr","type","number")

	vec.total<-tapply(tab.chr$number,tab.chr$chr,sum)

	tab.chr$total.chr<-vec.total[match(tab.chr$chr,names(vec.total))]
	tab.chr$perc<-tab.chr$number/tab.chr$total.chr*100

	# Bars order
	tab.chr$type<-factor(paste0(tab.chr$type),levels=rev(names(col.ab)))
	tab.chr$chr<-droplevels(factor(paste0(tab.chr$chr),levels=c(paste0("chr",1:22),"chrX","chrY")))

	plot<-NULL
	
	if(type.plot=="Number"){

		plot<-ggplot(data=tab.chr,aes(x=chr,y=number))+
			geom_bar(aes(fill=type),stat="identity",position = "stack",width=0.8)+
			theme_classic()+
			scale_fill_manual(values=col.ab)+
			theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6),legend.position="none")+
			scale_y_continuous(breaks=(0:10)*100,expand=c(0,0,0.1,0))+
			xlab("")+ylab("Number of bins")#+labs(subtitle=paste0("Number of bins for each type of compartment (",reference,")") )
	}

	if(type.plot=="Percentage"){

		plot<-ggplot(data=tab.chr,aes(x=chr,y=perc))+
			geom_bar(aes(fill=type),stat="identity",position = "stack",width=0.8)+
			geom_text(data=tab.chr,aes(x=chr,y=+Inf,label=format(total.chr,big.mark=",")),angle=90,size=1.8,hjust=1.2,fontface=3,colour="darkgray")+
			theme_classic()+
			scale_fill_manual(values=col.ab)+
			theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6),legend.position="none")+
			scale_y_continuous(breaks=(0:10)*10,expand=c(0,0,0.15,0))+
			xlab("")+ylab("Percentage of bins")#+labs(subtitle=paste0("Percentage of bins for each type of compartment (",reference,")") )

	} 

	if(length(plot)==0) { cat("Only 'Percentage' or 'Number' are valid 'type.plot' options") }

	#ggsave(plot,file=paste0(dir.plot,"ABdifference_resume_",binsize,".pdf"),width=13,height=3)

	return(plot)

}


plot.compartment.summary<-function(E.all,type.plot="Number",col.ab,reference){

	#E.n<-#rbindlist(lapply(E.all,function(x){ return(x$n) }))
	#E.n<-E.n[-which(E.n$fraction==reference),]

	tab.chr<-as.data.table(table(E.all[,c(reference),with=FALSE]))
	colnames(tab.chr)<-c("type","number")

	vec.total<-sum(tab.chr$number)

	tab.chr$total.chr<-vec.total
	tab.chr$perc<-tab.chr$number/tab.chr$total.chr*100

	# Bars order
	tab.chr$chr<-"Genome"
	tab.chr$type<-factor(paste0(tab.chr$type),levels=rev(names(col.ab)))
	
	plot<-NULL
	
	if(type.plot=="Number"){

		plot<-ggplot(data=tab.chr,aes(x=chr,y=number))+
			geom_bar(aes(fill=type),stat="identity",position = "stack",width=0.8)+
			theme_classic()+
			scale_fill_manual(values=col.ab)+
			theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6),legend.position="none")+
			scale_y_continuous(breaks=(0:10)*100,expand=c(0,0,0.1,0))+
			xlab("")+ylab("Number of bins")#+labs(subtitle=paste0("Number of bins for each type of compartment (",reference,")") )
	}

	if(type.plot=="Percentage"){

		plot<-ggplot(data=tab.chr,aes(x=chr,y=perc))+
			geom_bar(aes(fill=type),stat="identity",position = "stack",width=0.8)+
			geom_text(data=tab.chr,aes(x=chr,y=+Inf,label=format(total.chr,big.mark=",")),angle=90,size=1.8,hjust=1.2,fontface=3,colour="darkgray")+
			theme_classic()+
			scale_fill_manual(values=col.ab)+
			theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6),legend.position="none")+
			scale_y_continuous(breaks=(0:10)*10,expand=c(0,0,0.15,0))+
			xlab("")+ylab("Percentage of bins")#+labs(subtitle=paste0("Percentage of bins for each type of compartment (",reference,")") )

	} 

	if(length(plot)==0) { cat("Only 'Percentage' or 'Number' are valid 'type.plot' options") }

	#ggsave(plot,file=paste0(dir.plot,"ABdifference_resume_",binsize,".pdf"),width=13,height=3)

	return(plot)

}



plot.compartment.alluvional<-function(E.all,col.ab,select){

	
	
  	A<-table(E.all[,select,with=FALSE])
	B<-as.data.table(A)
	B$Type.ref<-B[[select[2]]]
	B$Freq<-B$N/sum(B$N)*100
	#B$Freq<-B$Freq/sum(B$Freq)*100
	D<-to_lodes_form(B,axes=1:3)
	
  	plot<-ggplot(data=D,aes(x = x, y=Freq,stratum = stratum, alluvium = alluvium)) +
  		geom_flow(aes(fill = Type.ref),alpha=0.4,width=1/5,size=0)+
  		geom_stratum(aes(fill=stratum),width=1/5,size=0)+
  		theme(legend.position = "bottom") +
  		scale_fill_manual(values=col.ab)+
  		theme_classic()+
  		scale_x_discrete(limits = select, expand = c(.1, .1))+
  		scale_y_continuous(breaks=(0:10)*10,expand=c(0,0,0.15,0))+
  		theme(legend.position="none")+xlab("")+ylab("Percentage of bins")

	#ggsave(plot,file=paste0(prefix.path,"Plot/Alluvional_all8.pdf"),width=3,height=3)

	return(plot)

}



plot.compartment.shift<-function(E.all,select.major,major.shift,col.shift,grouped=FALSE,col.grouped){

	data<-E.all
	data$label<-paste(data$DMSO,data$JQ1,sep="->")

	## Grouped by major categories of shifts
	data$class<-""
	for(i in 1:length(major.shift)){  data$class[data$label %in% major.shift[[i]]]<-names(major.shift)[i] }
	data$class<-factor(paste0(data$class),levels=names(major.shift))

	br<-c(0,1,5,10,25,50,100,250,500,1000)
	
	if(!grouped){
		plot<-ggplot(data=data[data$class %in% select.major],aes(x=(Peaks.bp.Mb+1),colour=label,fill=label) )+
			geom_density(aes(y=..scaled..),trim = TRUE)+
			scale_fill_manual(values=alpha(col.shift,0.15))+
			scale_colour_manual(values=col.shift)+
			scale_x_continuous(limits=c(1,1000),trans="log2",breaks=br+1,labels=br,expand=c(0,0.1,0,0.1))+
			facet_grid(class~.)+
			theme_classic()+xlab("H3K27ac peaks coverage on Mb")+ylab("Scaled density")+
			theme(legend.position="none",panel.grid.major.y = element_line(color = "gray",size = 0.1),
				panel.grid.major.x = element_line(color = "gray",size = 0.1))
	} else {

		plot<-ggplot(data=data[data$class %in% select.major],aes(x=(Peaks.bp.Mb+1),colour=class,fill=class) )+
			geom_density(aes(y=..scaled..),trim = TRUE)+
			scale_fill_manual(values=alpha(col.group,0.15))+
			scale_colour_manual(values=col.group)+
			scale_x_continuous(limits=c(1,1000),trans="log2",breaks=br+1,labels=br,expand=c(0,0.1,0,0.1))+
			facet_grid(class~.)+
			theme_classic()+xlab("H3K27ac peaks coverage on Mb")+ylab("Scaled density")+
			theme(legend.position="none",panel.grid.major.y = element_line(color = "gray",size = 0.1),
				panel.grid.major.x = element_line(color = "gray",size = 0.1))

	}

	return(plot)

}












