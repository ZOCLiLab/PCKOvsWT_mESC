require(R.filesets)
require(Seurat) #v3.1.1
require(Matrix)
require(foreach)
require(doParallel)
require(future)
require(dplyr)
require(data_mat)

`%+%` <- paste0

workDir = "./PCKOvsWT_mESCs/"
setwd(workDir)

options(stringsAsFactors=F)
options(future.globals.maxSize= 2048 * 1024^2)

create_dir = function(adirpath){
  if (!dir.exists(adirpath)){
    dir.create(adirpath, recursive = TRUE)
  }
}

load_10X_data = function(adirpath){
  tmp_files = list.files(adirpath)
  if (length(which(tmp_files %in% c("barcodes.tsv.gz","features.tsv.gz","matrix.mtx.gz")))==3){
    tmp_mm = Read10X(adirpath)
    return(tmp_mm)
  } else {
    print(adirpath %+% " contains no 10X counts.")
    return(NULL)
  }
}

pretreat_cellFiltration = function(raw_expr_mat, report=TRUE, report_dir="."){
  ## filter potential cells of library with simple components/poor library/doublets
  cS = data.frame(libSize=colSums(raw_expr_mat), geneDetect=apply(raw_expr_mat, 2, function(X) sum(X>0)))
  tB = with(cS, geneDetect < (libSize*0.1))
  
  p_hi <- 1e-3
  p_lo <- 1e-2
  ## library size threshold
  tempFitLS <- fitdistr(cS$libSize[!tB], "negative binomial")
  c_hiL <- qnbinom(p_hi, size=tempFitLS$estimate["size"], mu=tempFitLS$estimate["mu"], lower.tail=F)
  c_loL <- qnbinom(p_lo, size=tempFitLS$estimate["size"], mu=tempFitLS$estimate["mu"], lower.tail=T)
  ## gene detection threshold
  tempFitGD <- fitdistr(cS$geneDetect[!tB], "negative binomial")
  c_hiG <- qnbinom(p_hi, size=tempFitGD$estimate["size"], mu=tempFitGD$estimate["mu"], lower.tail=F)
  c_loG <- qnbinom(p_lo, size=tempFitGD$estimate["size"], mu=tempFitGD$estimate["mu"], lower.tail=T)
  ## marker dangerous cells
  temp_doublets <- (cS$libSize[!tB] > c_hiL) | (cS$geneDetect[!tB] > c_hiG)
  temp_crapLibs <- (cS$libSize[!tB] < c_loL) | (cS$geneDetect[!tB] < c_loG)
  ## count matrix filtration
  raw_expr_mat_f <- raw_expr_mat[,!tB]
  raw_expr_mat_f <- raw_expr_mat_f[,!(temp_doublets | temp_crapLibs)]
  raw_expr_mat_f <- raw_expr_mat_f[rowSums(raw_expr_mat_f)>0,]
  ## report QC results
  if (report){
    results = list(cS=cS, tB=tB, rmat=raw_expr_mat, rmat_f=raw_expr_mat_f, temp_doublets=temp_doublets, temp_crapLibs=temp_crapLibs)
    pcontaminated_glist = report_qc_figures(results=results, output_dir=report_dir, check_contaminated=TRUE)
    return(list(expr_mat = raw_expr_mat_f, pc_glist=pcontaminated_glist))
  }
  return(raw_expr_mat_f)
}

pretreat_mitoFiltration = function(raw_expr_mat, report=TRUE, report_dir=".", mito_glist=NULL, ratio=0.4){
  if (is.null(mito_glist)){
    mito_gind = grepl("^MT-",rownames(raw_expr_mat))
  } else {
    mito_gind = which(rownames(raw_expr_mat) %in% mito_glist )
  }
  cSf <- data.frame(libSize=colSums(raw_expr_mat),
                    geneDetect=apply(raw_expr_mat,2,function(X) sum(X>0)),
                    mitoPct=colSums(raw_expr_mat[mito_gind,])/colSums(raw_expr_mat))
  
  drop_mitoMads <- 4
  drop_mitoCut <- median(cSf$mitoPct) + mad(cSf$mitoPct)*drop_mitoMads ## quantile
  if (drop_mitoCut > ratio) { drop_mitoCut <- ratio}
  drop_mito <- cSf$mitoPct > drop_mitoCut

  raw_expr_mat_f <- raw_expr_mat[,!drop_mito]
  raw_expr_mat_f <- raw_expr_mat_f[apply(raw_expr_mat_f,1,sum) > 0,]
  
  if(report){
    results = list(cSf=cSf, drop_mito=drop_mito, 
                   drop_mitoMads=drop_mitoMads, 
                   drop_mitoCut=drop_mitoCut, 
                   rmat=raw_expr_mat, 
                   rmat_f=raw_expr_mat_f)
    report_mito_figure(results, output_dir=report_dir)
  }
  return(raw_expr_mat_f)
}

report_qc_figures = function(results, output_dir, check_contaminated=FALSE){
  
  cS = results$cS
  tB = results$tB
  rmat = results$rmat
  rmat_f = results$rmat_f
  temp_doublets = results$temp_doublets
  temp_crapLibs = results$temp_crapLibs
  
  ## cells filtration
  QC_filtration_file = paste(output_dir,"QualityControl","QC_cellsFiltration.pdf",sep="/")
  create_dir(dirname(QC_filtration_file))
  pdf(file = QC_filtration_file, width = 8.4, height = 4.2)
  layout(cbind(matrix(c(2,1,0,3),2),matrix(c(5,4,0,6),2)),widths=c(3.5,.7,3.5,.7),heights=c(.7,3.5))
  ## figure1-library per cell vs detected genes per cell
  par(mar=c(3,3,0,0),mgp=2:0)
  plot(geneDetect~libSize,data=cS[!tB,],
       xlim=range(cS$libSize),ylim=range(cS$geneDetect),
       pch=21,col=alpha("black",0.2),bg=alpha("black",0.1),cex=1.2,
       xlab="Library Size",ylab="Genes Detected")
  points(geneDetect~libSize,data=cS[tB,],
         pch=21,col=alpha("red",0.2),bg=alpha("red",0.1),cex=1.2)
  points(geneDetect~libSize,data=cS[!tB,][temp_doublets | temp_crapLibs,],
         pch=4,cex=1.2,col="red")
  mtext(paste("cell stats"),side=3,line=-1.5,font=2,cex=1.2)
  legend("topleft",bty="n",inset=c(-.02,.05),
         legend=c(paste("Total genes:",nrow(rmat)),
                  paste("Post-filter:",nrow(rmat_f)),
                  paste("Total cells:",ncol(rmat))))
  legend("bottomright",bty="n",
         pch=c(4,4,21,21),
         col=c("red","red","red","black"),
         pt.bg=alpha(c(NA,NA,"red","black"),0.3),
         y.intersp=1.3, 
         legend=c(paste(sep="\n","Poor-quality libraries",
                  paste0("(p<",p_lo,"): ",sum(temp_crapLibs))),
                  paste(sep="\n","Predicted doublets",
                  paste0("(p<",p_hi,"): ",sum(temp_doublets))),
                  paste("Outlier cells:",sum(tB)),
                  paste("Remain cells:",ncol(rmat_f))))
  ## distribution for library per cell
  par(mar=c(0,3,1,0))
  tempD <- density(rnbinom(10000,size=tempFitLS$estimate["size"],mu=tempFitLS$estimate["mu"]))
  hist(cS$libSize,breaks=100,freq=F,col="grey",main=NULL,xaxt="n",ylab="Density")
  lines(tempD,lwd=2,col=alpha("red",0.5))
  abline(v=c_hiL,lty=2,lwd=2,col="darkred")
  abline(v=c_loL,lty=2,lwd=2,col="darkred")
  ## distribution for detected genes per cell
  par(mar=c(3,0,0,1))
  tempD <- density(rnbinom(10000,size=tempFitGD$estimate["size"],mu=tempFitGD$estimate["mu"]))
  tempH <- hist(cS$geneDetect,breaks=100,plot=F)
  tempB <- barplot(tempH$density,horiz=T,space=0,col="grey",main=NULL,xlab="Density")
  tempSF <- (max(tempB) - min(tempB)) / (max(tempH$mids) - min(tempH$mids))
  lines(y=tempD$x * tempSF + (min(tempB) - min(tempH$mids) * tempSF),
        x=tempD$y,lwd=2,col=alpha("red",0.5))
  abline(h=c_hiG * tempSF + (min(tempB) - min(tempH$mids) * tempSF),lty=2,lwd=2,col="darkred")
  abline(h=c_loG * tempSF + (min(tempB) - min(tempH$mids) * tempSF),lty=2,lwd=2,col="darkred")
  ## figure2-library per cell(log scale) vs detected genes per cell 
  par(mar=c(3,3,0,0),mgp=2:0)
  plot(geneDetect~libSize,data=cS[!tB,],log="x",
       xlim=range(cS$libSize),ylim=range(cS$geneDetect),
       pch=21,col=alpha("black",0.2),bg=alpha("black",0.1),cex=1.2,
       xlab="Library Size (log scale)",ylab="Genes Detected")
  points(geneDetect~libSize,data=cS[tB,],
         pch=21,col=alpha("red",0.2),bg=alpha("red",0.1),cex=1.2)
  points(geneDetect~libSize,data=cS[!tB,][temp_doublets | temp_crapLibs,],
         pch=4,cex=1.2,col="red")
  mtext(paste("cell stats"),side=3,line=-1.5,font=2,cex=1.2)
  legend("topleft",bty="n",y.intersp=1.3,inset=c(0,.05),
         pch=c(4,4,21),col="red",pt.bg=alpha(c(NA,NA,"red"),0.3),
         legend=c(paste(sep="\n","Poor-quality libraries",
                  paste0("(p<",p_lo,"): ",sum(temp_crapLibs))),
                  paste(sep="\n","Predicted doublets",
                  paste0("(p<",p_hi,"): ",sum(temp_doublets))),
                  paste("Outlier cells:",sum(tB))))
  ## library per cell (log scale), histogram plot
  par(mar=c(0,3,1,0))
  hist(log10(cS$libSize),breaks=100,freq=F,col="grey",main=NULL,xaxt="n",ylab="Density")
  abline(v=log10(c_hiL),lty=2,lwd=2,col="darkred")
  abline(v=log10(c_loL),lty=2,lwd=2,col="darkred")
  ## detected genes per cell, histogram plot
  par(mar=c(3,0,0,1))
  tempH <- hist(cS$geneDetect,breaks=100,plot=F)
  tempB <- barplot(tempH$density,horiz=T,space=0,col="grey",main=NULL,xlab="Density")
  tempSF <- (max(tempB) - min(tempB)) / (max(tempH$mids) - min(tempH$mids))
  abline(h=c_hiG * tempSF + (min(tempB) - min(tempH$mids) * tempSF),lty=2,lwd=2,col="darkred")
  abline(h=c_loG * tempSF + (min(tempB) - min(tempH$mids) * tempSF),lty=2,lwd=2,col="darkred")
  dev.off()

  if (check_contaminated){
    ### pretreatment_contamFilt
    ## DR: detected ratio, MDTC: mean detected tx per cell, MTC: mean tx per cell
    gsOut <- data.frame(DR=apply(rmat[,tB],1,function(X) sum(X > 0)/length(X)),
                        MDTC=apply(rmat[,tB],1,function(X) mean(X[X>0])),
                        MTC=apply(rmat[,tB],1,mean))
    rownames(gsOut) <- rownames(rmat)
    gsIn <- data.frame(DR=apply(rmat_f,1,function(X) sum(X > 0)/length(X)),
                       MDTC=apply(rmat_f,1,function(X) mean(X[X>0])),
                       MTC=apply(rmat_f,1,mean))
    rownames(gsIn) <- rownames(rmat_f)
    ## contamination filtration
    QC_contamination_file = paste(output_dir,"QualityControl","QC_contaminationFiltration.pdf",sep="/")
    create_dir(dirname(QC_contamination_file))
    pdf(file = QC_contamination_file, width = 8.4, height = 4.2)
    ## figure1-outlier cells
    par(mfrow=c(1,2),mar=c(3,3,2,1),mgp=2:0)
    tempMGC <- 1e2
    plot(log10(MDTC)~DR,data=gsOut,
         pch=21,cex=1.2,xlab="Proportion of cells detecting gene",
         ylab=expression(Log[10]~"Mean non-zero gene count"),
         col=alpha(c("black","red"),0.3)[(gsOut$MTC > tempMGC)+1],
         bg=alpha(c("black","red"),0.1)[(gsOut$MTC > tempMGC)+1],
         main=paste("outlier cells"))
    text(log10(MDTC)~DR,data=gsOut[gsOut$MTC > tempMGC,],
         labels=rownames(gsOut)[gsOut$MTC > tempMGC],pos=2,col="darkred",cex=0.6)
    tempMGCnames <- rownames(gsOut[gsOut$MTC > tempMGC,])
    ## figure2-gene distribution of high expressed genes detected in outlier cells in hi-QC cells 
    plot(log10(MDTC)~DR,data=gsIn[!rownames(gsIn) %in% tempMGCnames,],
         pch=21,cex=1.2,col=alpha("black",0.3),bg=alpha("black",0.1),
         xlab="Proportion of cells detecting gene",
         ylab=expression(Log[10]~"Mean non-zero gene count"),
         main=paste("main cell population"))
    points(log10(MDTC)~DR,data=gsIn[tempMGCnames,],
           pch=21,cex=1.2,col=alpha("red",0.5),bg=alpha("red",0.3))
    text(log10(MDTC)~DR,data=gsIn[tempMGCnames,],pos=4,col="red",cex=1.2,
         labels=rownames(gsIn)[rownames(gsIn) %in% tempMGCnames])
    dev.off()
    ## potential contaminated gene list
    QC_contaminationGene_file = paste(output_dir, "QualityControl", "potential_contamination_geneList.txt", sep="/")
    write.table(tempMGCnames, file=QC_contaminationGene_file, sep="\t", col.names=TRUE, row.names=TRUE)
    return(tempMGCnames)
  }
}

report_mito_figure = function(results, output_dir="."){
  
  cSf = results$cSf
  drop_mito = results$drop_mito
  drop_mitoMads = results$drop_mitoMads
  drop_mitoCut = results$drop_mitoCut
  rmat = results$rmat
  rmat_f = results$rmat_f
  
  QC_mitochondrion_file = paste(output_dir,"QualityControl","QC_mitoFiltration.pdf",sep="/")
  create_dir(dirname(QC_mitochondrion_file))
  pdf(file = QC_mitochondrion_file, width = 9.4, height = 4.7)
  layout(matrix(c(1,2,4,3,0,5),2),c(4.1,4.1,1),c(0.7,3.8))
  par(mar=c(0,3,3,1),mgp=2:0)
  plot.new()
  title(main=paste("Filtered on mitochondrial\ngene expression proportion"))
  ## figure1-highlight cells with high mitochondrial gene ratio
  par(mar=c(3,3,0,1))
  plot(geneDetect~libSize,data=cSf,log="x",
       pch=21,cex=1.2,col=alpha("black",0.2),bg=alpha("black",0.1),
       xlab="Library Size (log scale)",ylab="Genes detected",
       xlim=range(cSf$libSize),ylim=range(cSf$geneDetect))
  points(geneDetect~libSize,data=cSf[drop_mito,],pch=4,cex=1.2,col=alpha("red",.5))
  legend("topleft",bty="n",pch=c(4,NA,NA,NA,NA),col=alpha(c("red",NA,NA,NA,NA),0.5),
         legend=c(paste(drop_mitoMads,"MADs above median"),
                  paste(sum(drop_mito),"cells removed"),
                  paste(ncol(rmat_f),"cells remain"),
                  paste(nrow(rmat)-nrow(rmat_f),"genes removed"),
                  paste(nrow(rmat_f),"genes remain")))
  ## figure2-highlight cells with high mitochondrial reads
  par(mar=c(3,3,0,0))
  plot(mitoPct~libSize,data=cSf,log="x",
       pch=21,cex=1.2,col=alpha("black",0.2),bg=alpha("black",0.1),
       xlab="Library Size (log scale)",ylab="Mitochondrial Transcript Proportion")
  abline(h=drop_mitoCut,lwd=2,lty=2,col=alpha("red",0.5))
  legend("topright",bty="n",lty=c(2,NA,NA,NA,NA),lwd=c(2,NA,NA,NA,NA),col=alpha(c("red",NA,NA,NA,NA),0.5),
         legend=c(paste(drop_mitoMads,"MADs above median"),
                  paste(sum(drop_mito),"cells removed"),
                  paste(ncol(rmat_f),"cells remain"),
                  paste(nrow(rmat)-nrow(rmat_f),"genes removed"),
                  paste(nrow(rmat_f),"genes remain")))
  ## figure3-distribution of library
  par(mar=c(0,3,0,0))
  hist(cSf$libSize,freq=T,breaks=50,col="grey",main=NULL,xaxt="n")
  ## figure4-distribution of mitochondrial proportion
  par(mar=c(3,0,0,0))
  barplot(hist(cSf$mitoPct,breaks=50,plot=F)$counts,
          horiz=T,space=0,col="grey",main=NULL,xlab="Frequency")
  dev.off()
}
