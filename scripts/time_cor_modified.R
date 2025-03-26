
require(rdist)




# modified script for computing time correlation
compute_time <- function(cds, bulk_expr, 
                         name = "BestTime", 
                         ref.time = NULL,
                         pseudoCountSingle = 1,
                         pseudoCountBulk = 1,
                         minAutoCor = 0.5, minSD = 1,
                         log=T, 
                         detect_low = 0, 
                         detect_high = Inf,
                         return_cds=T,
                         plotfile=NULL, plotidx = 1:50, span=.5) {
  
  if(is.null(ref.time) || length(ref.time) != ncol(bulk_expr)) {
    stop("Time reference required.")
  } 
  
  sc_expr <- exprs(cds)
  
  sharedg <- intersect(rownames(bulk_expr), rownames(sc_expr))
  message(paste0("Keep ", length(sharedg), " shared genes."))
  bulk_expr <- bulk_expr[sharedg,]
  sc_expr <- sc_expr[sharedg,]
  
  if(log){
    log_bulk = log2(bulk_expr + pseudoCountBulk)
    log_sc = log2(sc_expr + pseudoCountSingle)
  } else{
    log_bulk = tc_data_scGenes
    log_sc = sc_expr
  }
  
  #Filter keeps only genes with high variance in tc data, high autocor (suggesting temporal continuity, if appropriate) 
  #and expressed in a decent number of single cells
  tc_autocor = apply(log_bulk,MARGIN=1,function(x){cor(x,x[c(2:length(x),1)])})
  tc_sd=apply(log_bulk,MARGIN=1,sd)
  Keep = tc_sd>minSD & tc_autocor>minAutoCor 
  Time_Genes = rownames(log_bulk)[Keep]
  message(paste0(length(Time_Genes), " genes passed bulk variance and autocorrlation filtering"))
  
  log_bulk <- log_bulk[Time_Genes,]
  log_sc <- log_sc[Time_Genes,]
  
  message("finding best matches")
  
  warn_thresh <- 50 # Warn if too few genes were used for computing correlation
  warn_num <- 0
  cor_mtx <- sapply(1:ncol(log_sc), function(i){
    if(i%%1000==0){print(paste0("Computed ",i," cells"))}
    x <- log_sc[,i]
    use_gidx <- x>=detect_low & x<=detect_high
    scexpr<- x[use_gidx]
    if(length(scexpr) < warn_thresh) {
      warn_num <<- warn_num + 1
    }
    cor(scexpr, log_bulk[use_gidx,], method="pearson")
  })
  
  message(paste0(warn_num, " Cells with too few (less than ", warn_thresh, ") genes for computing correlation."))
  
  if(!is.null(plotfile)){
    graphics.off()
    pdf(plotfile)
    for(i in plotidx) {
      test_plot <- as.data.frame(cor_mtx[,i, drop=F])
      test_plot$time <- ref.time
      colnames(test_plot) <- c("Cor", "Time")
      peak_plot(x=test_plot$Time, y=test_plot$Cor, span=span, peak.def=max, colors = c("Gray", "Black", "Red"), title="Loess fit of Correlation")
    }
    dev.off()
  }
  
  bt_idx = apply(cor_mtx, MARGIN = 2, SearchPeak, searchfun = which.max, t = ref.time, span = span)
  SDs = apply(cor_mtx,MARGIN=2,sd)
  DynamicRange = apply(cor_mtx,MARGIN=2,max)-apply(cor_mtx,MARGIN=2,min)
  BestTimes = ref.time[bt_idx]
  
  if(return_cds) {
    pData(cds)[,name]<-BestTimes
    pData(cds)[,paste(name,"_SD",sep="")]<-SDs
    pData(cds)[,paste(name,"_DR",sep="")]<-DynamicRange
    return(cds)
  } else {
    return(BestTimes)
  }
}







# Function deprecated, this scheme is little bit too complicated 
# compute time based on angular distance of expressed genes between single cell and bulk under estimated single-bulk scaling ratio
compute_time_deprecated <- function(cds, bulk_expr, 
                                    name = "BestTime", 
                                    ref.time = NULL,
                                    pseudoCountSingle = 1,
                                    pseudoCountBulk = 1,
                                    minAutoCor = 0.5, minSD = 1,
                                    log=T, 
                                    detect_low = 0, # Detection threshold for a gene to be used for computing correlation for A CELL
                                    detect_high = Inf, # At high expression level, single cell may not correlate well with bulk
                                    plotfile=NULL, plotidx = 1:50,
                                    return_cds=T, span=0.5,
                                    regress_fun = lm,
                                    example_idx = NULL, ...) {
  
  if(is.null(ref.time) || length(ref.time) != ncol(bulk_expr)) {
    stop("Time reference required.")
  } 
  
  sc_expr <- exprs(cds)
  
  sharedg <- intersect(rownames(bulk_expr), rownames(sc_expr))
  message(paste0("Keep ", length(sharedg), " shared genes."))
  bulk_expr <- bulk_expr[sharedg,]
  sc_expr <- sc_expr[sharedg,]
  
  if(log){
    log_bulk = log2(bulk_expr + pseudoCountBulk)
    log_sc = log2(sc_expr + pseudoCountSingle)
  } else{
    log_bulk = tc_data_scGenes
    log_sc = sc_expr
  }
  
  #Filter keeps only genes with high variance in tc data, high autocor (suggesting temporal continuity, if appropriate) 
  #and expressed in a decent number of single cells
  tc_autocor = apply(log_bulk,MARGIN=1,function(x){cor(x,x[c(2:length(x),1)])})
  tc_sd=apply(log_bulk,MARGIN=1,sd)
  Keep = tc_sd>minSD & tc_autocor>minAutoCor 
  Time_Genes = rownames(log_bulk)[Keep]
  message(paste0(length(Time_Genes), " genes passed bulk variance and autocorrlation filtering"))
  
  log_bulk <- log_bulk[Time_Genes,]
  log_sc <- log_sc[Time_Genes,]
  
  message("Estimating single cell - bulk scaling factor...")
  
  warn_thresh <- 20 # Warn if too few genes were used for computing correlation
  warn_num <- 0
  
  fit_mtx <- sapply(1:ncol(log_sc),function(i){
    if(i %% 100 == 0) print(paste0(i," cells processed with regression fit"))
    x = log_sc[,i]
    use_gidx <- x>=detect_low & x<=detect_high
    scexpr<- x[use_gidx]
    if(length(scexpr) < warn_thresh) {
      warn_num <<- warn_num + 1
    }
    apply(log_bulk[use_gidx,], 2, function(y){
      lmfit <- regress_fun(y ~ scexpr + 0, ...)
      lmfit$coefficients[1]
    })
  })
  message(paste0(warn_num, " Cells with too few (less than ", warn_thresh, " genes for computing angular distance."))
  #BestTime = ref.time[apply(fit_mtx,MARGIN=2,SearchPeak, searchfun = which.max, t=ref.time, span=span)]
  fit_max = apply(fit_mtx,MARGIN=2,SearchPeak, searchfun = max, t=ref.time, span=span)
  
  res_mtx <- sapply(1:ncol(log_sc),function(i){
    if(i %% 100 == 0) print(paste0(i," cells computed angular distance"))
    x <- log_sc[,i]
    use_gidx <- x>=detect_low & x<=detect_high
    x<- x[use_gidx]
    fit_coef <- fit_max[i]
    yhat <- x * fit_coef
    apply(log_bulk[use_gidx,], 2, function(y) {
      sum((y-yhat)^2)
    })
    #as.numeric(cdist(t(as.matrix(yhat)), t(log_bulk[use_gidx,]), metric = "angular"))
  })
  
  if(!is.null(plotfile)){
    graphics.off()
    if(!is.null(example_idx)) {
      pdf(paste0(plotfile, "_cell", example_idx, ".pdf"), width = 7, height=7)
      x <- log_sc[,example_idx]
      use_gidx <- x>=detect_low & x<=detect_high
      scexpr<- x[use_gidx]
      for(j in 1: ncol(log_bulk)) {
        dat <- as.data.frame(cbind(y=log_bulk[use_gidx,j], x=scexpr))
        p<-ggplot(dat, aes(x=x, y=y)) + 
          geom_point()+
          geom_smooth(method=regress_fun, formula = y~x+0, se=T, fullrange=T) +
          xlim(0,min(10,detect_high)) + 
          theme_bw()
        print(p)
      }
      dev.off()
    }
    pdf(paste0(plotfile, "_peak_plot.pdf"), width=8, height=4)
    for(i in plotidx) {
      test_plot <- as.data.frame(res_mtx[,i, drop=F])
      test_plot$time <- ref.time
      test_plot$FitCoef <- fit_mtx[,i]
      colnames(test_plot) <- c("Res", "Time", "FitCoef")
      par(mfrow=c(1,2))
      peak_plot(x=test_plot$Time, y=test_plot$FitCoef, span=span, peak.def=max, colors = c("Gray", "Black", "Blue"), title="Loess fit of LM Coefficient")
      peak_plot(x=test_plot$Time, y=test_plot$Res, span=span, peak.def=min, colors = c("Gray", "Black", "Orange"), title="Loess fit of Residual")
    }
    dev.off()
  }
  
  res_idx = apply(res_mtx,MARGIN=2,SearchPeak, searchfun = which.min, t=ref.time, span=span)
  
  BestTime = ref.time[res_idx]
  
  if(return_cds) {
    pData(cds)[,paste0(name)]<-BestTime
    return(cds)
  } else {
    return(
      BestTime
    )
  }
}



peak_plot <- function(x,y,span, peak.def = max, colors = c("Grey", "Black", "Red"), title=NULL) {
  y.smooth <- loess(y ~ x, span=span)$fitted
  idx <- which(y.smooth == peak.def(y.smooth))
  plot(x, y, cex=0.75, col=colors[1], main=title)
  lines(x, y.smooth,  lwd=2, col=colors[2]) #$
  y.min <- min(y)
  lines(c(x[idx],x[idx]), c(y.min, y.smooth[idx]),col=colors[3], lty=2)
  points(x[idx],  y.smooth[idx], col=colors[3], pch=19, cex=1.25)
}



SearchPeak <- function(x, t, span, searchfun = which.max) {
  # cat("Length of x:", length(x), " Length of t:", length(t), "\n")
  y.smooth <- loess(x ~ t, span=span)$fitted
  idx <- searchfun(y.smooth)
  if(length(idx) > 1){
    idx <- idx[1]
  }
  return(idx)
}


















