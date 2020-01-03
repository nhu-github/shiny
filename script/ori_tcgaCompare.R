
tcgaCompare_ori = function(var, capture_size = NULL, tcga_capture_size = 50, cohortName = NULL, tcga_cohorts = NULL, primarySite = FALSE, col = c('gray70', 'black'), bg_col = c('#EDF8B1', '#2C7FB8'), medianCol = 'red', logscale = TRUE, rm_hyper = FALSE){
  
  tcga.cohort = system.file('extdata', 'tcga_cohort.txt.gz', package = 'maftools')
  
  if(Sys.info()[['sysname']] == 'Windows'){
    tcga.cohort.gz = gzfile(description = tcga.cohort, open = 'r')
    tcga.cohort <- suppressWarnings( data.table(read.csv( file = tcga.cohort.gz, header = TRUE, sep = '\t', stringsAsFactors = FALSE)) )
    close(tcga.cohort.gz)
  } else{
    tcga.cohort = data.table::fread(cmd = paste('zcat <', tcga.cohort), sep = '\t', stringsAsFactors = FALSE)
  }
  
  if(primarySite){
    tcga.cohort = tcga.cohort[,.(Tumor_Sample_Barcode, total, site)]
    colnames(tcga.cohort)[3] = 'cohort'
  }else{
    tcga.cohort = tcga.cohort[,.(Tumor_Sample_Barcode, total, cohort)]
  }
  
  if(!is.null(tcga_cohorts)){
    tcga.cohort = tcga.cohort[cohort %in% tcga_cohorts]
    if(nrow(tcga.cohort) == 0){
      stop("Something went wrong. Provide correct names for 'tcga_cohorts' arguments")
    }
  }
  # if(length(maf) == 1){
  #   maf = list(maf)
  # }
  # maf.mutload = lapply(maf, function(m){
  #   x = getSampleSummary(m)[,.(Tumor_Sample_Barcode, total)]
  #   x
  # })
  # change 1
  mafmy <-  unique(var[,c("ORDER_ID","TMB")])
  colnames(mafmy) <- c("Tumor_Sample_Barcode","total")
  mafmy$total <- as.numeric(mafmy$total)
  mafmy$Tumor_Sample_Barcode <- as.factor(mafmy$Tumor_Sample_Barcode)
  mafmy <- list(mafmy)
  mafmy[[1]] <- data.table(mafmy[[1]])
  maf.mutload<- mafmy

  maf <- list(var)
  if(is.null(cohortName)){
    cohortName = paste0('Input', seq_len(length(maf)))
  }else if(length(cohortName) != length(maf)){
    stop("Please provide names for all input cohorts")
  }
  
  names(maf.mutload) = cohortName
  maf.mutload = data.table::rbindlist(l = maf.mutload, idcol = "cohort")
  
  #maf.mutload[,cohort := cohortName]
  tcga.cohort$total = as.numeric(as.character(tcga.cohort$total))
  maf.mutload$total = as.numeric(as.character(maf.mutload$total))
  
  if(rm_hyper){
    #Remove outliers from tcga cohort
    tcga.cohort = split(tcga.cohort, as.factor(tcga.cohort$cohort))
    tcga.cohort = lapply(tcga.cohort, function(x){
      xout = boxplot.stats(x = x$total)$out
      if(length(xout) > 0){
        message(paste0("Removed ", length(xout), " outliers from ", x[1,cohort]))
        x = x[!total %in% xout]
      }
      x
    })
    tcga.cohort = data.table::rbindlist(l = tcga.cohort)
    #Remove outliers from input cohort
    xout = boxplot.stats(x = maf.mutload$total)$out
    if(length(xout) > 0){
      message(paste0("Removed ", length(xout), " outliers from Input MAF"))
      maf.mutload = maf.mutload[!total %in% xout]
    }
  }
  
  if(!is.null(capture_size)){
    maf.mutload[,total := total/capture_size]
    tcga.cohort[,total := total/tcga_capture_size]
  }
  
  tcga.cohort = rbind(tcga.cohort, maf.mutload)
  
  tcga.cohort.med = tcga.cohort[,.(.N, median(total)),cohort][order(V2, decreasing = TRUE)]
  
  tcga.cohort$cohort = factor(x = tcga.cohort$cohort,levels = tcga.cohort.med$cohort)
  colnames(tcga.cohort.med) = c('Cohort', 'Cohort_Size', 'Median_Mutations')
  
  tcga.cohort$TCGA = ifelse(test = tcga.cohort$cohort %in% cohortName,
                            yes = 'Input', no = 'TCGA')
  #return(tcga.cohort)
  
  tcga.cohort = split(tcga.cohort, as.factor(tcga.cohort$cohort))
  plot.dat = lapply(seq_len(length(tcga.cohort)), function(i){
    x = tcga.cohort[[i]]
    x = data.table::data.table(rev(seq(i-1, i, length.out = nrow(x))),
                               x[order(total, decreasing = T), total],
                               x[,TCGA])
    x
  })
  names(plot.dat) = names(tcga.cohort)
  if(logscale){
    y_lims = range(log10(data.table::rbindlist(l = plot.dat)[,V2]))
  }else{
    y_lims = range(data.table::rbindlist(l = plot.dat)[,V2])
  }
  
  #y_lims = range(log10(unlist(lapply(plot.dat, function(x) max(x[,V2], na.rm = TRUE)))))
  y_max = ceiling(max(y_lims))
  y_min = floor(min(y_lims))
  #change 2
  if(y_min <0){
    y_min=0
  }
  y_lims = c(y_min, y_max)
  y_at = pretty(y_lims)
  
  par(mar = c(4, 3, 2, 1))
  plot(NA, NA, xlim = c(0, length(plot.dat)), ylim = y_lims, axes = FALSE, xlab = NA, ylab = NA)
  rect(xleft = seq(0, length(plot.dat)-1, 1), ybottom = min(y_lims), xright = seq(1, length(plot.dat), 1),
       ytop = y_max, col = grDevices::adjustcolor(col = bg_col, alpha.f = 0.2),
       border = NA)
  abline(h = pretty(y_lims), lty = 2, col = "gray70")
  #abline(v = seq(1, length(plot.dat)), lty = 1, col = "gray70")
  lapply(seq_len(length(plot.dat)), function(i){
    x = plot.dat[[i]]
    if(x[1, V3] == "Input"){
      if(logscale){
        points(x$V1, log10(x$V2), pch = 16, cex = 0.4, col = col[2])
      }else{
        points(x$V1, x$V2, pch = 16, cex = 0.4, col = col[2])
      }
    }else{
      if(logscale){
        points(x$V1, log10(x$V2), pch = 16, cex = 0.4, col = col[1])
      }else{
        points(x$V1, x$V2, pch = 16, cex = 0.4, col = col[1])
      }
    }
  })
  axis(side = 2, at = y_at, las = 2, line = -1, tick = FALSE)
  samp_sizes = lapply(plot.dat, nrow)
  axis(side = 1, at = seq(0.5, length(plot.dat)-0.5, 1), labels = names(plot.dat),
       las = 2, tick = FALSE, line = -1)
  axis(side = 3, at = seq(0.5, length(plot.dat)-0.5, 1), labels = unlist(samp_sizes),
       las = 2, tick = FALSE, line = -1.2, font = 3)
  if(logscale){
    if(is.null(capture_size)){
      mtext(text = "TMB", side = 2, line = 1.2)
    }else{
      mtext(text = "TMB log10(per MB)", side = 2, line = 1.2)
    }
  }else{
    if(is.null(capture_size)){
      mtext(text = "TMB", side = 2, line = 1.2)
    }else{
      mtext(text = "TMB log10", side = 2, line = 1.6)
    }
  }
  
  
  tcga.cohort.med[, Median_Mutations_log10 := log10(Median_Mutations)]
  if(logscale){
    lapply(seq_len(nrow(tcga.cohort.med)), function(i){
      segments(x0 = i-1, x1 = i, y0 = tcga.cohort.med[i, Median_Mutations_log10],
               y1 = tcga.cohort.med[i, Median_Mutations_log10], col = medianCol)
    })
  }else{
    lapply(seq_len(nrow(tcga.cohort.med)), function(i){
      segments(x0 = i-1, x1 = i, y0 = tcga.cohort.med[i, Median_Mutations],
               y1 = tcga.cohort.med[i, Median_Mutations], col = medianCol)
    })
  }
  
  tcga.cohort.med
}

#rm(list = ls())
#dat <- openxlsx::read.xlsx("E:/song/project/20191106biomarker/shiny/demo/pipeline/example/example_mutation.xlsx")
#tcgaCompare_ori(var = dat,cohortName = 'ori',medianCol = "gray70",col = c( "black","black"))
#tcgaCompare_ori(var = dat,cohortName = 'ori')





