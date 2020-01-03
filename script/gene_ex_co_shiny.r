library(data.table)

plot_gene_exco <- function(mut_df,gene_col,var_col,sample_col,outprefix,top = 30,pvalue = c(0.05, 0.01),
                           returnAll = FALSE,findPathways = TRUE,kMax = 3,fontSize = 0.8,verbose = TRUE,
                           genes = NULL,width=8,height=6){

  gene_df <- mut_df[,c(gene_col,var_col,sample_col)]
  colnames(gene_df) <- c("Hugo_Symbol","VAR_TYPE","SMP_ID")
  
  if(is.null(genes)){
  #genes = getGeneSummary(x = maf)[1:top, Hugo_Symbol]
  
  gene_var <- table(gene_df$Hugo_Symbol,gene_df$VAR_TYPE)
  gene_var <- as.data.frame.matrix(gene_var)
  gene_var$Hugo_Symbol <- rownames(gene_var)
  gene_var <- data.table(gene_var)
  #head(gene_var)
  gene_var$SamplesNum <- rowSums(subset(gene_var, select =-Hugo_Symbol))
  gene_var <- plyr::arrange(gene_var,SamplesNum,decreasing = TRUE)
  genes = unique(gene_var$Hugo_Symbol)[1:top]
  }

  if(length(genes) < 2){
    stop("Minimum two genes required!")
  }

  #om = createOncoMatrix(m = maf, g = genes)
  #all.tsbs = as.character(getSampleSummary(x = maf)[,Tumor_Sample_Barcode])
  all.tsbs = unique(gene_df$SMP_ID)
  
  mutMat = mycreateOncoMatrix(m = gene_df, g = genes)
  missing.tsbs = all.tsbs[!all.tsbs %in% rownames(mutMat)]
  
  if(nrow(mutMat) < 2){
    stop("Minimum two genes required!")
  }
  
  if(length(missing.tsbs) > 0){
    missing.tsbs = as.data.frame(matrix(data = 0, nrow = length(missing.tsbs), ncol = ncol(mutMat)),
                                 row.names = missing.tsbs)
    colnames(missing.tsbs) = colnames(mutMat)
    mutMat = rbind(mutMat, missing.tsbs)
  }
  
  
  #### from maftools package, somaticInteractions
  interactions = sapply(1:ncol(mutMat), function(i) sapply(1:ncol(mutMat), function(j) {f<- try(fisher.test(mutMat[,i], mutMat[,j]), silent=TRUE); if(class(f)=="try-error") NA else ifelse(f$estimate>1, -log10(f$p.val),log10(f$p.val))} ))
  oddsRatio <- oddsGenes <- sapply(1:ncol(mutMat), function(i) sapply(1:ncol(mutMat), function(j) {f<- try(fisher.test(mutMat[,i], mutMat[,j]), silent=TRUE); if(class(f)=="try-error") f=NA else f$estimate} ))
  rownames(interactions) = colnames(interactions) = rownames(oddsRatio) = colnames(oddsRatio) = colnames(mutMat)
  
  if(returnAll){
    sigPairs = which(x = 10^-abs(interactions) < 1, arr.ind = TRUE)
  }else{
    sigPairs = which(x = 10^-abs(interactions) < max(pvalue), arr.ind = TRUE)
  }
  
  if(nrow(sigPairs) < 1){
    stop("No meaningful interactions found.")
  }
  
  sigPairsTbl = data.table::rbindlist(
    lapply(X = seq_along(1:nrow(sigPairs)), function(i) {
      x = sigPairs[i,]
      g1 = rownames(interactions[x[1], x[2], drop = FALSE])
      g2 = colnames(interactions[x[1], x[2], drop = FALSE])
      tbl = as.data.frame(table(apply(X = mutMat[,c(g1, g2), drop = FALSE], 1, paste, collapse = "")))
      combn = data.frame(t(tbl$Freq))
      colnames(combn) = tbl$Var1
      pval = 10^-abs(interactions[x[1], x[2]])
      fest = oddsRatio[x[1], x[2]]
      d = data.table::data.table(gene1 = g1,
                                 gene2 = g2,
                                 pValue = pval, oddsRatio = fest)
      d = cbind(d, combn)
      d
    }), fill = TRUE)
  
  sigPairsTbl = sigPairsTbl[!gene1 == gene2] #Remove doagonal elements
  sigPairsTbl$Event = ifelse(test = sigPairsTbl$oddsRatio > 1, yes = "Co_Occurance", no = "Mutually_Exclusive")
  sigPairsTbl$pair = apply(X = sigPairsTbl[,.(gene1, gene2)], MARGIN = 1, FUN = function(x) paste(sort(unique(x)), collapse = ", "))
  sigPairsTblSig = sigPairsTbl[order(as.numeric(pValue))][!duplicated(pair)]
  
  #Source code borrowed from: https://www.nature.com/articles/ncomms6901
  if(nrow(interactions) >= 2){
    interactions[10^-abs(interactions) > max(pvalue)] = 0
    diag(interactions) <- 0
    m <- nrow(interactions)
    n <- ncol(interactions)
    
    interactions[interactions < -4] = -4
    interactions[interactions > 4] = 4
    
    interactions[lower.tri(x = interactions)] = NA
    
    
    pdfname <- paste(outprefix,".gene.ex.co.pdf",sep="")
    #pdf(pdfname,width=width,height=height)
    
    par(bty="n", mgp = c(2,.5,0), mar = c(2, 4, 3, 5)+.1, las=2, tcl=-.33)
    image(x=1:n, y=1:m, interactions, col=RColorBrewer::brewer.pal(9,"PiYG"),
          breaks = c(-4:0-.Machine$double.eps,0:4), xaxt="n", yaxt="n",
          xlab="",ylab="", xlim=c(0, n+4), ylim=c(0, n+2))
    abline(h=0:n+.5, col="white", lwd=.5)
    abline(v=0:n+.5, col="white", lwd=.5)
    
    mtext(side = 2, at = 1:m, text = colnames(interactions), cex = fontSize, font = 3)
    mtext(side = 3, at = 1:n, text = colnames(interactions),cex = fontSize,line=-1, font = 3,las = 2)
    
    w = arrayInd(which(10^-abs(interactions) < min(pvalue)), rep(m,2))
    points(w, pch="*", col="black")
    w = arrayInd(which(10^-abs(interactions) < max(pvalue)), rep(m,2))
    points(w, pch=".", col="black")
    #image(y = 1:8 +6, x=rep(n,2)+c(2,2.5)+1, z=matrix(c(1:8), nrow=1), col=brewer.pal(8,"PiYG"), add=TRUE)
    image(y = seq(0.5*nrow(interactions), 0.9*nrow(interactions), length.out = 8), x=rep(n,2)+c(2,2.5)+1, z=matrix(c(1:8), nrow=1), col = RColorBrewer::brewer.pal(8,"PiYG"), add=TRUE)
    #axis(side = 4, at = seq(1,7) + 6.5,  tcl=-.15, label=seq(-3, 3), las=1, lwd=.5)
    atLims = seq(0.5*nrow(interactions), 0.9*nrow(interactions), length.out = 7)
    axis(side = 4, at = atLims,  tcl=-.15, labels =c(3:1, 0, 1:3), las=1, lwd=.5)
    mtext(side=4, at = median(atLims), "-log10 (p-value)", las=3, cex = 0.9, line = 2, font = 2)
    
    par(xpd=NA)
    text(x=n+2.2, y= max(atLims)+1.2, "Co-occurance", pos=4, cex = 0.9, font = 2)
    text(x=n+2.2, y = min(atLims)-1.2, "Exclusive", pos=4, cex = 0.9, font = 2)
    
    points(x = n+1, y = 0.2*n, pch = "*", cex = 2)
    text(x = n+1, y = 0.2*n, paste0(" p < ", min(pvalue)), pos=4, cex = 0.9, font = 2)
    points(x = n+1, y = 0.15*n, pch = ".", cex = 4)
    text(x = n+1, y = 0.15*n, paste0("p < ", max(pvalue)), pos=4, cex = 0.9, font = 2)
    
    #dev.off()
    }
}

###################

mycreateOncoMatrix <- function(m, g = NULL){
submatrix <- m[m$Hugo_Symbol %in% g,]
submatrix$Hugo_Symbol <- as.character(submatrix$Hugo_Symbol )
submatrix <- data.table(submatrix)

oncomat = data.table::dcast(data = submatrix[,.(Hugo_Symbol, VAR_TYPE, SMP_ID)], formula = Hugo_Symbol ~ SMP_ID,
                            fun.aggregate = function(x){
                              x = unique(as.character(x))
                              if(length(x)>0){
                                x = ifelse(test = length(x) > 1, yes = 'Multi_Hit', no = x)
                              }
                              return(x)
                            } , value.var = 'VAR_TYPE', fill = '', drop = FALSE)

#convert to matrix
data.table::setDF(oncomat)
rownames(oncomat) = oncomat$Hugo_Symbol
oncomat = as.matrix(oncomat[,-1, drop = FALSE])

variant.classes = as.character(unique(submatrix[,VAR_TYPE]))
variant.classes = c('',variant.classes, 'Multi_Hit')
names(variant.classes) = 0:(length(variant.classes)-1)


oncomat.copy <- oncomat
#Make a numeric coded matrix
for(i in 1:length(variant.classes)){
  oncomat[oncomat == variant.classes[i]] = names(variant.classes)[i]
}


#If maf has only one gene
if(nrow(oncomat) == 1){
  mdf  = t(matrix(as.numeric(oncomat)))
  rownames(mdf) = rownames(oncomat)
  colnames(mdf) = colnames(oncomat)
  return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes))
}

#convert from character to numeric
mdf = as.matrix(apply(oncomat, 2, function(x) as.numeric(as.character(x))))
rownames(mdf) = rownames(oncomat.copy)


#If MAF file contains a single sample, simple sorting is enuf.
if(ncol(mdf) == 1){
  sampleId = colnames(mdf)
  mdf = as.matrix(mdf[order(mdf, decreasing = TRUE),])
  colnames(mdf) = sampleId
  
  oncomat.copy = as.matrix(oncomat.copy[rownames(mdf),])
  colnames(oncomat.copy) = sampleId
  
  return(numericMatrix = mdf)
} else{
  #Sort by rows as well columns if >1 samples present in MAF
  #Add total variants per gene
  mdf = cbind(mdf, variants = apply(mdf, 1, function(x) {
    length(x[x != "0"])
  }))
  #Sort by total variants
  mdf = mdf[order(mdf[, ncol(mdf)], decreasing = TRUE), ]
  #colnames(mdf) = gsub(pattern = "^X", replacement = "", colnames(mdf))
  nMut = mdf[, ncol(mdf)]
  
  mdf = mdf[, -ncol(mdf)]
  
  #mdf.temp.copy = mdf #temp copy of original unsorted numeric coded matrix
  
  mdf[mdf != 0] = 1 #replacing all non-zero integers with 1 improves sorting (& grouping)
  tmdf = t(mdf) #transposematrix
  mdf = tmdf[do.call(order, c(as.list(as.data.frame(tmdf)), decreasing = TRUE)), ] #sort
  
  # mdf.temp.copy = mdf.temp.copy[rownames(mdf),] #organise original matrix into sorted matrix
  # mdf.temp.copy = mdf.temp.copy[,colnames(mdf)]
  # mdf = mdf.temp.copy
  
  #organise original character matrix into sorted matrix
  # oncomat.copy <- oncomat.copy[,colnames(mdf)]
  # oncomat.copy <- oncomat.copy[rownames(mdf),]
  
  #return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes))
  return(mdf)
  
  }
}

