library(data.table)
library(ggplot2)

drive_gene <- function(mut_df,gene_col,sample_col,AA_change_col,Var_class_col,outprefix,var_type_col=NULL,
                       nBgGenes = 100,pvalMethod = 'zscore',minMut=3,pvalue=0.05,fdr=TRUE,width=8,height=6,WriteTable=T,plotpdf=T){
  
  if(!is.null(var_type_col)){
    prot_data <- mut_df[,c(gene_col,sample_col,AA_change_col,Var_class_col,var_type_col)]
    
    
    
    colnames(prot_data) <- c("Hugo_Symbol","SMP_ID","AA_CHANGE","Variant_Classification","var_type")
    prot_data <- prot_data[pro_data$var_type == "SNV",]
  }else{
    prot_data <- mut_df[,c(gene_col,sample_col,AA_change_col,Var_class_col)]
    colnames(prot_data) <- c("Hugo_Symbol","SMP_ID","AA_CHANGE","Variant_Classification")
  }
  
  
  prot.tmp = strsplit(x = as.character(prot_data$AA_CHANGE), split = '.', fixed = TRUE)
  prot.tmp2 = sapply(sapply(prot.tmp, function(x) x[length(x)]), '[', 1)
  
  prot_data$conv <- prot.tmp2
  
  #If conversions are in HGVSp_long (default HGVSp) format, we will remove strings Ter followed by anything (e.g; p.Asn1986GlnfsTer13)
  pos = gsub(pattern = 'Ter.*', replacement = '',x = prot_data$conv)
  
  #Following parsing takes care of most of HGVSp_short and HGVSp_long format
  pos = gsub(pattern = '[[:alpha:]]', replacement = '', x = pos)
  pos = gsub(pattern = '\\*$', replacement = '', x = pos) #Remove * if nonsense mutation ends with *
  pos = gsub(pattern = '^\\*', replacement = '', x = pos) #Remove * if nonsense mutation starts with *
  pos = gsub(pattern = '\\*.*', replacement = '', x = pos) #Remove * followed by position e.g, p.C229Lfs*18
  
  pos = suppressWarnings( as.numeric(sapply(strsplit(x = pos, split = '_', fixed = TRUE), '[', 1)) )
  
  prot_data$pos <- pos
  prot_data = prot_data[!is.na(pos),] #Remove NA's
  
  judge_syn <- strsplit(x = as.character(prot_data$conv), split = prot_data$pos, fixed = TRUE)
  
  aa1 <- sapply(judge_syn, '[', 1)
  aa2 <- sapply(judge_syn, '[', 2)
  prot_data$aa1 <- aa1
  prot_data$aa2 <- aa2
  
  prot_data$syn <- ifelse((prot_data$aa1 == prot_data$aa2) ,"yes","no")
  prot_data = prot_data[!is.na(prot_data$syn),]
  syn_data <- prot_data[prot_data$syn == "yes",]
  
  #########################################################################
  #### remove Splice_Site
  #all.prot.dat = prot_data[Variant_Classification != 'Splice_Site']
  ##########################################################################
  
  #number of samples in maf
  numSamples = length(unique(prot_data$SMP_ID))
  
 # gl <- read.table("D:/work/scientific_research_platform/gene_length.txt",sep="\t",header=1)
  gl <- read.table("./data/gene_length.txt",sep="\t",header=1)
  
  ####syn variants for background
  if(nrow(syn_data)==0){
    message('No syn mutations found! Skipping background estimation. Using predefined values. (Mean = 0.279; SD = 0.13)')
    bg.mean = 0.279
    bg.sd = 0.13
  }else if(nrow(syn_data) < nBgGenes){
    message("Not enough genes to build background. Using predefined values. (Mean = 0.279; SD = 0.13)")
    bg.mean = 0.279
    bg.sd = 0.13  
  }else{
    syn.bg.scores = parse_prot(dat = syn_data,  gl, m = minMut, calBg = TRUE, nBg = nBgGenes)
    if(is.null(syn.bg.scores)){
      message("Not enough genes to build background. Using predefined values. (Mean = 0.279; SD = 0.13)")
      bg.mean = 0.279
      bg.sd = 0.13
    }else{
      bg.mean = mean(syn.bg.scores$clusterScores)
      bg.sd = sd(syn.bg.scores$clusterScores)
      message(paste('Estimated background mean: ', bg.mean))
      message(paste('Estimated background SD: ', bg.sd))
    }
  }
  
  ## nonsys 
  # non.syn.maf = maf@data
  # nonsyn.scores = parse_prot(dat = non.syn.maf, AACol = AACol, gl = gl, m = minMut, calBg = FALSE, nBg = nBgGenes)
  non.syn.data <- prot_data[!rownames(prot_data) %in% rownames(syn_data),]
  
  nonsyn.scores = parse_prot(dat = non.syn.data, gl = gl, m = minMut, calBg = FALSE, nBg = nBgGenes)
  
  ## calculate pvalue
  if(!is.null(nonsyn.scores)){
    #gene_summary = getGeneSummary(maf)
    gene_summary <- as.data.frame.matrix(table(prot_data$Hugo_Symbol,prot_data$Variant_Classification))
    gene_summary$total <- rowSums(gene_summary)
    gene_summary$Hugo_Symbol <- rownames(gene_summary)
    if(pvalMethod == 'combined'){
      message('Comapring with background model and estimating p-values..')
      nonsyn.scores$zscore = (nonsyn.scores$clusterScores - bg.mean) / bg.sd
      nonsyn.scores$tPval = 1- pnorm(nonsyn.scores$zscore)
      nonsyn.scores$tFdr = p.adjust(nonsyn.scores$tPval, method = 'fdr')
      
      nonsyn.scores = merge(gene_summary, nonsyn.scores, by = 'Hugo_Symbol')
      nonsyn.scores <- as.data.table(nonsyn.scores)
      nonsyn.scores[,fract_muts_in_clusters := muts_in_clusters/total]
      
      counts.glm = glm(formula = total ~ protLen+clusters, family = poisson(link = identity), data = nonsyn.scores) #Poisson model
      nonsyn.scores$Expected = counts.glm$fitted.values #Get expected number of events (mutations) from the model
      
      observed_mut_colIndex = which(colnames(nonsyn.scores) == 'total')
      expected_mut_colIndex = which(colnames(nonsyn.scores) == 'Expected')
      
      #Poisson test to caluclate difference (p-value)
      nonsyn.scores$poissonPval = apply(nonsyn.scores, 1, function(x) {
        poisson.test(as.numeric(x[observed_mut_colIndex]), as.numeric(x[expected_mut_colIndex]))$p.value
      })
      
      nonsyn.scores$poissonFdr = p.adjust(nonsyn.scores$poissonPval, method = 'fdr')
      nonsyn.scores = nonsyn.scores[order(poissonFdr)]
      
      nonsyn.scores$fdr = apply(nonsyn.scores[,.(tFdr, poissonFdr)], MARGIN = 1, FUN = min)
    } else if(pvalMethod == 'zscore'){
      #Oncodrive clust way of caluclating pvalues
      #Calculate z scores; compare it to bg scores and estimate z-score, pvalues, corrected pvalues (fdr) (assumes normal distribution)
      message('Comapring with background model and estimating p-values..')
      nonsyn.scores$zscore = (nonsyn.scores$clusterScores - bg.mean) / bg.sd
      nonsyn.scores$pval = 1- pnorm(nonsyn.scores$zscore)
      nonsyn.scores$fdr = p.adjust(nonsyn.scores$pval, method = 'fdr')
      
      nonsyn.scores = merge(gene_summary, nonsyn.scores, by = 'Hugo_Symbol')
      nonsyn.scores <- as.data.table(nonsyn.scores)
      nonsyn.scores[,fract_muts_in_clusters := muts_in_clusters/total]
      #nonsyn.scores[,fract_MutatedSamples := MutatedSamples/numSamples]
      nonsyn.scores = nonsyn.scores[order(fdr)]
    }else{
      #Assuming poisson distribution of mutation counts
      #Now model observed number of mutations as a function of number of clusters and protein length. Calculate expected number of events based on poisson distribution.
      nonsyn.scores = merge(gene_summary, nonsyn.scores, by = 'Hugo_Symbol')
      nonsyn.scores <- as.data.table(nonsyn.scores)
      nonsyn.scores[,fract_muts_in_clusters := muts_in_clusters/total]
      
      counts.glm = glm(formula = total ~ protLen+clusters, family = poisson(link = identity), data = nonsyn.scores) #Poisson model
      nonsyn.scores$Expected = counts.glm$fitted.values #Get expected number of events (mutations) from the model
      
      observed_mut_colIndex = which(colnames(nonsyn.scores) == 'total')
      expected_mut_colIndex = which(colnames(nonsyn.scores) == 'Expected')
      
      #Poisson test to caluclate difference (p-value)
      nonsyn.scores$pval = apply(nonsyn.scores, 1, function(x) {
        poisson.test(as.numeric(x[observed_mut_colIndex]), as.numeric(x[expected_mut_colIndex]))$p.value
      })
      
      nonsyn.scores$fdr = p.adjust(nonsyn.scores$pval, method = 'fdr')
      nonsyn.scores = nonsyn.scores[order(fdr)]
    }
    
    ## plot
    
    theme <- theme(
      #axis.title.y=element_blank(),
      plot.title=element_text(hjust=0.5),
      panel.background = element_blank(),
      #axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      panel.border=element_rect(colour='grey60', fill=NA, size=1),
      #panel.grid.major = element_line(colour="grey90"), 
      #panel.grid.minor =  element_line(colour="grey80"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      #legend.position="none"
    )
    
    nonsyn.scores$logfdr<- -log10(nonsyn.scores$fdr)
    if(fdr){
      nonsyn.scores <- nonsyn.scores[nonsyn.scores$fdr < pvalue,]
    }else{
      nonsyn.scores <- nonsyn.scores[nonsyn.scores$pval < pvalue,] 
    }
    if(plotpdf==T){
      pdfname <- paste(outprefix,".drivegene.pdf",sep="")
      #pdf(pdfname,width=width,height=height)
      p = ggplot(nonsyn.scores,aes_string(x="Hugo_Symbol",y="logfdr"))
      p=p + geom_point(aes(size=clusterScores,color=logfdr)) +
        scale_color_gradient(low="green",high = "red")+ 
        labs(color=expression(-log[10](fdr)),y="-log10(fdr)",x="Gene") + 
        theme+geom_text(aes(Hugo_Symbol, logfdr+0.25, label=Hugo_Symbol))
      print(p)
      #dev.off()
    }

    if(WriteTable==T){
      data_table <- nonsyn.scores[,c("Hugo_Symbol","clusters","muts_in_clusters", "clusterScores", "protLen","zscore","pval","fdr")]
      filename <- paste(outprefix,"_drivegene.xlsx",sep="")
      #openxlsx::write.xlsx(data_table,file=filename,quote = FALSE, sep = "\t",row.names = FALSE)
      #print(data_table)
      return(data_table)
    }
    # data_table <- nonsyn.scores[,c("Hugo_Symbol","clusters","muts_in_clusters", "clusterScores", "protLen","zscore","pval","fdr")]
    # filename <- paste(outprefix,".drivegene.xls",sep="")
    # write.table(data_table,file=filename,quote = FALSE, sep = "\t",row.names = FALSE)
  }else{
    if(WriteTable==T){
      data_table <- c("no drive gene")
      #openxlsx::write.xlsx(data_table,file=filename,quote = FALSE, sep = "\t",row.names = FALSE)
      #print(data_table)
      return(data_table)
    }
    #data_table <- c("no drive gene")
    #write.table(data_table,file=filename,quote = FALSE, sep = "\t",row.names = FALSE,col.names =FALSE)
  }
  
}




###function

parse_prot <- function(dat,  gl, m, calBg = FALSE, nBg){
  all.prot.dat = dat[,c("Hugo_Symbol", "Variant_Classification", "AA_CHANGE")]
  #all.prot.dat = all.prot.dat[Variant_Classification != 'Splice_Site']
  #parse AAchanges to get postion
  prot.spl = strsplit(x = as.character(all.prot.dat$AA_CHANGE), split = '.', fixed = TRUE)
  prot.conv = sapply(sapply(prot.spl, function(x) x[length(x)]), '[', 1)
  
  all.prot.dat$conv <- prot.conv
  all.prot.dat = all.prot.dat[!all.prot.dat$conv == 'NULL',]
  
  #If conversions are in HGVSp_long (default HGVSp) format, we will remove strings Ter followed by anything (e.g; p.Asn1986GlnfsTer13)
  pos = gsub(pattern = 'Ter.*', replacement = '',x = all.prot.dat$conv)
  
  #Following parsing takes care of most of HGVSp_short and HGVSp_long format
  pos = gsub(pattern = '[[:alpha:]]', replacement = '', x = pos)
  pos = gsub(pattern = '\\*$', replacement = '', x = pos) #Remove * if nonsense mutation ends with *
  pos = gsub(pattern = '^\\*', replacement = '', x = pos) #Remove * if nonsense mutation starts with *
  pos = gsub(pattern = '\\*.*', replacement = '', x = pos) #Remove * followed by position e.g, p.C229Lfs*18
  
  pos = suppressWarnings( as.numeric(sapply(strsplit(x = pos, split = '_', fixed = TRUE), '[', 1)) )
  all.prot.dat$pos <- pos
  
  all.prot.dat = all.prot.dat[!is.na(pos),] #Remove NA's
  
  #gene.sum = summarizeMaf(maf = dat, chatty = FALSE)$gene.summary
  gene.sum = as.data.frame.matrix(table(dat$Hugo_Symbol,dat$Variant_Classification))
  gene.sum$total <- rowSums(gene.sum)
  gene.sum$Hugo_Symbol <- rownames(gene.sum)
  gene.sum = merge(x = gene.sum, y = gl, by="Hugo_Symbol", all.x = TRUE)
  gene.sum = gene.sum[!is.na(gene.sum$aa.length),]
  
  num_mut_colIndex = which(colnames(gene.sum) == 'total')
  aalen_colIndex = which(colnames(gene.sum) == 'aa.length')
  
  #Get background threshold
  gene.sum$th = apply(gene.sum, 1, function(x) get_threshold(gene_muts = as.numeric(x[num_mut_colIndex]), gene_length = as.numeric(x[aalen_colIndex])))
  #use only genes with atleast 2 (or m ) mutations.
  gene.sum = gene.sum[gene.sum$total >= m,]
  if(calBg){
    if(nrow(gene.sum) < nBg){
      #message("Not enough genes to build background. Using predefined values. (Mean = 0.279; SD = 0.13)")
      return(NULL)
    } else{
      syn.res = c()
      pb <- txtProgressBar(min = 0, max = nrow(gene.sum), style = 3) #progress bar
      
      for(i in 1:nrow(gene.sum)){
        prot.dat = all.prot.dat[all.prot.dat$gene %in% gene.sum[i, "gene"],]
        syn.res = rbind(syn.res, cluster_prot(prot.dat = prot.dat, gene = gene.sum[i, "gene"], th = gene.sum[i, "th"], protLen = gene.sum[i,"aa.length"]))
        setTxtProgressBar(pb, i)
      }
      return(syn.res)
    }
  } else{
    nonsyn.res = c()
    pb <- txtProgressBar(min = 0, max = nrow(gene.sum), style = 3) #progress bar
    
    for(i in 1:nrow(gene.sum)){
      hs = gene.sum[i, "Hugo_Symbol"]
      #print(hs)
      prot.dat = all.prot.dat[all.prot.dat$Hugo_Symbol %in% hs,]
      nonsyn.res = rbind(nonsyn.res, cluster_prot(prot.dat = prot.dat, Hugo_Symbol = hs, th = gene.sum[gene.sum$Hugo_Symbol %in% hs, "th"], protLen = gene.sum[gene.sum$Hugo_Symbol %in% hs, "aa.length"]))
      setTxtProgressBar(pb, i)
    }
    return(nonsyn.res)
  }
}

## get_threshold
get_threshold = function(gene_muts, gene_length){
  th = which(unlist(lapply(X = 2:gene_muts, FUN = function(x) dbinom(x = x, size = gene_muts, prob = 1/gene_length) )) < 0.01)[1]
  return(th+1)
}


## cluster_prot

cluster_prot = function(prot.dat, Hugo_Symbol, th, protLen,mergeDist = 5){
  
  mergeDist = mergeDist #hard coded inter event distance.
  
  #Summarise counts per position
  prot.dat <- as.data.table(prot.dat)
  pos.counts = prot.dat[,.N,pos]
  pos.counts = pos.counts[order(pos)]
  
  #classify position as meaningful if its greater than background threshhold.
  pos.counts$cluster = ifelse(test = pos.counts$N >= th, yes = 'meaningful', no = 'nonMeaningful')
  
  #Just choose meaningful positions
  clust.tbl = pos.counts[cluster %in% 'meaningful']
  nonclust.tbl = pos.counts[cluster %in% 'nonMeaningful']
  
  if(nrow(clust.tbl) == 0){
    #message(paste('No meaningful positions found for', gene, sep=' '))
    return(NULL)
  }
  
  clust.tbl$distance = c(0,diff(clust.tbl$pos)) #calculate inter event distance.
  
  #If more than one meaningful positions are found within a 5 aa distance, join them to form a cluster.
  if(nrow(clust.tbl) > 1){
    
    #initialize variables.
    cstart = end = clust.tbl[1,pos]
    n = clust.tbl[1,N]
    cdf = c()
    cluster = 1
    
    #Go through entire table and update variables.
    for(i in 2:nrow(clust.tbl)){
      pos = clust.tbl[i,pos]
      
      d = clust.tbl[i,distance]
      
      if(d < mergeDist){
        end = pos
        n = n + clust.tbl[i,N]
      }else{
        tempdf = data.frame(cluster = paste('cluster', cluster, sep='_'), start = cstart, end = end ,N = n)
        cdf = rbind(cdf, tempdf)
        cstart = end = pos
        n = clust.tbl[i,N]
        cluster = cluster + 1
      }
    }
    cdf = rbind(cdf, data.frame(cluster = paste('cluster', cluster, sep='_'), start = cstart, end = end ,N = n))
  } else {
    cdf = data.frame(cluster = 'cluster_1', start = clust.tbl$pos, end = clust.tbl$pos ,N = clust.tbl$N)
  }
  
  #merge adjacent variants to clusters.
  for(i in 1:nrow(cdf)){
    tempcdf = cdf[i,]
    nonclust.tbl$startDist = nonclust.tbl$pos - tempcdf$start
    nonclust.tbl$endDist = nonclust.tbl$pos - tempcdf$end
    
    merge.adj.to.start = nonclust.tbl[startDist >= -5 & startDist <= 0]
    if(nrow(merge.adj.to.start) > 0){
      tempcdf$start = merge.adj.to.start[which(merge.adj.to.start$startDist == min(merge.adj.to.start$startDist)),pos]
      tempcdf$N = tempcdf$N + sum(merge.adj.to.start$N)
    }
    
    merge.adj.to.end = nonclust.tbl[endDist <= 5 & endDist >= 0]
    if(nrow(merge.adj.to.end) > 0){
      tempcdf$end = merge.adj.to.end[which(merge.adj.to.end$endDist == max(merge.adj.to.end$endDist)),pos]
      tempcdf$N = tempcdf$N + sum(merge.adj.to.end$N)
    }
    cdf[i,] = tempcdf
  }
  cdf$Hugo_Symbol = Hugo_Symbol
  
  #Calcluate cluster score.
  
  total.muts = nrow(prot.dat) #total variants for this gene.
  clusterScores = c()
  
  for(i in 1:nrow(cdf)){
    temp.prot.dat = prot.dat[pos >= as.numeric(cdf$start[i]) & pos <= as.numeric(cdf$end[i])]
    temp.prot.dat.summary = temp.prot.dat[,.N, pos]
    temp.prot.dat.summary[,fraction:= N/total.muts]
    
    peak = temp.prot.dat.summary[N == max(N), pos]
    
    posVector = as.numeric(temp.prot.dat.summary[,pos])
    fractionMutVector = unlist(lapply(posVector, FUN = function(x) temp.prot.dat.summary[pos == x, fraction]))
    distanceVector = suppressWarnings(abs(posVector - peak))
    
    clusterScores = c(clusterScores,  sum( fractionMutVector / (sqrt(2)^ distanceVector)))
    
  }
  cdf$clusterScore = clusterScores
  
  gene.clust.res = data.frame(Hugo_Symbol = Hugo_Symbol, clusters = nrow(cdf), muts_in_clusters = sum(cdf$N), clusterScores = sum(cdf$clusterScore), protLen = protLen)
  return(gene.clust.res)
}

