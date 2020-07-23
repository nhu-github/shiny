# 6 colnames,ORDER_ID,GENE,GENOMIC , DNA_CHANGE, VAR_TYPE ,REPORT_OR_VUS
# NN is the cutoff of mut frequency
# mergenamesbin,merge the genenames based on distance on chromosome
# splitnames, 4 genenames is in a line


chromosomal_changes <- function(var,NN=2,mergenamesbin=15,splitnames=4,
                                expand=NULL,pdfnames="chromosomal_changes",
                                pdfoutput= T,tableoutput=F,
                                pdfwidth=7.5,pdfheight=10,
                                genenamessize=0.6,freqsize=0.8
                                
                                ){
  
  
  mergegenenames<- function(gene1,BIN){
    count <- 1
    for(i in 1:nrow(gene1)){
      if (count ==i){
        #print(count)
        siteN <- gene1$site[count] + BIN
        NN1 <- which(!gene1$site > siteN)
        gene1$GENE[NN1] <- paste0(gene1$GENE[NN1],collapse = ",")
        gene1$FRE[NN1] <- max(gene1$FRE[NN1])
        gene1$ST[NN1] <- (gene1$ST[min(NN1)] + gene1$ST[max(NN1)])/2
        gene1$ED[NN1] <- (gene1$ED[min(NN1)] + gene1$ED[max(NN1)])/2
        gene1$site[NN1] <- max(gene1$site)+1000
        count <- max(NN1)+1
      }
    }
    ge <- unique(gene1[,c(1,2,4,5,6,7)])
    return(ge)
  }
  
  splitgenename <- function(gene1,nn){
    for(i in 1:nrow(gene1)){
      gs <- unlist(strsplit(gene1$GENE[i],","))
      if(length(gs)> nn){
        NNN <- ceiling(length(gs)/nn)
        gs0 <- paste0(na.omit(gs[1:nn]),collapse = ',')
        gs <- gs[-(1:nn)]
        for(x in 1:(NNN-1)){
          gs1 <- paste0(na.omit(gs[1:nn]),collapse = ',')
          gs0 <-  paste0(gs0,"\n",gs1)
          gs <- gs[-(1:nn)]
        }
        gene1$GENE[i]<-gs0
      }
    }
    return(gene1)
  }
  
  # plot
  draw_plot <- function(gene,gene1,gene2,expand=expand,pdfnames=pdfnames,
                        pdfwidth=pdfwidth,pdfheight=pdfheight,
                        genenamessize=genenamessize,freqsize=freqsize,pdfoutput=pdfoutput
  ){
    
    if(!is.null(expand)){
      gene$FRE <-  gene$FRE*expand
      gene1$FRE <- gene1$FRE*expand
      gene2$FRE <- gene2$FRE*expand
    }
    if(pdfoutput==T){
      pdf(paste0(pdfnames,".pdf"), width = pdfwidth, height =pdfheight)
    }
    plot(c(1:100),xlim=c(0,100),ylim=c(0,1100),type="n",xlab = "",ylab = "",axes = F)
    for(i in 1:22){
      if(i==1){yst=0}else{yst=cyto[i-1,]$LOC}
      yed=cyto[i,]$LOC
      cen=cyto[i,]$CEN_L
      if(cyto[i,]$WZ=="z"){xst=50-1.5;xed=50;text(50-5,(yst+yed)/2,labels=i,cex=1);}else{xst=50;xed=50+1.5;text(50+5,(yst+yed)/2,labels=i,cex=1);}
      rect(xst,yed,xed,yst,col = "grey")
      rect(xst,cen,xed,cen,border="red")
    }
    
    # v1 
    # for(i in 1:nrow(gene)){
    #   yst=gene[i,]$ST
    #   yed=gene[i,]$ED
    #   if(gene[i,]$WZ == "z"){
    #     xst=50-10-gene[i,]$FRE;xed=50-10;rect(xst,yed,xed,yst,border="blue",lwd=2);
    #   if(gene[i,]$GENE %in% txgn$GENE){
    #     text(30,(yst+yed)/2,labels=gene[i,]$GENE,col="blue",cex=0.5)
    #     }
    #   }else{
    #     xst=50+10;xed=50+10+gene[i,]$FRE;rect(xst,yed,xed,yst,border="red",lwd=2);
    #     if(gene[i,]$GENE %in% txgn$GENE){
    #       text(70,(yst+yed)/2,labels=gene[i,]$GENE,col="red",cex=0.5);
    #       }
    #   }
    # 
    # }
    
    # z is left,w is right
    for(i in 1:nrow(gene)){
      yst=gene[i,]$ST
      yed=gene[i,]$ED
      if(gene[i,]$WZ == "z"){
        xst=50-10-gene[i,]$FRE;xed=50-10;rect(xst,yed,xed,yst,border="blue",lwd=2);
      }else{
        xst=50+10;xed=50+10+gene[i,]$FRE;rect(xst,yed,xed,yst,border="red",lwd=2);
      }
    }
    
    for(i in 1:nrow(gene1)){
      yst=gene1[i,]$ST
      yed=gene1[i,]$ED
      Len1 <- length(unlist(strsplit(gene1$GENE[i],',')))
      Len2 <- round(gene1$FRE[i])+1
      text(40 -(Len2+Len1*2.8),(yst+yed)/2,labels=gene1[i,]$GENE,col="blue",cex=genenamessize)
    }
    
    for(i in 1:nrow(gene2)){
      yst=gene2[i,]$ST
      yed=gene2[i,]$ED
      Len1 <- length(unlist(strsplit(gene2$GENE[i],',')))
      Len2 <- round(gene2$FRE[i])+1
      text(60 +(Len2+Len1*2.8),(yst+yed)/2,labels=gene2[i,]$GENE,col="red",cex=genenamessize)
    }
    
    # for(i in 1:nrow(gene2)){
    #   yst=gene2[i,]$ST
    #   yed=gene2[i,]$ED
    #   Len <- length(unlist(strsplit(gene2$GENE[i],',')))
    #   if(Len <3){
    #     text(70,(yst+yed)/2,labels=gene2[i,]$GENE,col="red",cex=0.5);
    #   }else if(Len >=3 & Len <8){
    #     text(75,(yst+yed)/2,labels=gene2[i,]$GENE,col="red",cex=0.5)
    #   }else {
    #     text(85,(yst+yed)/2,labels=gene2[i,]$GENE,col="red",cex=0.5)
    #   }
    # }
    
    ##add baseline of cnv
    rect(50-10,0,50-10,1000,border="blue",lwd=2)
    rect(50+10,0,50+10,1000,border="red",lwd=2)
    ##add fre axis
    axismax<- ceiling(max(gene$FRE))
    for(i in pretty(0:axismax,2)){
      rect(50-10-i,1020,50-10-i,1030,border = "black")
      rect(50+10+i,1020,50+10+i,1030,border = "black")
      if(!is.null(expand)){
        text(50-10-i,1070,labels=paste(i/expand,"%",""),cex=freqsize,srt=90,font=0.5)
        text(50+10+i,1070,labels=paste(i/expand,"%",""),cex=freqsize,srt=90,font=0.5)
      }else{
        # text(50-10-i,1070,labels=paste(i,"%",""),cex=0.8,srt=90,font=0.5)
        text(50-10-i,1070,labels=paste(i,"%",""),cex=freqsize,srt=90,font=0.5)
        text(50+10+i,1070,labels=paste(i,"%",""),cex=freqsize,srt=90,font=0.5)
      }
      
    }
    axismax1 <- max(pretty(0:axismax,2))
    rect(50-10,1020,50-10-axismax1,1020,border = "black")
    rect(50+10,1020,50+10+axismax1,1020,border = "black")
    
    if(pdfoutput==T){
      dev.off()
    }
  }
  
  
  

  hg19 <- read.xlsx("./data/hg19_chr.xlsx")
  cyto <- data.frame("CHRO"=paste0("chr",1:22),
                     "LOC"=floor(hg19$cumlen/max(hg19$cumlen[1:23])*1000)[2:23],
                     "CEN_L"=round((hg19$cumlen+hg19$centromerEnd)/max(hg19$cumlen[1:23])*1000)[1:22],
                     "WZ"=rep(c("z","y"),11)
  )
  cs <- max(hg19$cumlen[1:23])
  #gene=read.delim("files/Gene.sumLoc.norm", stringsAsFactors = F)
  #cyto=read.delim("files/chr1_22.sumLoc.norm",stringsAsFactors = F)
  #txgn=read.delim("files/text.genelsts.uniq", stringsAsFactors = F)
  var.dec <- var[grepl('CNV', var$VAR_TYPE),c("ORDER_ID","GENE","GENOMIC", "DNA_CHANGE", "VAR_TYPE","REPORT_OR_VUS")]
  var.dec <- var.dec[var.dec$REPORT_OR_VUS=="report",]
  var.dec <- tidyr::separate(var.dec, GENOMIC, c('Chr', 'Start', 'End'), sep = '-|:')
  var.dec$Chr <- gsub('chr', '', var.dec$Chr)
  var.dec$Chr <-  as.numeric(var.dec$Chr)
  var.dec$Start <- as.numeric(var.dec$Start)
  var.dec$End <- as.numeric(var.dec$End)
  var.dec$type[grepl("AMP|Amp|amp|扩增",var.dec$DNA_CHANGE)] <- "AMP"
  var.dec$type[grepl("DEL|Del|del|缺失",var.dec$DNA_CHANGE)] <- "DEL"
  
  #data.table::setnames(var.dec,"基因","GENE")
  var.dec$WZ <- ifelse(var.dec$Chr %%2==1,"z","y")
  var.dec <- na.omit(var.dec)
  
  hg <- hg19[,c("chrom","cumlen")]
  var.dec1 <- merge(var.dec,hg,by.x="Chr",by.y="chrom")
  var.dec1$ST <- (var.dec1$Start +var.dec1$cumlen)/cs *1000
  var.dec1$EN <- (var.dec1$End +var.dec1$cumlen)/cs *1000
  tb <- unique(var.dec1[,c("ORDER_ID","GENE","type")])
  df <- data.frame(table(tb[,2:3]))
  var.dec1 <- merge(var.dec1,df,by=c("GENE","type"))
  var.dec1 <- unique(var.dec1[,c("GENE","type","Freq","ST","EN","WZ")])
  var.dec1$FRE <- var.dec1$Freq/length(unique(var$ORDER_ID))*100
  var.dec1 <- var.dec1[order(-var.dec1$FRE),]
  data.table::setnames(var.dec1,c("type","Freq","EN"),c("TYPE","CNT","ED"))
  
  CN_information <- var.dec1[,c(1,2,3,7)]
  if(tableoutput==T){
    return(CN_information)
  }
  #write.xlsx(CN_information,"CN_information.xlsx")
  gs <- as.data.frame(var.dec1$GENE[var.dec1$FRE >NN])
  colnames(gs) <- "GENE"
  txgn <- gs
  gene <- var.dec1[var.dec1$GENE %in% txgn$GENE,]
  gene$site <- (gene$ST+gene$ED)/2
  gene1 <- gene[gene$WZ=="z",]
  gene1 <- gene1[order(gene1$site),]
  gene2 <- gene[!gene$WZ=="z",]
  gene2 <- gene2[order(gene2$site),]
  
  gene1 <- mergegenenames(gene1,mergenamesbin)
  gene2 <- mergegenenames(gene2,mergenamesbin)
  
  gene1 <- splitgenename(gene1,splitnames)
  gene2 <- splitgenename(gene2,splitnames)
  
  draw_plot(gene,gene1,gene2,expand=expand,pdfnames=pdfnames,
            pdfwidth=pdfwidth,pdfheight=pdfheight,
            genenamessize=genenamessize,freqsize=freqsize,pdfoutput=pdfoutput
  )
}

# library(openxlsx)
# var <- read.xlsx("E:/song/project/Himalaya/example/example_mutation.xlsx")
# chromosomal_changes(var,pdfnames="chromosomal_changes_example1",pdfoutput = T,tableoutput = T)
# head(var)
# 
# var <- read.xlsx("Variants.xlsx")
# chromosomal_changes(var,pdfnames="chromosomal_origin")
# tt <- chromosomal_changes(var,expand=2,pdfnames="chromosomal_expand",pdfoutput = F,tableoutput = T)
# write.xlsx(tt,"CN_information.xlsx")







