


plot_CNV_seqment <- function(cnv_df,cli_df=NULL,sample_col,gene_col,genomic_col,cnv_type_col,outprefix,sample_num=NULL,ref.build='hg19',
                             cytobandTxtSize=0.6,genesize= 0.6,lablecut=0,top=10,y_lim=NULL,def_ystep=NULL,all_gene=FALSE,
                             width=8,height=6,Writetable=F,Plotpdf=T){

  cnv_data <- cnv_df[,c(sample_col,gene_col,genomic_col,cnv_type_col)]  
  colnames(cnv_data) <- c("SMP_ID","gene","GENOMIC","CNV_type")
  cnv_data <- cnv_data[!duplicated(cnv_data[,c("SMP_ID","gene","CNV_type")]),]
  
  ref_hg <- ref.build
  if(is.null(sample_num)){
    sample_num <- length(unique(cnv_data$SMP_ID))
  }
  color = c("royalblue", "maroon")
  
  if(!is.null(cli_df)){
    samples <- unique(as.character(cli_df[,c(sample_col)]))
    sample_num <- length(samples)
    cnv_data <- cnv_data[cnv_data[,c("SMP_ID")] %in% samples,]
  }

  ###### plot segment CNV

  new_cnv_tmp <- tidyr::separate(cnv_data,"GENOMIC",sep=":",c("chr","position"))
  new_cnv <- tidyr::separate(new_cnv_tmp,"position",sep="-",c("chr_start","chr_end"))
  new_cnv$type <- "others"
  has_amp <- grep("Ampl|ampl",new_cnv$CNV_type)
  has_del <- grep("Dele|dele",new_cnv$CNV_type)
  if(length(has_amp)>0){
    new_cnv[grep("Ampl|ampl",new_cnv$CNV_type),]$type <- "Amp"
  }
  if(length(has_del)>0){
    new_cnv[grep("Dele|dele",new_cnv$CNV_type),]$type <- "Del"
  }
  #g[,End_Position := sapply(strsplit(x = g$loc, split = '-'), '[', 2)]
  
  
  ########### use_patient_num
  new_cnv1 <- unique(new_cnv[,c("SMP_ID","gene","type")])
  amp_del <- as.data.frame.matrix(table(new_cnv$gene,new_cnv$type))
  
  
  if(length(has_amp)>0){
    amp_del$Amp <- amp_del$Amp/sample_num
  }
  if(length(has_del)>0){
    amp_del$Del <- amp_del$Del/sample_num
  }
  
  ylim_up <- 1
  ylim_down <- -1
  if("Amp" %in% colnames(amp_del)){
    ylim_up <- max(amp_del[,"Amp"])
  }else{
    ylim_up <- 0.1
  }
  if("Del" %in% colnames(amp_del)){
    ylim_down <- -max(amp_del[,"Del"])
  }else{
    ylim_down <- -0.1
  }
  
  
  y_value <- min(c(ylim_up,abs(ylim_down)))
  if(y_value<1){
    ytmp1 <- round(ylim_down,1)- 0.1
    ytmp2 <- round(ylim_up,1)+ 0.1
    ystep <- round((ytmp2-ytmp1)/8,2)
    y_at <- seq(ytmp1,ytmp2,by=ystep)
    text_add <- ylim_up/5
  }
  amp_del$gene <- rownames(amp_del)
  final_table <- merge(new_cnv,amp_del,by="gene",sort=FALSE)
  final_table$chr_start <- as.integer(final_table$chr_start)
  final_table$chr_end <- as.integer(final_table$chr_end)
  
  ###### plot 
  
  
  if(ref.build == 'hg19'){
    chr.lens = c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663,
                 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540,
                 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566,
                 155270560, 59373566)
  } else if(ref.build == 'hg18'){
    chr.lens = c(247249719, 242951149, 199501827, 191273063, 180857866, 170899992,
                 158821424, 146274826, 140273252, 135374737, 134452384, 132349534,
                 114142980, 106368585, 100338915, 88827254, 78774742, 76117153,
                 63811651, 62435964, 46944323, 49691432, 154913754, 57772954)
  } else if(ref.build == 'hg38'){
    chr.lens = c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
                 159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
                 114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
                 58617616, 64444167, 46709983, 50818468, 156040895, 57227415)
  }else{
    stop("ref.build can only be hg18, hg19 or hg38")
  }
  
  
  chr.lens.cumsum = cumsum(chr.lens)
  #nchrs = length(unique(gis.scores$Chromosome))
  chr.labels= c(1:22, 'X', 'Y')
  chr.tbl = data.table::data.table(chr = chr.labels, start = c(1, chr.lens.cumsum[1:length(chr.lens.cumsum)-1]), end = chr.lens.cumsum)
  chr.tbl$color = rep(c('black','white'), length=nrow(chr.tbl))
  
  chr.tbl$chrname <- paste("chr",chr.tbl$chr,sep="")
  plot_table <- merge(final_table,chr.tbl,by.x="chr",by.y="chrname",sort=FALSE)
  plot_table$chr_start <- as.integer(plot_table$chr_start)
  plot_table$chr_end <- as.integer(plot_table$chr_end)
  plot_table$total_start <- plot_table$start + plot_table$chr_start -1
  plot_table$total_end <- plot_table$start + plot_table$chr_end -1
  
  min_y1 <- min(plot_table$Amp[plot_table$Amp >0])
  min_y2 <- min(plot_table$Del[plot_table$Del >0])
  if(min_y1>0 & min_y2>0){
    min_y <- min(c(min_y1,min_y2))
  }else if (min_y1>0){
    min_y <- min_y1
  }else if(min_y2>0){
    min_y <- min_y2
  }
  
  rect_dis1 <- ylim_up/8
  rect_dis2 <- min_y*0.9
  rect_dis <- min(c(rect_dis1,rect_dis2))
  if(rect_dis < ylim_up/30){
    rect_dis <- ylim_up/15
  }
  
  if(file_test("-d",outprefix)){
    if(!file.exists(file.path(outprefix))){dir.create(outprefix,recursive = TRUE)}
    pdfname <- paste(outprefix,"/cnv.seqment.pdf",sep="")
  }else{
    path_split <- strsplit(outprefix,split="/")[[1]]
    file_dir <- paste(path_split[1:length(path_split)-1],collapse="/")
    if(!file.exists(file.path(file_dir))){dir.create(file_dir,recursive = TRUE)}
    pdfname <- paste(outprefix,".cnv.seqment.pdf",sep="")
  }
  
  if(Plotpdf==T){
    #pdf(pdfname,width=width,height=height)
    
    par(mar = c(1, 4, 2, 1))
    if(is.null(y_lim)){
      y_expand1 <- ytmp1*1.5
      y_expand2 <- ytmp2*1.5
      y_lim <- c(y_expand1,y_expand2)
    }else{
      y_lim <- strsplit(y_lim,",")[[1]]
      y_lim <- as.numeric(y_lim)
      
    }
    
    
    if(is.null(def_ystep)){
      def_ystep <- ystep
    }
    plot(NA, NA, xlim = c(0, chr.lens.cumsum[length(chr.lens.cumsum)]), ylim = y_lim,
         axes = FALSE, xlab = NA, ylab = NA)
    #abline(v = chr.tbl$end, h = y_lims, lty = 2, col = grDevices::adjustcolor("gray70", 0.25))
    #axis(side = 2, at =seq(floor(ylim_down),ceiling(ylim_up),by=ystep) , las = 2)
    axis(side = 2, at = round(seq(y_lim[1],y_lim[2],by=def_ystep),3), las = 2)
    axis(side = 2, at =0 , las = 2)
    mtext(text = "mutation percent", side = 2, line = 3, cex = 1.2)
    
    out_table <-  data.frame(chr = character(0), gene = character(0), freq = numeric(0))
    all_gene_table <- data.frame(chr = character(0), gene = character(0), freq = numeric(0))
    
    if("Amp" %in% colnames(final_table)){
      amp_table <- unique(plot_table[plot_table$type=="Amp",c("chr","gene","Amp","total_start","total_end")])
      amp_table <- plyr::arrange(amp_table,Amp,decreasing=TRUE)
      
      if(dim(amp_table)[1]>top){
        segments(x0 = amp_table[1:top,]$total_start, y0 = rect_dis, x1 = amp_table[1:top,]$total_end,
                 y1 = amp_table[1:top,]$Amp + rect_dis, col = color[2])
        
      }else{
        segments(x0 = amp_table$total_start, y0 = rect_dis, x1 = amp_table$total_end,
                 y1 = amp_table$Amp + rect_dis, col = color[2])
      }

      #segments(x0 = amp_table$total_start, y0 = rect_dis, x1 = amp_table$total_end,
      #         y1 = amp_table$Amp + rect_dis, col = color[2])
      tmp_table <- amp_table
      colnames(tmp_table) <- c("chr","gene","freq")
      all_gene_table <- rbind(all_gene_table,tmp_table)
      if(lablecut==0){
        if(dim(amp_table)[1]>top){
          text_plot <- amp_table[1:top,]
        }else{
          text_plot <- amp_table
        }
        wordcloud::textplot(x = text_plot$total_start, y = text_plot$Amp+text_add,
                            words = text_plot$gene, new = FALSE, font = 3, cex = genesize)
        out_df <- text_plot[,c("chr","gene","Amp")]
        colnames(out_df) <- c("chr","gene","freq")
        out_table <- rbind(out_table,out_df)
      }else{
        text_plot <- amp_table[amp_table$Amp >lablecut,]
        wordcloud::textplot(x = text_plot$total_start, y = text_plot$Amp+text_add,
                            words = text_plot$gene, new = FALSE, font = 3, cex = genesize)  
        out_df <- text_plot[,c("chr","gene","Amp")]
        colnames(out_df) <- c("chr","gene","freq")
        out_table <- rbind(out_table,out_df)
      }
    }
    if("Del" %in% colnames(final_table)){
      del_table <- unique(plot_table[plot_table$type=="Del",c("chr","gene","Del","total_start","total_end")])
      del_table <- plyr::arrange(del_table,Del,decreasing=TRUE)
      
      if(dim(del_table)[1]>top){
        segments(x0 =  del_table[1:top,]$total_start, y0 = -rect_dis, x1 =  del_table[1:top,]$total_end,
                 y1 = -del_table[1:top,]$Del - rect_dis, col = color[1])
      }else{
        segments(x0 = del_table$total_start, y0 = -rect_dis, x1 = del_table$total_end,
                 y1 = - del_table$Del-rect_dis, col = color[1])
      }

      #segments(x0 = del_table$total_start, y0 = -rect_dis, x1 = del_table$total_end,
      #        y1 = - del_table$Del-rect_dis, col = color[1])

      tmp_table <- del_table
      colnames(tmp_table) <- c("chr","gene","freq")
      tmp_table$freq <- -tmp_table$freq
      all_gene_table <- rbind(all_gene_table,tmp_table)
      if(lablecut==0){
        if(dim(del_table)[1]>top){
          text_plot <- del_table[1:top,]
        }else{
          text_plot <- del_table
        }
        wordcloud::textplot(x = text_plot$total_start, y = -text_plot$Del-text_add,
                            words = text_plot$gene, new = FALSE, font = 3, cex = genesize)
        
        out_df <- text_plot[,c("chr","gene","Del")]
        colnames(out_df) <- c("chr","gene","freq")
        out_df$freq <- -out_df$freq
        out_table <- rbind(out_table,out_df)
        
      }else{
        text_plot <- del_table[del_table$Del >lablecut,]
        colnames(text_plot) <- c("chr","gene","freq")
        out_table <- rbind(out_table,text_plot)
        wordcloud::textplot(x = text_plot$total_start, y = -text_plot$Del-text_add,
                            words = text_plot$gene, new = FALSE, font = 3, cex = genesize)   
        out_df <- text_plot[,c("chr","gene","Del")]
        colnames(out_df) <- c("chr","gene","freq")
        out_df$freq <- - out_df$freq
        out_table <- rbind(out_table,out_df)
      }
    }
    
    
    rect(xleft = chr.tbl$start, ybottom = -rect_dis, xright = chr.tbl$end,
         ytop = rect_dis, col = chr.tbl$color)
    
    text(y = 0, x = apply(chr.tbl[,2:3], 1, mean), labels = chr.tbl$chr,
         cex = cytobandTxtSize, col = c('white', 'black'))
    #dev.off()
  }
  
  
  if(Writetable==T){
    if(!all_gene){
      filename <- paste(outprefix,".cnv.top.xls",sep="")
      #write.table(out_table,file=filename,sep="\t",row.names = FALSE)
      return(out_table)
    }else{
      filename <- paste(outprefix,".cnv.allgene.xls",sep="")
      #write.table(all_gene_table[,c("chr","gene","freq")],file=filename,sep="\t",row.names = FALSE) 
      return(all_gene_table[,c("chr","gene","freq")])
    }
  }
  message("CNV seqment finish...")
}










