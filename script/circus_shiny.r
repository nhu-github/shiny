library(circlize)
#BiocManager::install("chimeraviz")
library(RCircos)
library(RColorBrewer)


plot_fusion <- function(fus_df,outprefix){
  
  fusion_data <- mydata[mydata$VAR_TYPE %in% c("FUS","FUS2","FUS3"),]
  
  gene_pair <- fus_df[,"GENE_PAIR",drop=F]
  
  gene_table <- separate(gene_pair[,"GENE_PAIR",drop=F],c("GENE_PAIR"),c("gene1","gene2"),sep="-")
  #### remove one gene and Intergenic
  gene_table <- gene_table[!is.na(gene_table$gene2),]
  gene_table <- gene_table[!gene_table$gene1 %in% c("Intergenic","基因间区"),]
  gene_table <- gene_table[!gene_table$gene2 %in% c("Intergenic","基因间区"),]
  
  ###### plot gene pair circos
  circos_table <- table(gene_table$gene1,gene_table$gene2)
  circos_table <- as.data.frame.matrix(circos_table)
  circos_matrix <- as.matrix(circos_table)
  circos_matrix[1,] <- 2
  circos_matrix[2,] <- 3
  circos_matrix[3,] <- 4
  circos_matrix[4,] <- 5
  
  col_set <- RColorBrewer::brewer.pal(n = 12, name = "Set3")
  grid.col[colnames(circos_matrix)] <- col_set[1:dim(circos_matrix)[1]]
  grid.col[rownames(circos_matrix)] <- "grey"
  circos.par(gap.degree = c(rep(2, ncol(circos_matrix)-1), 10, rep(2, nrow(circos_matrix)-1), 10),
             start.degree = 90)
  
  chordDiagram(circos_matrix, directional = TRUE, diffHeight = 0.06, grid.col = grid.col, transparency = 0.5)
  circos.clear()
 
}



get_track_position <- function(track.num){
#### get track position
inside.pos=NULL
outside.pos=NULL 
side = "in" 
boundary <- RCircos.Get.Plot.Boundary(track.num=track.num, side, inside.pos,
                                      outside.pos, FALSE);
outerPos <- boundary[1]-0.03;
innerPos  <- boundary[2]-0.03;
result <- list(outerPos=outerPos,innerPos=innerPos)
return(result)
}

  
plot_circos_gene <- function(mut_df,genomic_col,mut_type_col,gene_col,outprefix,plot_mut_type=c("snv","cnv","fusion"),width=8,height=8,
                             gene_pair_col=NULL,amp_del_col=NULL,chr.exclude=NULL,single=NULL,log=TRUE,lable_size=0.3){
  
  #plot_types <- strsplit(plot_mut_type,split=",")[[1]]
  
  plot_types <- plot_mut_type
  num_plot <- length(plot_types)
  if(length(setdiff(plot_types,c("snv","cnv","fusion")))>0){
    stop("plot_mut_type must be in snv,cnv,fusion")
  }
  
  ## set total circos tracks
  if(num_plot==3){
    tracks.inside=7
    tracks.inside_list = c(1,2,3,4,6,7)
    names(tracks.inside_list) <- c("cnv_g","cnv","snv_g","snv","fusion_g","fusion")
  }else if(num_plot==2){
    tracks.inside=5
    tracks.inside_list = c(2,3,6,7)
    if(length(setdiff(plot_types,c("cnv","snv")))==0){
      names(tracks.inside_list) <- c("cnv_g","cnv","snv_g","snv")
    }else if(length(setdiff(plot_types,c("cnv","fusion")))==0){
      names(tracks.inside_list) <- c("cnv_g","cnv","fusion_g","fusion")
    }else if(length(setdiff(plot_types,c("snv","fusion")))==0){
      names(tracks.inside_list) <- c("snv_g","snv","fusion_g","fusion")
    }
  }else if(num_plot==1){
    tracks.inside=3
    tracks.inside_list = c(2,3)
    g_name <- paste(plot_types,"g",sep="_")
    name <- plot_types
    names(tracks.inside_list) <- c(g_name,name)
  }
  tracks.outside=0
  
  data(UCSC.HG19.Human.CytoBandIdeogram)
  
  gene_location <- read.table("./data/gene_location.txt",sep="\t",header=FALSE,stringsAsFactors=FALSE)
  colnames(gene_location) <- c("chr","start","end","strain","name")
  

  # chr.exclude<-NULL #设置不显示的染色体，如 c(1,3)
  # single <- c("chr1") # 只显示一条染色体
  all_chr<-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)
  all_chr <- paste("chr",all_chr,sep="")
  all_chr <- c(all_chr,"chrX","chrY")
  if(!is.null(single)){
    chr.exclude <- all_chr[!all_chr %in% single]
  }
  
  gene_location <- gene_location[gene_location$chr %in% all_chr,]  
  
  ######## SNV ########
  if("snv" %in% plot_types){
    #print("SNV>>>>>>>>>>>")
    SNv_df <- mut_df[mut_df[,mut_type_col] %in% c("SNV","LONG","BLD"),]
    split1 <- tidyr::separate(SNv_df[,genomic_col,drop=F],genomic_col,c("chr","position"),sep=":")
    split2 <- tidyr::separate(split1,c("position"),c("Start_Position","End_Position"),sep="-")
    SNv_df[,c("chr","Start_Position","End_Position")]<- split2 
    
    gene_count <- table(SNv_df[,gene_col])
    gene_count <- as.data.frame(gene_count)
    
    SNV_merge<- merge(gene_count,gene_location,by.x="Var1",by.y="name",sort = FALSE)
    SNV_merge<- SNV_merge[SNV_merge$chr %in% all_chr,]
    
    SNV_hist <- SNV_merge[,c("chr","start","end","Freq","Var1")]
    SNV_hist$start <- as.integer(SNV_hist$start)
    SNV_hist$end <- as.integer(SNV_hist$end)
    SNV_hist$Freq <- as.numeric(SNV_hist$Freq)
    
    if(!is.null(single)){
      SNV_hist <- SNV_hist[SNV_hist$chr == single,]
    }
    if(!is.null(chr.exclude)){
      SNV_hist <- SNV_hist[!SNV_hist$chr  %in% chr.exclude,]
    }
    
    ### add high snv gene
    mean_snv <- mean(SNV_hist$Freq)
    snv_plot_gene <- SNV_hist[SNV_hist$Freq >mean_snv,]
  }
  
  #################### CNV
  if("cnv" %in% plot_types){
    if(is.null(amp_del_col)){
      stop("when plot cnv,amp_del_col name with content Amp and Del must be given!")
    }
    
    CNV_data <- mut_df[mut_df[,mut_type_col] == "CNV",c(gene_col,amp_del_col)]
    #head(CNV_data)
    colnames(CNV_data) <- c("Hugo_Symbol","DNA_CHANGE")
    CNV_num <- table(CNV_data$Hugo_Symbol,CNV_data$DNA_CHANGE)
    CNV_num <- as.data.frame.matrix(CNV_num)
    CNV_value <- as.data.frame(matrix(numeric(0),ncol=1))
    colnames(CNV_value) <- "value"
  
    amp_col <- colnames(CNV_num)[grep("Amp",colnames(CNV_num))]
    del_col <- colnames(CNV_num)[grep("Del",colnames(CNV_num))]
    
    for(eachamp in amp_col){
      CNV_tmp <- CNV_num[CNV_num[eachamp]!=0,eachamp,drop=FALSE]
      colnames(CNV_tmp) <- c("value")
      CNV_value <- rbind(CNV_value,CNV_tmp)
    }
    for(eachdel in del_col){
      CNV_tmp <- CNV_num[CNV_num[eachdel]!=0,eachdel,drop=FALSE]
      colnames(CNV_tmp) <- c("value")
      CNV_tmp$value <- -CNV_tmp$value
      CNV_value <- rbind(CNV_value,CNV_tmp)
    }
    
    CNV_value$gene <- rownames(CNV_value)
    
    CNV_hist <- merge(CNV_value,gene_location,by.x="gene",by.y="name",sort=FALSE)
    CNV_hist <- CNV_hist[,c("chr","start","end","value","gene")]
    CNV_hist$start <- as.integer(CNV_hist$start)
    CNV_hist$end <- as.integer(CNV_hist$end)
    CNV_hist$Freq <- as.numeric(CNV_hist$value)
    CNV_hist<- CNV_hist[CNV_hist$chr %in% all_chr,]
    
    if(!is.null(single)){
      CNV_hist <- CNV_hist[CNV_hist$chr == single,]
    }
    if(!is.null(chr.exclude)){
      CNV_hist <- CNV_hist[!CNV_hist$chr  %in% chr.exclude,]
    }
    
    mean_cnv <- mean(abs(CNV_hist$Freq))
    cnv_plot_gene <- CNV_hist[abs(CNV_hist$Freq) >mean_cnv,]
  }
  
  
  ##### fusion
  if("fusion" %in% plot_types){
    
    if(is.null(gene_pair_col)){
      stop("when plot fusion,gene_pair_col name with content like ALK-EML4 must be given!")
    }
    
    fusion_data <- mut_df[mut_df[,mut_type_col] %in% c("FUS","FUS2","FUS3","Fusion","FUSION","rearrangement","fusion"),]
    
    gene_pair <- fusion_data[,gene_pair_col,drop=F]
    
    gene_table <- tidyr::separate(gene_pair[,gene_pair_col,drop=F],gene_pair_col,c("gene1","gene2"),sep="-")
    #### remove one gene and Intergenic
    gene_table <- gene_table[!is.na(gene_table$gene1),]
    gene_table <- gene_table[!is.na(gene_table$gene2),]
    gene_table <- gene_table[!gene_table$gene1 %in% c("Intergenic","基因间区","intergenic"),]
    gene_table <- gene_table[!gene_table$gene2 %in% c("Intergenic","基因间区","intergenic"),]
    
  
    
    # gene_table <- tidyr::separate(gene_table,"gene1",c("gene1","tmp1"),sep="@",remove=FALSE)
    # gene_table <- gene_table[is.na(gene_table$tmp1),]
    # gene_table <- tidyr::separate(gene_table,"gene2",c("gene2","tmp2"),sep="@",remove=FALSE)
    # gene_table <- gene_table[is.na(gene_table$tmp2),]
    
    gene_table$gene1 <-  gsub(pattern = '@', replacement = '', x = gene_table$gene1, fixed = TRUE)
    gene_table$gene2 <-  gsub(pattern = '@', replacement = '', x = gene_table$gene2, fixed = TRUE)
    
    dim(gene_table)
  
    
    c_merge <- cbind(gene_table["gene1"],gene_table["gene2"])
    
    fusion_table <- merge(c_merge,gene_location,by.x="gene1",by.y="name",sort = FALSE)
    fusion_table <- merge(fusion_table,gene_location,by.x="gene2",by.y="name",sort = FALSE)
    
    
    colnames(fusion_table) <- c("gene2","gene1","chr1","start1","end1","strain1","chr2","start2","end2","strain2")
    fusion_link <- fusion_table[,c("chr1","start1","end1","chr2","start2","end2")]
    fusion_link$start1 <- as.integer(fusion_link$start1)
    fusion_link$end1 <- as.integer(fusion_link$end1)
    fusion_link$start2 <- as.integer(fusion_link$start2)
    fusion_link$end2 <- as.integer(fusion_link$end2)
    
    #RCircos.Link.Plot(fusion_link,track.num=6,by.chromosome=TRUE,is.sorted=FALSE) 
    
    
    f1_df <- fusion_table[,c("chr1","start1","end1","gene1")]
    colnames(f1_df) <- c("chr","start","end","gene")
    f2_df <- fusion_table[,c("chr2","start2","end2","gene2")]
    colnames(f2_df) <- c("chr","start","end","gene")
    fusion_plot_gene <- rbind(f1_df,f2_df)
    fusion_plot_gene$start <- as.integer(fusion_plot_gene$start )
    fusion_plot_gene$end <- as.integer(fusion_plot_gene$end )
    
    # fusion_plot_gene <- tidyr::separate(fusion_plot_gene,"chr",c("chr_str","chr_num"),sep="r",remove=FALSE)
    # fusion_plot_gene$chr_num <- 
    # plyr::arrange(fusion_plot_gene,chr_num)
  }
  
  # if(file_test("-d",outprefix)){
  #   if(!file.exists(file.path(outprefix))){dir.create(outprefix,recursive = TRUE)}
  #   pdfname <- paste(outprefix,"/circos.pdf",sep="")
  # }else{
  #   path_split <- strsplit(outprefix,split="/")[[1]]
  #   file_dir <- paste(path_split[1:length(path_split)-1],collapse="/")
  #   if(!file.exists(file.path(file_dir))){dir.create(file_dir,recursive = TRUE)}
  #   pdfname <- paste(outprefix,".circos.pdf",sep="")
  # }
  # 
  
  #pdf(file=pdfname, height=height, width=width, compress=TRUE);
  ## initialized RCircos
  print("initialized>>>>>>>>>>>")
  RCircos.Set.Core.Components(cyto.info=UCSC.HG19.Human.CytoBandIdeogram,chr.exclude=chr.exclude,tracks.inside=tracks.inside,tracks.outside=tracks.outside)
  #CairoPDF(outF,width = 8,height = 8)
  RCircos.Set.Plot.Area()
  RCircos.Chromosome.Ideogram.Plot()
  
  mycolor <- rgb(251,255,185,alpha=90,max=255)
  ###plot cnv
  if("cnv" %in% plot_types){
    print("plot cnv")
    #mycolor <- rgb(251,255,185,alpha=90,max=255)
      
    genomic.columns <- 3
    data.col <- 4
    track.num <- tracks.inside_list[names(tracks.inside_list)=="cnv"]
    add_histogram(hist.data=CNV_hist,track.num=track.num,track.colors=mycolor,genomic.columns=genomic.columns,data.col=data.col,plot_type = "CNV")  
    
    ### add_cnv gene
    name.col <- 5 # 
    side <- "in"
    track.num <- tracks.inside_list[names(tracks.inside_list)=="cnv_g"]
    #RCircos.Gene.Name.Plot(cnv_plot_gene, name.col,track.num,side,is.sorted=FALSE) ## 基因
    my_gene_name(cnv_plot_gene, name.col,track.num,side,is.sorted=FALSE,lable_size=lable_size)
  }
  
  ###plot snv
  if("snv" %in% plot_types){
    print("plot SNV")
    genomic.columns <- 3
    data.col <- 4
    track.num <- tracks.inside_list[names(tracks.inside_list)=="snv"]
    
    add_histogram(hist.data=SNV_hist,track.num=track.num,genomic.columns=genomic.columns,data.col=data.col) 
    
    name.col <- 5 # 
    side <- "in"
    track.num <- tracks.inside_list[names(tracks.inside_list)=="snv_g"]
    #RCircos.Gene.Name.Plot(snv_plot_gene, name.col,track.num,side,is.sorted=FALSE) ## 基因
    my_gene_name(snv_plot_gene, name.col,track.num,side,is.sorted=FALSE,lable_size=lable_size)
  }
  
  ####plot fusion
  if("fusion" %in% plot_types){
    print("plot fusion")
    track.num= tracks.inside_list[names(tracks.inside_list)=="fusion"]
    RCircos.Link.Plot(fusion_link,track.num=track.num,by.chromosome=TRUE,is.sorted=FALSE) 
    
    name.col <- 4 # 
    side <- "in"
    track.num <- tracks.inside_list[names(tracks.inside_list)=="fusion_g"]
    #RCircos.Gene.Name.Plot(fusion_plot_gene, name.col,track.num,side,is.sorted=FALSE) ## 基因
    my_gene_name(unique(fusion_plot_gene), name.col,track.num,side,is.sorted=FALSE,lable_size=lable_size-0.05)
    
  }
  
  legend("topright",bty = "n",legend=c("CNV","SNV","fusion-intrachromosomal","fusion-interchromosomal"), col=c(mycolor,"#CCFFFF","red","blue"), 
         pch=c(15,15,NA,NA),lty=c(0,0,1,1),pt.bg = "black")  
  
  #dev.off()
  
}


##RCircos.Gene.Name.Plot()
my_gene_name <-function (gene.data = NULL, name.col = NULL, track.num = NULL, 
          side = "in", inside.pos = NULL, outside.pos = NULL, 
          genomic.columns = 3, is.sorted = FALSE,lable_size=0.3) 
{
  if (is.null(gene.data)) 
    stop("Genomic data missing in RCircos.Gene.Name.Plot().\n")
  if (is.null(genomic.columns)) 
    stop("Missing number of columns for genomic position.\n")
  if (is.null(name.col) || name.col <= genomic.columns) 
    stop("Data column must be ", genomic.columns + 
           1, " or bigger.\n")
  RCircos.Pos <- RCircos.Get.Plot.Positions()
  RCircos.Par <- RCircos.Get.Plot.Parameters()
  textColors <- RCircos.Get.Plot.Colors(gene.data, RCircos.Par$text.color)
  boundary <- RCircos.Get.Plot.Boundary(track.num, side, inside.pos, 
                                        outside.pos, FALSE)
  gene.data <- RCircos.Get.Single.Point.Positions(gene.data, 
                                                  genomic.columns)
  gene.data <- RCircos.Get.Gene.Label.Locations(gene.data, 
                                                genomic.columns, is.sorted)
  rightSide <- nrow(RCircos.Pos)/2
  thePoints <- as.numeric(gene.data[, ncol(gene.data)])
  if (side == "in") {
    labelPos <- boundary[1]
    textSide <- rep(4, nrow(gene.data))
    textSide[thePoints <= rightSide] <- 2
  }
  else {
    labelPos <- boundary[2]
    textSide <- rep(2, nrow(gene.data))
    textSide[thePoints <= rightSide] <- 4
  }
  for (aText in seq_len(nrow(gene.data))) {
    geneName <- as.character(gene.data[aText, name.col])
    rotation <- RCircos.Pos$degree[thePoints[aText]]
    text(RCircos.Pos[thePoints[aText], 1] * labelPos, RCircos.Pos[thePoints[aText], 
                                                                  2] * labelPos, label = geneName, pos = textSide[aText], 
         cex = lable_size, srt = rotation, offset = 0, 
         col = textColors[aText])
  }
}



add_histogram <- function(hist.data,track.num,track.colors="#CCFFFF",genomic.columns=3,data.col=4,
                          max.value=NULL,min.value=NULL,log=TRUE,plot_type="n_CNV"){

  #
  
  track_positon <- get_track_position(track.num=track.num)
  outerPos <- track_positon$outerPos
  innerPos <- track_positon$innerPos
  if(plot_type == "CNV"){
  outerPos <- track_positon$outerPos + 0.03
  innerPos  <- track_positon$innerPos + 0.03
  }
  
  RCircos.Track.Outline(outerPos, innerPos, 1,track.colors=track.colors);  
  
  hist.data <- RCircos.Get.Single.Point.Positions(hist.data,genomic.columns);
  
  hist.data <- hist.data[hist.data$Freq !=0,]
  
  #locations <- RCircos.Get.Start.End.Locations(hist.data,RCircos.Par$hist.width)
  RCircos.Par <- RCircos.Get.Plot.Parameters()
  locations <- get_locations(hist.data,RCircos.Par$hist.width)
  
  #    histgram colors and height
  #    =========================================================
  #
  #histColors <- RCircos.Get.Plot.Colors(hist.data, RCircos.Par$hist.color);
  hist.data$color <- "red"
  histColors <- ifelse(hist.data$Freq < 0, "green", "red")
  
  
  histValues <- as.numeric(hist.data[, data.col]);
  if(log){
    histValues <- log2(abs(histValues)) + 0.8
  }
  if(is.null(max.value) || is.null(min.value)) {
    max.value <- max(abs(histValues));
    min.value <- min(abs(histValues));
  } else {
    if(min.value > max.value) stop("min.value > max.value.")
  }
  plot_position <- outerPos-innerPos
  # if(plot_type=="CNV"){
  #   plot_position <- (outerPos-innerPos)/2
  # }
  histHeight <- RCircos.Get.Data.Point.Height(histValues, 0, max.value, plot.type="points", plot_position);
  RCircos.Pos <- RCircos.Get.Plot.Positions()
  
  for(aPoint in seq_len(nrow(hist.data)))
  {
    #print(aPoint)
    #if(plot_type=="CNV"){
    #height <- innerPos+(outerPos-innerPos)/2 + histHeight[aPoint];

    height <- innerPos + histHeight[aPoint];
    
    theStart <- locations[aPoint, 1];
    theEnd <- locations[aPoint, 2];
    #print(height)
    
    #    Plot rectangle with specific height for each data point
    #    =========================================================
    #
    polygonX <- c(RCircos.Pos[theStart:theEnd,1]*height, RCircos.Pos[theEnd:theStart,1]*innerPos);
    polygonY <- c(RCircos.Pos[theStart:theEnd,2]*height, RCircos.Pos[theEnd:theStart,2]*innerPos);
    polygon(polygonX, polygonY, col=histColors[aPoint], border=NA);
  } 
}


################## public function error, modify
get_locations <- function (plot.data, plot.width) 
{
  RCircos.Cyto <- RCircos.Get.Plot.Ideogram()
  dataChroms <- as.character(plot.data[, 1])
  chromosomes <- unique(dataChroms)
  cyto.chroms <- as.character(RCircos.Cyto$Chromosome)
  point.loc <- as.numeric(plot.data$Location)
  locations <- cbind(point.loc - plot.width, point.loc + plot.width)
  for (aChr in seq_len(length(chromosomes))) {
    cyto.rows <- which(cyto.chroms == chromosomes[aChr])
    chr.start <- min(RCircos.Cyto$StartPoint[cyto.rows])
    chr.end <- max(RCircos.Cyto$EndPoint[cyto.rows])
    data.rows <- which(dataChroms == chromosomes[aChr])
    start.outliers <- which(locations[data.rows, 1] < chr.start)
    if (length(start.outliers) > 0) 
      locations[data.rows[start.outliers], 1] <- chr.start
    end.outliers <- which(locations[data.rows, 2] > chr.end)
    if (length(end.outliers) > 0) 
      #locations[data.rows[start.outliers], 2] <- chr.end
      locations[data.rows[end.outliers], 2] <- chr.end
  }
  return(locations)
}



