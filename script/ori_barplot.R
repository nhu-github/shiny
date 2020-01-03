#' @include CountVariants.R ori_tools.R
#' @import ggplot2
#' @import ggsci
#' @import tidyr
#' @import ggpubr
NULL
#' Drawing bar plot
#'
#' @inheritParams CountVariants
#' @param varorder c('Fusion/Rearrangement', 'Substitution/Indel', 'Gene Amplification', 'Gene Homozygous Deletion', 'Truncation')
#' @param byorder Orders of column \code{\link{by}}.
#' @param gs The gene sets displayed in the plot. Allowed value are numberic or vector.
#' @param color The color palette to be used for coloring or filling by groups. Allowed values include "grey" for grey color palettes; brewer palettes e.g. "RdBu", "Blues", ...; or custom color palette e.g. c("blue", "red"); and scientific journal palettes from ggsci R package, e.g.: "npg", "aaas", "lancet", "jco", "ucscgb", "uchicago", "simpsons" and "rickandmorty". Can be also a numeric vector of length(groups); in this case a basic color palette is created using the function palette.
#' @param ytype y axix type to be show. Allower value include "counts" and "percentage".
#' @param outpdf filename
#'
#' @return A ggplot2 object or pdf
#' @examples
#'
#'
#'
#' @export

ori_barplot <- function(file,
                    id = 'ORDER_ID',
                    gene = 'GENE',
                    vartype = 'VarType',
                    varorder = c('Fusion/Rearrangement', 'Substitution/Indel', 'Gene Amplification', 'Gene Homozygous Deletion', 'Truncation'),
                    by = NULL,
                    byorder = NULL,
                    gs = 20,
                    color = NULL,
                    ytype = 'percentage',
                    outpdf = NULL
                    ){
  if(is.null(color)){
    color <- c("Fusion/Rearrangement" = "#EEEE00",      ##yellow
             "Substitution/Indel" = "#228B22",      ##green
             "Gene Amplification" = "#EE0000",        ##red
             "Gene Homozygous Deletion" = "#0000EE",  ##blue
             "Truncation" = "#8E388E")                ##purple
             #'Splicing' = '#B2DF8A')
  }
  if(length(setdiff(varorder, names(color)) > 0)){
    stop('Please check variant type or specify color argument(color)')
  }

  oridata <- CountVariants(xlsxfile = file, id = id, gene = gene,
                        vartype = vartype, by = by, outxlsx = 'list')
  if(is.null(byorder)) byorder <- names(oridata)[2:length(oridata)]

  # plot bar
  topdata <- function(cvlist, gkeep = gs){
    #cvlist, statistical table by CountVariants.R
    pfre <- cvlist[[1]]

    if(is.numeric(gs)){
      gkeep <- pfre[[1]][1:gs]
    }else if(is.vector(gs)){
      gkeep <- gs
    }else{
      gkeep <- pfre[[1]][1:gs]
    }

    pdata <- cvlist[[2]]

    if(ytype == 'percentage'){
      pdata <- pdata[,grepl('gene$|%$', colnames(pdata), ignore.case = T)]
      colnames(pdata) <- gsub('%$', '', colnames(pdata))
      pdata[,2:ncol(pdata)] <- sapply(pdata[,2:ncol(pdata)], function(x){x/100})
    }else if(ytype == 'counts'){
      pdata <- pdata[,!(grepl('%$', colnames(pdata)))]
    }
    pdata <- pdata[pdata[[1]] %in% gkeep,]
    pdata[['Gene']] <- factor(pdata[['Gene']], levels = gkeep)
    pdata <- tidyr::gather(pdata, key = 'Variants', value = 'y', -1)
    pdata[['Variants']] <- factor(pdata[['Variants']], levels = names(color))
    return(pdata)
  }

  if(is.null(by)){
    pdata <- topdata(oridata)
    ggcolor <- color[names(color) %in% unique(pdata[[2]])]
    p <- ori_ggbar(pdata, 'Gene', 'y', 'Variants', ggcolor)
    if(!(is.null(outpdf))){
      ggsave(plot = p, filename = outpdf, height = 4, width = 6)
    }
  }else{
    #
    pdata.whole <- topdata(oridata[[1]])
    ggcolor <- color[names(color) %in% unique(pdata.whole[[2]])]
    p.whole <- ori_ggbar(pdata.whole, 'Gene', 'y', 'Variants', ggcolor)

    #bar by group
    pdata.g1 <- topdata(oridata[[byorder[1]]], pdata.whole[['Gene']])
    pdata.g2 <- topdata(oridata[[byorder[2]]], pdata.whole[['Gene']])
    pdata.g2[['y']] <- -1*pdata.g2[['y']]
    pdata.g <- rbind(pdata.g1, pdata.g2)
    p <- ori_ggbar(pdata.g, 'Gene', 'y', 'Variants', ggcolor)

    p <- ggpubr::ggarrange(p.whole, p, nrow = 2, legend = 'bottom', common.legend = T)
    if(!(is.null(outpdf))){
      ggsave(plot = p, filename = outpdf, height = 9, width = 6)
    }
  }
  return(p)
}
