#' @include ori_tools.R
#' @import tidyr
#' @import openxlsx
NULL

# var <- read.delim('../Data/data/variation_somatic.txt', stringsAsFactors = F)
#' Count mutation frequency of gene by ID

#' @param xlsxfile An xlsx file which the data are to be read from.
#' @param id Column name of ID.
#' @param gene Column name of Gene.
#' @param vartype Column name of Variant type, eg: VAR_TYPE_SX.
#' @param varorder Priority of variants, eg:c('Fusion/Rearrangement', 'Substitution/Indel', 'Gene Amplification', 'Gene Homozygous Deletion', 'Truncation')
#' @param by Column names which used to split data into groups.
#' @param outxlsx xlsx file name. If value is 'list', a list object will be return.
#'
#' @return A list of table
#' @examples
#' CountVariants(data, id = 'ORDER_ID', gene = 'GENE', vartype = 'VarType')
#'
#' @export

CountVariants <- function(xlsxfile,
                          id = 'ID',
                          gene = 'Gene',
                          vartype = 'VarType',
                          varorder = c('Fusion/Rearrangement', 'Substitution/Indel', 'Gene Amplification', 'Gene Homozygous Deletion', 'Truncation'),
                          by = NULL,
                          outxlsx = "frequency.xlsx"
                          ){
  if(!(is.vector(varorder))){
    warning('No Variant priority provided, and default will be used')
    varorder = c('Fusion/Rearrangement', 'Substitution/Indel', 'Gene Amplification', 'Gene Homozygous Deletion', 'Truncation')
  }

  #oridata <- openxlsx::read.xlsx(xlsxfile, sheet = 1, skipEmptyRows = TRUE,
  #                           na.strings = c("NA", "NaN"), detectDates = TRUE)
  oridata <- xlsxfile

  if(is.null(by)){
    oridata <- oridata[,c(id, gene, vartype)]
  }else{
    oridata <- oridata[,c(id, gene, vartype, by)]
  }

  oridata <- rmvar(oridata, varorder)

  keggpath <-readRDS('./data/keggpath.rds')

  #
  CVsub <- function(CVdata){

    uid <- unique(CVdata[[1]])
    colnames(CVdata) <- c('ID', 'Gene', 'Var')

    # Gene mutation frequency
    CVdata.genefre <- aggregate(ID~Gene, CVdata, length)
    CVdata.genefre <- as.data.frame(CVdata.genefre, stringAsFactors = F)
    colnames(CVdata.genefre) <- c('Gene', 'Frequency')
    CVdata.genefre$TotalNumber <- length(uid)
    CVdata.genefre$Percentage <- round(CVdata.genefre$Frequency/length(uid) * 100, digits = 2)
    CVdata.genefre <- CVdata.genefre[order(CVdata.genefre$Frequency, decreasing = T),]

    # Variants mutation frequency
    CVdata.varfre <- aggregate(ID~Gene + Var, CVdata, length)
    CVdata.varfre <- tidyr::spread(CVdata.varfre, key = 'Var', value = 'ID')
    CVdata.varfre[is.na(CVdata.varfre)] <- 0
    CVdata.varfre <- as.data.frame(CVdata.varfre, stringAsFactors = F)
    CVdata.varfre <- CVdata.varfre[,c('Gene', varorder[varorder %in% unique(CVdata[[3]])])]
    CVdata.varfre <- CVdata.varfre[match(CVdata.genefre$Gene, CVdata.varfre$Gene),]
    vcol <- colnames(CVdata.varfre)[2:ncol(CVdata.varfre)]
    CVdata.varfre[,paste0(vcol, '%')] <- sapply(CVdata.varfre[,vcol], function(x){round(x/length(uid) * 100, digits = 2)})

    # Pathway alteration frequency
    CVdata.pathfre <- lapply(names(keggpath$kg.genesymbolsets),
                           function(x){
                             xgene <- keggpath$kg.genesymbolsets[[x]];
                             dmut <- CVdata[CVdata[[2]] %in% xgene,]
                             mutg <- paste0(unique(dmut[[2]]), collapse = ',')
                             mutid <- length(unique(dmut[[1]]));
                             return(c(x, paste0(xgene, collapse = ','), mutg, mutid))})
    CVdata.pathfre <- do.call(rbind, CVdata.pathfre)
    CVdata.pathfre <- as.data.frame(CVdata.pathfre, stringsAsFactors = F)
    colnames(CVdata.pathfre) <- c('Pathway', 'Items', 'Mutant items', 'Frequency')
    CVdata.pathfre$Frequency <- as.numeric(CVdata.pathfre$Frequency)
    CVdata.pathfre <- CVdata.pathfre[order(CVdata.pathfre$Frequency, decreasing = T),]
    CVdata.pathfre$Percentage <- round(CVdata.pathfre$Frequency/length(uid) * 100, digits = 2)

    return(list(GeneFre = CVdata.genefre, VarFre = CVdata.varfre, PathwayFre = CVdata.pathfre))
  }

  #
  if(is.null(by)){
    oridata.fre <- CVsub(oridata[,1:3])
  }else{
    ugg <- table(oridata[[4]])
    ugg <- names(sort(ugg, decreasing = T))
    ugg <- gsub('\\s+', '_', ugg)
    oridata.freAll <- CVsub(oridata[,1:3])
    oridata.fre <- sapply(ugg,
                       function(x){
                         xfre <- CVsub(oridata[oridata[[4]] == x, 1:3]);
                         names(xfre) <- paste(x, names(xfre), sep = '_');
                         return(xfre)}, simplify = F, USE.NAMES = T)
    oridata.fre <- c(list(WholeDataset = oridata.freAll), oridata.fre)
    if(outxlsx != 'list') oridata.fre <- Reduce(c, oridata.fre)
  }
  return(oridata.fre)
}
