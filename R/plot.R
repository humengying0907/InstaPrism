
#' Compare ground truth and estimated cell type fractions
#' @description Calculate Pearson correlation between the ground truth and the estimated cell type fractions for each cell type.
#'    When cell type names in the two fraction matrices do not match exactly, this function automatically identifies the most likely
#'    correspondences between cell types by finding pairs with the highest pairwise correlation
#' @param truth_frac ground truth cell-type fraction with samples in rows and cell-types in columns
#' @param estimated_frac cell-type fraction estimation with samples in rows and cell-types in columns
#'
#' @return a list with summary statistics of the correlation analysis
#' @export
#'
frac_evalu = function(truth_frac,estimated_frac){

  if(nrow(truth_frac)!=nrow(estimated_frac)){
    stop('truth_frac and testimated_frac must correspond to the same samples in their rows.')
  }

  if(any(colSums(truth_frac)==0) | any(colSums(estimated_frac)==0)){
    stop('Please exclude any cell types that have zero fractions across all samples.')
  }

  if(is.null(colnames(truth_frac))|is.null(colnames(estimated_frac))){
    stop('The provided fraction matrix must include cell type names as its column.')
  }

  if(!requireNamespace("tibble", quietly = TRUE)){
    stop('The "tibble" package is not installed. Please install it using install.packages("tibble") and try again')
  }else{
    library(tibble)
  }

  Y = truth_frac
  E = estimated_frac

  Y = Y[,order(colnames(Y)),drop = F]
  E = E[,order(colnames(E)),drop = F]

  if(ncol(Y) == ncol(E)){
    if(all(colnames(Y) == colnames(E))){
      E = E
    }else{
      # find cell types that exhibit the highest correlation with true_frac in Y
      maxCorName = c()
      for(ct in colnames(Y)){
        id = apply(cor(E,Y[,ct],use = 'pairwise.complete.obs'),2,which.max)
        maxCorName = c(maxCorName, colnames(E)[id])
      }

      E = E[,maxCorName,drop = F]
      colnames(E) = make.unique(colnames(E))

    }
  }else{
    maxCorName = c()
    for(ct in colnames(Y)){
      id = apply(cor(E,Y[,ct],use = 'pairwise.complete.obs'),2,which.max)
      maxCorName = c(maxCorName, colnames(E)[id])
    }

    E = E[,maxCorName,drop = F]
    colnames(E) = make.unique(colnames(E))
  }

  m1 = gather(Y %>% as.data.frame(),cell_type,true_frac)
  m2 = gather(E %>% as.data.frame(),maxCorName,estimate)
  M = cbind(m1,m2)

  # add summary statistics
  summ <- M %>%
    group_by(cell_type) %>%
    summarise(
      RMSE = caret::RMSE(true_frac, estimate,na.rm = T),
      cor = cor(true_frac,estimate,method = 'pearson',use = 'pairwise.complete.obs')) %>%
    mutate_if(is.numeric, round, digits=2) %>% as.data.frame() %>% column_to_rownames('cell_type')

  summ$maxCorName = M$maxCorName[match(rownames(summ),M$cell_type)]

  return(list(M=M,
              summ = summ))
}


#' Visualize deconvolution performance
#' @description Visualize per-cell type correlation between true fraction and estimates
#' @param truth_frac ground truth cell-type fractions with samples in rows and cell-types in columns
#' @param estimated_frac cell-type fraction estimates with samples in rows and cell-types in columns
#' @param title title for the plot
#' @param nrow number of rows in the plot grid. Default = 1
#' @param xlabel xlabel of the plot. Default = 'true fraction'
#' @param ylabel ylabel of the plot. Default = 'estimates'
#' @export
#'
deconv_performance_plot<-function(truth_frac,estimated_frac,title=NULL,nrow=1,xlabel='true fraction',ylabel='estimates'){
  require(ggplot2)
  require(ggpmisc)
  require(tibble)

  if(nrow(truth_frac)!=nrow(estimated_frac)){
    stop('truth_frac and testimated_frac must correspond to the same samples in their rows.')
  }

  if(any(colSums(truth_frac)==0) | any(colSums(estimated_frac)==0)){
    stop('Please exclude any cell types that have zero fractions across all samples.')
  }

  if(is.null(colnames(truth_frac))|is.null(colnames(estimated_frac))){
    stop('The provided fraction matrix must include cell type names as its column.')
  }

  l = frac_evalu(truth_frac,estimated_frac)

  M=l$M
  summ=l$summ
  summ = summ %>% rownames_to_column('cell_type')
  summ = summ[order(summ$cell_type),] # order alphabetically

  p=ggplot(M,aes(true_frac,estimate))+
    geom_point()+
    geom_abline(intercept = 0, slope = 1,linetype='dotted',color='red')+
    facet_wrap(~cell_type,nrow=nrow)+
    theme(aspect.ratio=1)+
    ggtitle(title)+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))+
    xlab(xlabel)+
    ylab(ylabel)

  p + geom_table_npc(data = summ, label = lapply(split(summ, summ$cell_type),
                                                 FUN = function(entry) {subset(entry, select = -cell_type)}),
                     npcx = 0.00, npcy = 1, hjust = 0, vjust = 1, size=3,
                     table.theme = ttheme_gtlight)


}



