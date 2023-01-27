

#' Title
#'
#' @param est a
#' @param true a
#' @param title a
#' @param nrow a
#' @export

deconv_performance_plot<-function(est,true,title=NULL,nrow=2,xlabel='true fraction',ylabel='estimates'){
  est=t(est)
  true=t(true)

  est=est[,colnames(true)]
  m1=gather(true %>% as.data.frame(),cell_type,true_frac)
  m2=gather(est %>% as.data.frame(),maxCorName,estimate)
  M=cbind(m1,m2)

  summ <- M %>%
    group_by(cell_type) %>% # note that in summ, cell_type is ordered alphabetically
    summarise(RMSE = caret::RMSE(true_frac, estimate,na.rm = T),
              cor = cor(true_frac,estimate,method = 'pearson',use = 'pairwise.complete.obs')) %>%
    mutate_if(is.numeric, round, digits=2)

  summ=summ[order(summ$cell_type),]

  p=ggplot(M,aes(true_frac,estimate))+
    geom_point()+
    geom_abline(intercept = 0, slope = 1,linetype='dotted',color='red')+
    facet_wrap(~cell_type,nrow=nrow)+  # scatter plot will be arranged alphabetically
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
