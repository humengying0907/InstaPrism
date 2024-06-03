commonRows=function(data1, data2){
  intersect(rownames(data1), rownames(data2))
}

# function to normalize expression matrix, s.t. it sum up to one for each row, with the zero entries = pseudo.min ((adapted from BayesPrism))
norm.to.one <- function(ref,
                        pseudo.min){

  G <- ncol(ref)

  phi <- ref/rowSums(ref) * (1-pseudo.min*G) + pseudo.min

  #if the minimum value is greater than zero. simply normalize by total depth
  min.value <- apply(ref,1,min)
  which.row <- min.value>0
  if(any(which.row)){
    phi[which.row,] <- ref[which.row,,drop=F]/rowSums(ref[which.row,,drop=F])
  }


  return(phi)

}

#' function to filter bulk outliers (adapted from BayesPrism)
#' @param mixture the bulk RNA-seq matrix (#of samples * # of genes).
#' @param outlier.cut & outlier.fraction: Filter genes in mixture whose expression fraction is greater than outlier.cut
#'		in more than outlier.fraction of bulk data. Removal of outlier genes will ensure that the inference will not be dominated by outliers.
#' @keywords internal
#' @noRd
filter.bulk.outlier <- function(mixture,
                                outlier.cut,
                                outlier.fraction){

  mixture.norm <- mixture / rowSums(mixture)

  outlier.idx <- colSums(mixture.norm > outlier.cut) / nrow(mixture.norm) > outlier.fraction
  mixture <- mixture[, !outlier.idx, drop=F]
  # cat("Number of outlier genes filtered from mixture =", sum(outlier.idx),"\n")

  return(mixture)
}

#' Get sub-clustering labels
#' @description Given single-cell expression and cell-type labels, return sub-clustering labels for each cell type. Currently supported
#'    subclustering-methods include scran::quickCluster() and SCISSORS::ReclusterCells().
#'
#' @param scExpr single cell expression matrix with genes in rows and cells in columns
#' @param cell_type_labels a vector indicating cell-type level annotations
#' @param subcluster_method subclustering methods to use. Available options include 'scran' and 'SCISSORS'. Default = 'scran'
#' @param min.subcluster.size min.size parameter passed to scran::quickCluster(). Default = 100
#' @param CalculateSilhouette a logical variable determine whether Silhouette score need to be calculated to determine which cell-types to be reclustered. Default = T. When set to FALSE, will
#'    recluster every cell type present in 'cell_type_labels'
#' @param SilhouetteScores.cutoff a cutoff determining which cell-types to be reclustered with 'SCISSORS' method. Cell-types with SilhouetteScores greater than this
#'    cutoff will not be reclustered
#' @param ... additional parameters pass to SCISSORS::ReclusterCells()
#'
#' @return a vector of subclustering labels, which can be pass to `refPrepare()` function for reference construction
#' @export
#'
get_subcluster = function(scExpr,cell_type_labels,
                          subcluster_method = 'scran',
                          min.subcluster.size = 100,
                          CalculateSilhouette = T,
                          SilhouetteScores.cutoff = 0.7,...){

  cell_type_labels = as.vector(cell_type_labels)

  if (subcluster_method == 'scran'){
    if(!requireNamespace("scran", quietly = TRUE)){
      stop("The 'scran' package is not installed. Please install the scran package using the following command and try again: \n",'BiocManager::install("scran")')
    }else{
      library(scran)
    }

    stopifnot(ncol(scExpr)==length(cell_type_labels))

    group = list()
    for(i in unique(cell_type_labels)){
      group[[i]] <- which(cell_type_labels %in% i)
    }

    get_sublabels = function(ct){
      cl = scran::quickCluster(scExpr[,group[[ct]]],min.size = min(min.subcluster.size,length(group[[ct]])))
      labels = paste0(ct,'_',cl)
      return(labels)
    }

    f = lapply(names(group),get_sublabels)
    names(f) = names(group)

    subcluster_IDtable = data.frame(id = do.call(c,group),
                                    label = do.call(c,f))

    subcluster_IDtable = subcluster_IDtable[order(subcluster_IDtable$id),]
    return(subcluster_IDtable$label)
  }else if(subcluster_method == 'SCISSORS'){

    if(!requireNamespace("SCISSORS", quietly = TRUE)){
      stop("The 'SCISSORS' package is not installed. Please install the SCISSORS package using the following command and try again: \n",'remotes::install_github("jr-leary7/SCISSORS")')
    }else{
      library(SCISSORS)
    }
    require(Seurat)
    require(SCISSORS)
    seurat_obj = CreateSeuratObject(counts = scExpr)
    seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
    seurat_obj <- FindVariableFeatures(seurat_obj,selection.method = "vst", nfeatures = 5000)
    seurat_obj <- ScaleData(seurat_obj)
    seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
    seurat_obj <- FindNeighbors(seurat_obj)
    seurat_obj$cell_type = cell_type_labels
    seurat_obj$seurat_clusters = as.factor(as.numeric(as.factor(cell_type_labels)))
    Idents(seurat_obj) = seurat_obj$seurat_clusters
    seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

    annot = levels(as.factor(cell_type_labels))

    if(length(annot) ==1){
      warning('only one cluster provided, SilhouetteScores is not available, will call ReclusterCells() on the entire scExpr')
      message(paste('recluster',annot))
      reclust_obj  = SCISSORS::ReclusterCells(seurat_obj,
                                              which.clust = 1,
                                              merge.clusters = F,
                                              n.PC = 15,
                                              use.parallel = FALSE,
                                              redo.embedding = T,...)
      subcluster_label = paste0(annot,'_',reclust_obj$seurat_clusters)
      gc()
      return(subcluster_label)
    }else{
      if(CalculateSilhouette){
        # decide which clusters to recluster
        SilhouetteScores = ComputeSilhouetteScores(seurat_obj,avg = T) %>% as.numeric()
        recluster_id = which(SilhouetteScores < SilhouetteScores.cutoff)

        gc()

        if(length(recluster_id)<1){
          message('no cell-type need to be reclustered under the given SilhouetteScores.cutoff, will return original cell type labels')
          return(cell_type_labels)
        }

      }else{
        recluster_id = seq(1,length(unique(cell_type_labels)))
      }

      seurat_obj$subcluster = seurat_obj$cell_type

      for(i in recluster_id){
        message(paste('recluster',annot[i]))
        reclust_obj  = SCISSORS::ReclusterCells(seurat_obj,
                                                which.clust = i,
                                                merge.clusters = F,
                                                n.PC = 15,
                                                use.parallel = FALSE,
                                                redo.embedding = T,...)
        id = match(colnames(reclust_obj),colnames(seurat_obj))
        subcluster_label = paste0(annot[i],'_',reclust_obj$seurat_clusters)

        seurat_obj$subcluster[id] = subcluster_label
        gc()

      }
      return(seurat_obj$subcluster)
    }
  }else{
    stop('please provide a valid subcluster method')
  }
}
