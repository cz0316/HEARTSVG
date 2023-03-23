#' @title  svg.cluster
#'
#' @description  Predictions of spatial domains based on SVG clustering results.
#' Cluster analysis on SVGs by the spatial expression similarity.
#'
#' @param data A data.frame with dimensions n*(p+2) contains gene coordinates and
#'gene expression counts. The first 2 columns should represent coordinates and
#'their column names must be 'row' and 'col', respectively.
#'The remaining p columns correspond to the expression levels of p genes.
#'There are n rows, representing n spots.
#'It is the same as the parameter 'data' in the function 'heartsvg()'.

#' @param svg A vector,indicating columns(genes) to select from the ST data.

#' @param method A character string for clustering method.
#' 'h' represents hierarchical cluster analysis.
#' 'k' represents k-means cluster analysis. If k-means method is chosen,
#' the parameter 'n.c' have to be provided.
#'
#' @param n.c A numerical value of the number of clusters.
#' For hierarchical clustering, n.c is a dispensable parameter.
#' If n.c is not provided in hierarchical clustering, our approach adopt
#' Yamamoto test to estimate the number of SVG cluters.
#' For k-means method, n.c must be provided.
#'
#' @return A data.frame with 2 columns, the first column ('gene') is the gene name.
#' The second column ('cluster') is the cluster of each gene.
#'
#' @export


svg.clust=function(data,svg,method='h',n.c=NULL){
  ## method =h/k
  coor=data[,1:2]
  ct=as.data.frame(simu.scale(data[,-c(1:2)]))
  counts=subset(ct,select =svg)

  if (method=='h'){
    co=dist(t(counts))
    # co[co=='NA']<-0
    hh=hclust(as.dist(co),method = 'complete')
    bp.od=yamamoto(hh$height)$breaks[length(yamamoto(hh$height)$breaks)]
    if ( is.null(n.c)){  ##未指定类的个数
      c6=cutree(hh,h=hh$height[bp.od])
      svg.cl=data.frame(gene=hh$labels,cl=c6)
    }
    if(!is.null(n.c)){   ##指定类的个数 =n.c
      c6=cutree(hh,k=n.c)
      svg.cl=data.frame(gene=hh$labels,cl=c6)
    }
  }
  if (method=='k'){
    if ( is.null(n.c)){  ##未指定类的个数
      stop('Missing the parameter n.c !')
      # quit()
    } else{
      km=kmeans(t(counts),centers = n.c)  ##kmeans默认对行聚类，所以要t()
      svg.cl=data.frame(gene=names(km$cluster),cl=km$cluster)
    }
  }
  colnames(svg.cl)<-c('gene','cluster')
  svg.cl

}
