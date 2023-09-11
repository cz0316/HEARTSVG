#' @title  svg.cluster
#' @description  Predictions of spatial domains based on SVG clustering results.
#' Cluster analysis on SVGs by the spatial expression similarity.
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
#' @param n.c A numerical value of the number of clusters.
#' For hierarchical clustering, n.c is a dispensable parameter.
#' If n.c is not provided in hierarchical clustering, our approach adopt
#' Yamamoto test to estimate the number of SVG cluters.
#' For k-means method, n.c must be provided.
#' @return A data.frame with 2 columns, the first column ('gene') is the gene name.
#' The second column ('cluster') is the cluster of each gene.
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

#' Title svg.seek
#'
#' @param data   A data.frame with dimensions n*(p+2) contains gene coordinates and
#'gene expression counts. The first 2 columns should represent coordinates and
#'their column names must be 'row' and 'col', respectively.
#' @param target A character string of the target gene for which you want to find nearby genes.
#' @param svg  A vector. This is a gene list used to narrow down the search range. 
#' The function will look for nearby genes within this list.
#' @param size This is an optional parameter that specifies how many nearby genes you want to return. 
#' If provided, the function will return the top size genes closest to the target gene. 
#' If this parameter is not provided, the function will determine the number of nearby genes
#'  based on some breakpoint 
#' @returnA data.frame with 2 columns, the first column ('gene') is the gene name.
#' The second column ('cluster') is the Euclidean distance.
#' @export
svg.seek=function(data,target,svg,size=NULL){
  ## size不是必需的
  ## svg (一个基因列表，用于缩小寻找范围)
  ## gene目标基因，一个
  coor=data[,1:2]
  ct=as.data.frame(simu.scale(data[,-c(1:2)]))
  
  # svg=str_replace(svg,'-','.')
  counts=subset(ct,select =c(target,intersect(svg,colnames(ct))))
  # svg2=str_replace(svg,'\\.','-')
  
  y=subset(counts,select=target)
  
  x=counts[colnames(counts)!=target]
  
  dis.seek=apply(x, 2, function(z){sqrt(sum((z-y)^2))})
  ll=dis.seek[order(dis.seek)]
  
  if (!is.null(size)){
    near.gene=data.frame(gene=names(ll[1:size]),dis=ll[1:size])
  }
  if (is.null(size)){
    
    bp=yamamoto(dis.seek)$breaks[1]
    
    near.gene=data.frame(gene=names(ll[1:bp]),dis=ll[1:bp])
  }
  
  dt=data.frame(neighbors=near.gene$gene,dist=near.gene$dis)
  return(dt)
}


