
#' @title  Enrichment analysis for SVGs.
#' @description  This function is based on the package gprofiler2
#'  for gene enrichment analysis.
#' @param svg  A vector or data frame containing gene names for enrichment analysis.
#' @param plot A logical value indicating whether to plot the
#' graphical results of the enrichment analysis. The default is True.
#' @param meta A logical value indicating whether to return the complete gprofiler2 results,
#' including annotations, background genes, etc. The default is False.
#' @param organism  Organism name. Default value is ‘hspainens'
#'  More details in the R package gprofiler2 or the g:GOSt web tool
#' @param sources A vector of data sources to use.
#' Currently, these include GO (GO:BP, GO:MF, GO:CC to select a particular GO branch),
#' KEGG, REAC, TF, MIRNA, CORUM, HP, HPA, WP.
#' Please see the R package gprofiler2 or the g:GOSt web tool for
#' the comprehensive list and details on incorporated data sources
#' @param user_threshold A numeric value.
#' The p-value threshold for significance. Default value is 0.05.
#' @param correction_method The algorithm used for multiple testing correction,
#' one of "gSCS" (synonyms: "analytical", "g_SCS"),"fdr" (synonyms: "false_discovery_rate"), "bonferroni".
#' @return A data.frame contains the enrichment analysis results.

#' @export
svg.enrichment=function(svg,plot=T,meta=F,organism = 'hsapiens',sources = NULL,user_threshold = 0.05,  correction_method = "g_SCS"){
  ##  data: a data.frame must have at least
  ##  two columns (gene and cl corresponding to cluster)
  ##  OR the result of svg.cluster,
  ##  organism:  organism name. Default value is ‘hspainens'
  ##  more details in the R package gprofiler2 or the g:GOSt web tool
  ##  sources: a vector of data sources to use.
  ##  Currently, these include GO (GO:BP, GO:MF, GO:CC to select a particular GO branch), KEGG, REAC, TF, MIRNA, CORUM, HP, HPA, WP.
  ##  Please see the g:GOSt web tool for the comprehensive list and details on incorporated data sources
  ##  user_threshold : A numeric value. The p-value threshold for significance.
  ##  Default value is 0.05.
  ##  correction_method:   the algorithm used for multiple testing correction,
  ##  one of "gSCS" (synonyms: "analytical", "g_SCS"),
  ##  "fdr" (synonyms: "false_discovery_rate"), "bonferroni"
  # svg=data['gene'][data['cl']==cluster]
  g1=gost(query=svg,
          organism = organism,sources = sources,
          user_threshold = user_threshold,  correction_method =correction_method,
          ordered_query = FALSE,
          multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
          measure_underrepresentation = FALSE, evcodes = FALSE,
          domain_scope = "annotated", custom_bg = NULL,
          numeric_ns = "",  as_short_link = FALSE)


  if(plot==T){
    print(gostplot(g1))
  }

  if(meta==T){result=g1}
  if(meta==F){result=g1$result}
  result

}
