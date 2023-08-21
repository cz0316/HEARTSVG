#' @title  heartsvg
#' @description Find spatially variable genes (SVGs).
#' @param data A data.frame with dimensions n*(p+2) contains gene coordinates and
#'gene expression counts. The first 2 columns should represent coordinates and
#'their column names must be 'row' and 'col', respectively.
#'The remaining p columns correspond to the expression levels of p genes.
#'There are n rows, representing n spots.
#' @param scale A logic value: T ot F. The default value is 'T'.
#' If scale==T, the expression counts matrix runs function scale.count to scaling expression.
#' @param noise  A logic value: T or F. The default value is 'T'.
#' If noise==F, the SVG-rank considers only the p-value and adjusted p-value.
#' If noise==T, the SVG-rank takes into account the p-value,adjusted p-valus
#'  and the distribution characteristics of the genes.
#' @param padj.m A character string. The default value is 'holm'.
#' The adjusted method of p-values.
#' The adjustment methods include "bonferroni", "holm", "hochberg",
#' "hommel","BH" or its alias "fdr", and "BY".
#' @return A data.frame with 4 columns, includes gene name("gene"),
#' p-value("pval"), adjusted p-values("p_adj"), and SVG-rank("rank").
#' @export
heartsvg=function (data, scale = T,
                   qh=0.985,noise = T, padj.m = "holm")
{
  if (scale == T) {
    new = cbind(data[, 1:2], scale.count(data[, -c(1:2)]))
  }
  if (scale == F) {
    new = data
  }
  new = as.data.frame(new)
  locus_in = new[, c(1, 2)]
  counts = new[, -c(1, 2)]
  l1 = "row"
  l2 = "col"
  z1.group = ceiling(log(diff(range(locus_in[l1]))))
  z2.group = ceiling(log(diff(range(locus_in[l2]))))
  locus_in1 = data.frame(locus_in, n.row = cut(x = as.matrix(locus_in[l1]),
                                               breaks = seq(from = min(locus_in[l1]), to = max(locus_in[l1]),
                                                            length = z1.group + 1), labels = 1:z1.group, include.lowest = T,
                                               right = T))
  locus_in2 = data.frame(locus_in, n.col = cut(x = as.matrix(locus_in[l2]),
                                               breaks = seq(from = min(locus_in[l2]), to = max(locus_in[l2]),
                                                            length = z2.group + 1), labels = 1:z2.group, include.lowest = T,
                                               right = T))
  new1 = cbind(locus_in1[colnames(locus_in1) != l1], counts)
  new2 = cbind(locus_in2[colnames(locus_in2) != l2], counts)
  colnames(new1)[1:2] <- c("x", "coor.z")
  colnames(new2)[1:2] <- c("x", "coor.z")
  new.row = aggregate(new1[, -c(1:2)], list(new1$x, new1$coor.z),
                      mean)
  new.col = aggregate(new2[, -c(1:2)], list(new2$x, new2$coor.z),
                      mean)
  zero.p = apply(new[, -c(1:2)], 2, function(y) {
    sum(y != 0, na.rm = T)/length(y)
  })
  mean = apply(new[, -c(1:2)], 2, function(y) {
    mean(y[y != 0], na.rm = T)
  })
  sum = data.frame(gene = names(mean), zero.p, mean)
  row.t = apply(new.row[, -c(1:2)], 2, function(y) {
    z = ifelse(sum(y != 0) == 0, 1, Box.test(y, lag = z1.group)$p.value)
    z
  })
  col.t = apply(new.col[, -c(1:2)], 2, function(y) {
    z = ifelse(sum(y != 0) == 0, 1, Box.test(y, lag = z2.group)$p.value)
    z
  })
  new.x = aggregate(new[, -c(1:2)], list(new$row), mean) ##  其实是实际上的列方向 by col
  new.y = aggregate(new[, -c(1:2)], list(new$col), mean)
  x.t = apply(new.x[, -1], 2, function(y) {
    Box.test(y, lag = z1.group)$p.value
  })
  y.t = apply(new.y[, -1], 2, function(y) {
    Box.test(y, lag = z2.group)$p.value
  })
  test = data.frame(row.t, col.t, x.t, y.t, gene = names(row.t))
  test[, 1:(ncol(test) - 1)] = sapply(test[, 1:(ncol(test) -
                                                  1)], function(y) {
                                                    z = ifelse(is.na(y) == T, 0.999, y)
                                                  })
  test$min = apply(test[, 1:(ncol(test) - 1)], 1, function(y) {
    poolr::stouffer(y)$p
  })
  test$adj_min = p.adjust(test$min, method = padj.m)
  test = merge(test, sum, by = "gene")
  test = subset(test, zero.p > 0)
  test$zero.p = test$zero.p/max(test$zero.p)
  test$mean = test$mean/max(test$mean)
  test$c2 = (test$zero.p + test$mean)/2
  if (noise == T) {
    data.table::setorder(test, adj_min, min, -zero.p, -mean)
    a1 = c("gene", "min", "adj_min")
    test2 = test[a1]
    b1 = c("gene", "pval", "p_adj")
    colnames(test2) <- b1
    test2$rank = 1:nrow(test2)
  }
  if (noise == F) {
    data.table::setorder(test, adj_min, min)
    a1 = c("gene", "min", "adj_min")
    test2 = test[a1]
    b1 = c("gene", "pval", "p_adj")
    colnames(test2) <- b1
    test2$rank = 1:nrow(test2)
  }
  test2
}
