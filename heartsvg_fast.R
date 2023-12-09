heartsvg_fast=function (data, scale = T, qh = 0.985, noise = T, padj.m = "holm") {
  
  scale.count=function(data,qh=0.985){
    counts_mat=apply(data, 1, function(y){z=y/max(y)})  
    ## 类似于调节测序深度
    counts_mat=apply(counts_mat, 2, function(y)
    {y[which(y>quantile(y,qh))]<-quantile(y,qh);
    y[which(y<quantile(y,0.25))]<-0;y})
    t(counts_mat)
  }
  
  
  new <- if (scale) {
    cbind(data[, 1:2], scale.count(data[, -c(1:2)]))
  } else {
    data
  }
  rm(data)
  gc()
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
  rm(new)
  gc()
  new = as.data.table(cbind(locus_in1[colnames(locus_in1) != 
                                        l1], locus_in2[colnames(locus_in2) != l2], counts))
  colnames(new)[c(2, 4)] <- c("v1", "v2")
  a2 = Sys.time()
  new.row <- new[, lapply(.SD, mean), by = .(col, v1), .SDcols = names(new)[-c(1:4)]]
  new.col <- new[, lapply(.SD, mean), by = .(row, v2), .SDcols = names(new)[-c(1:4)]]
  b2 = difftime(Sys.time(), a2, units = "mins")
  cat("new--", b2)
  zero.p <- colSums(new[, -c(1:4)] != 0, na.rm = TRUE)/nrow(new[, 
                                                                -c(1:4)])
  mean <- colMeans(new[, -c(1:4)][, lapply(.SD, function(y) mean(y[y != 
                                                                     0], na.rm = TRUE))])
  sum <- data.frame(gene = colnames(new[, -c(1:4)]), zero.p, 
                    mean)
  row.t = apply(new.row[, -c(1:2)], 2, function(y) {
    z = ifelse(sum(y != 0) == 0, 1, Box.test(y, lag = z1.group)$p.value)
    z
  })
  col.t = apply(new.col[, -c(1:2)], 2, function(y) {
    z = ifelse(sum(y != 0) == 0, 1, Box.test(y, lag = z2.group)$p.value)
    z
  })
  new.x = new[, lapply(.SD, mean), by = .(row), .SDcols = names(new)[-c(1:4)]]
  new.y = new[, lapply(.SD, mean), by = .(col), .SDcols = names(new)[-c(1:4)]]
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
