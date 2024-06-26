% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heartsvg.R
\name{heartsvg}
\alias{heartsvg}
\title{heartsvg}
\usage{
heartsvg(data, scale = T, qh = 0.985, noise = T, padj.m = "holm")
}
\arguments{
\item{data}{A data.frame with dimensions n*(p+2) contains gene coordinates and
gene expression counts. The first 2 columns should represent coordinates and
their column names must be 'row' and 'col', respectively.
The remaining p columns correspond to the expression levels of p genes.
There are n rows, representing n spots.}

\item{scale}{A logic value: T ot F. The default value is 'T'.
If scale==T, the expression counts matrix runs function scale.count to scaling expression.}

\item{noise}{A logic value: T or F. The default value is 'T'.
If noise==F, the SVG-rank considers only the p-value and adjusted p-value.
If noise==T, the SVG-rank takes into account the p-value,adjusted p-valus
 and the distribution characteristics of the genes.}

\item{padj.m}{A character string. The default value is 'holm'.
The adjusted method of p-values.
The adjustment methods include "bonferroni", "holm", "hochberg",
"hommel","BH" or its alias "fdr", and "BY".}
}
\value{
A data.frame with 4 columns, includes gene name("gene"),
p-value("pval"), adjusted p-values("p_adj"), and SVG-rank("rank").
}
\description{
Find spatially variable genes (SVGs).
}
