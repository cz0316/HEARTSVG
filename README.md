# HEARTSVG
HEARTSVG: a high-efficiency and robust method for identifying spatial expression patterns in large-scale spatial transcriptomic data
We proposed HEARTSVG, which builds upon serial correlation tests for rapidly and precisely detecting SVGs without relying on the underlying data-generative model or pre-defining spatial patterns.Furthermore, we developed the software HEARTSVG to provide an SVG auto-clustering module for predictions of spatial domain and visualizations, allowing users to discover novel biological phenomena. 

To be able to use the R package HEARTSVG successfully, we suggest manually installing the following packages and their dependent packages first: spatstat, Seurat, plotly, plot3D, TTR, stringr, fpp2, seasonal, cluster, poolr, reshape2, data.table, plyr, dplyr, ggplot2, MASS, EnvStats, BreakPoints, gprofiler2

```{r,warning=F,message=F,results='hide',error=F}
library(spatstat);library(Seurat);library(plot3D);library(plotly);library(TTR);
library(stringr);library(fpp2);library(seasonal);library(cluster);library(poolr);
library(reshape2);library(data.table);library(plyr);library(dplyr);library(ggplot2);
library(MASS);library(EnvStats);library(BreakPoints);library(gprofiler2)

```
## Required R modules
```{r ,warning=F}
R (>= 3.5.0)
```

## Dependencies
```{r ,warning=F}
fpp2 (>=2.4)
seasonal (>=1.9.0)
cluster (>=2.1.3)
poolr (>=1.1.1)
Seurat (4.1.1)
spatstat (>=2.3.4)
ggplot2 (>=3.4.2)
EnvStats (>=2.7.0)
BreakPoints (>=1.2)
plotly (>=4.10.0)
plot3D (>=1.4)
gprofiler2 (>=0.2.1)


```

## Install

```{r}
library(devtools)
devtools::install_github('https://github.com/cz0316/HEARTSVG',force = T)
library(HEARTSVG)
```
## 0 Load data

We use the ST data from the Wu et al. study as an example to introduce the usage of the R package HEARTSVG. The article and data download link are: https://doi.org/10.1158/2159-8290.CD-21-0316; http://www.cancerdiversity.asia/scCRLM/.


```{r ,warning=F}
stc2.dt=readRDS(file = '~/HEARTSVG/data/HEARTSVG_demo.RDS')
dim(stc2.dt)
## [1] 4174 15429
```
If your data is of type seuratObject, we provide two functions, **filter.gene** and **seurat.trans** , for converting the data into a format suitable for HEARTSVG.
We use the function **filter.gene** to filter the non-expressed genes and low expression proportion genes. Then we use the function **seurat.trans** to transform the data type, from SeuratObject to a data.frame.

```{r}
## stc2=filter.gene(stc2,thr=0.01) # thr=0.01 represents the minimum expression proportion of each gene across all cells.
## dim(stc2)
## [1] 8812 4174

## stc2.dt=seurat.trans(stc2)
## dim(stc2.dt)
## [1] 4174 8814
```


## 1 Find SVG

The input data of the function **heartsvg** is a data.frame with dimensions n*(p+2) contains gene coordinates and gene expression counts. 

The first 2 columns should represent coordinates and their column names must be 'row' and 'col', respectively. The remaining p columns correspond to the expression levels of p genes. There are n rows, representing n spots.
```{r}
demo=stc2.dt[,1:1002]
demo[1:6,1:6]
##                    row col SAMD11 NOC2L HES4 ISG15
## AAACAACGAATAGTTC-1   0  16      0     0    1     1
## AAACAAGTATCTCCCA-1  50 102      0     2    0     2
## AAACAATCTACTAGCA-1   3  43      0     1    0     0
## AAACACCAATAACTGC-1  59  19      1     0   11     6
## AAACAGAGCGACTCCT-1  14  94      0     0    1     2
## AAACAGCTTTCAGAAG-1  43   9      0     0    0     2
```


```{r}
## find SVG
result=heartsvg(demo)
head(result)
##       gene pval p_adj rank
## 729  RPS27    0     0    1
## 604   PIGR    0     0    2
## 746 S100A6    0     0    3
## 185  CSRP1    0     0    4
## 725  RPL11    0     0    5
## 732   RPS8    0     0    6
```

## 2 Visualization
Different visualizations of top SVGs.

### 2.1 2D plot
```{r}
svg.2Dplot(data = demo,gene=result$gene[1:3],method = 'raw')
```

![image](https://user-images.githubusercontent.com/57090974/227206644-dfdc02ab-94c6-480d-81a3-f92cb0dcd8a3.png)
![image](https://user-images.githubusercontent.com/57090974/227206706-d284478f-4c75-43d8-9855-c99b65c694be.png)
![image](https://user-images.githubusercontent.com/57090974/227206753-c0825a0e-8c85-4beb-8c73-fc0971ca9469.png)



### 2.2 Marginal expression plot of SVG

The marginal expression of genes are obtained by the semi-pooling process.

```{r}
svg.Mplot(data=demo,gene=result$gene[1:3],method = 'raw')
```

![image](https://user-images.githubusercontent.com/57090974/227206902-da9df34a-803b-4e73-a503-b3e10cc5d23f.png)
![image](https://user-images.githubusercontent.com/57090974/227206963-d036d030-827c-4791-8df3-79aa6a7a9e19.png)
![image](https://user-images.githubusercontent.com/57090974/227207038-98177ad3-46d6-48c8-8f82-9bfb5145f8e4.png)



### 2.3 3D plot

#### 2.3.1 3D scatter plot

```{r}
svg.3Dplot(data=demo,gene = c('RPS8'),type = 'scatter')

```

![newplot](https://user-images.githubusercontent.com/57090974/227207527-221c9d5c-8a7b-47da-8bd7-1f2686dd60b0.png)


#### 2.3.1 3D surface plot

```{r}
svg.3Dplot(data=demo,gene=c('RPS27'),type = 'surface')

```

![image](https://user-images.githubusercontent.com/57090974/227207699-86fa6205-3293-40fe-bbc0-3c74d8b3345e.png)


## 3 Prediction of spatial domains

We offers an auto-clustering module for SVGs that facilitates further analyses, including predicting spatial domains, conducting functional analyses, and visualization based on the set of SVG files.  
We categorize SVGs based on their similarity in spatial expression and adopt a data-driven approach to estimate the number of spatially coherent regions in both SVG expression and histology.

### 3.1 clustering for SVGs
```{r}
clust_demo=svg.clust(data=demo,svg=result$gene[1:500],method = 'h')
table(clust_demo$cluster)
## 
##  1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17 
## 403  16  59   6   1   1   2   1   3   1   1   1   1   1   1   1   1 

```

```{r}
km_demo=svg.clust(data=demo,svg=result$gene[1:500],method = 'k',n.c=6)
table(km_demo$cluster)
## 
##   1   2   3   4   5   6 
## 114 139 123  27  68  29
```


### 3.2 visualization of spatial domains

```{r}
svg.2Dplot(data=demo,gene=clust_demogene[clustdemogene[clust_democluster==1],method = 'mean',title = 'domain-1')
```
![image](https://user-images.githubusercontent.com/57090974/227209204-cb369b7f-273b-4568-bc21-5477bb2ed4f0.png)


```{r}
svg.2Dplot(data=demo,gene=clust_demogene[clustdemogene[clust_democluster==2],method = 'mean',title = 'domain-2')
```

cl1.en=svg.enrichment(svg=clust_demo$gene[clust_demo$cluster==1])

head(cl1.en)


## 4 Enrichment

We perform enrichment analysis of SVG sets based on the package gprofiler2.

```{r}
cl1.en=svg.enrichment(svg=clust_demogene[clustdemogene[clust_democluster==1])

head(cl1.en)
##     query significant      p_value term_size query_size intersection_size
## 1 query_1        TRUE 5.007490e-08      5861        405               174
## 2 query_1        TRUE 2.599282e-07      4088        405               132
## 3 query_1        TRUE 5.064712e-07      6348        405               181
## 4 query_1        TRUE 5.795109e-07      5657        405               166
## 5 query_1        TRUE 1.015270e-06      6446        405               182
## 6 query_1        TRUE 2.018340e-06      2664        405                95
##   precision     recall    term_id source
## 1 0.4296296 0.02968777 GO:0048856  GO:BP
## 2 0.3259259 0.03228963 GO:0042221  GO:BP
## 3 0.4469136 0.02851292 GO:0048518  GO:BP
## 4 0.4098765 0.02934418 GO:0048522  GO:BP
## 5 0.4493827 0.02823456 GO:0032502  GO:BP
## 6 0.2345679 0.03566066 GO:0010033  GO:BP
##                                   term_name effective_domain_size source_order
## 1          anatomical structure development                 21128        13923
## 2                      response to chemical                 21128        10676
## 3 positive regulation of biological process                 21128        13622
## 4   positive regulation of cellular process                 21128        13626
## 5                     developmental process                 21128         8073
## 6             response to organic substance                 21128         3939
##                              parents
## 1                         GO:0032502
## 2                         GO:0050896
## 3             GO:0008150, GO:0050789
## 4 GO:0009987, GO:0048518, GO:0050794
## 5                         GO:0008150
## 6                         GO:0042221

```

![Rplot](https://user-images.githubusercontent.com/57090974/227210739-74e236d1-ab4b-4032-9b92-50565af266ad.png)


## Citation
Yuan X, Ma Y, Gao R, et al. HEARTSVG: a fast and accurate method for spatially variable gene identification in large-scale spatial transcriptomic data[J]. bioRxiv, 2023: 2023.08. 06.552154.
doi: https://doi.org/10.1101/2023.08.06.552154





