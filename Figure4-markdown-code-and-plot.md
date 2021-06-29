Figure4 markdown code and plot
================
Wanglab
2021.2.19

``` r
library(ggplot2)
```

    ## Warning: package 'ggplot2' was built under R version 4.0.5

``` r
library(dplyr)
```

    ## Warning: package 'dplyr' was built under R version 4.0.5

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
set.seed(1)
range1 <- c(19,19,24,24,24,24,24,24,24,24,24,24,19,19)
range2 <- c(28,22,35,35,35,35,43,47,35,35,35,35,26,28)

path_DirEightSpieces = 'ContainedData/Meta2dResult/AllSpecies/'
  vector_FilenameAllTissues = c("Human_blood.csv","Human_skin.csv","Mouse_LIV.csv","Mouse_KID.csv","Fly_old.csv","Fly_young.csv","Yeast_high.csv","Yeast_low.csv","Neurospora.csv","Arabidopsis_SD LEAF.csv","Arabidopsis_SD M.csv",
                                "Arabidopsis_SD VA.csv","Chlamydomounas.csv","Cyanobacteria.csv")
  vector_NameAllTissues = c("Human_blood","Human_skin","Mouse_LIV","Mouse_KID","Fly_old","Fly_young","Yeast_high","Yeast_low","Neurospora","Arabidopsis_SD LEAF","Arabidopsis_SD M",
                            "Arabidopsis_SD VA","Chlamydomounas","Cyanobacteria")
  id_ShowNumAllTissues = c(1:14)
{
  par(mfrow = c(4,4))
  table2 = c()
  set.seed(1)
  for(i in id_ShowNumAllTissues)
  {
    {
      expr = read.csv(paste0(path_DirEightSpieces,vector_FilenameAllTissues[i]))
      ind_pc = expr$CycID != "#N/A"
      expr = expr[ind_pc,]   
      expr_mean = log(expr$meta2d_Base)
      expr_amp = log(expr$meta2d_AMP)
      expr_range <- apply(expr[,c(range1[i]:range2[i])],1, function(x){return(max(x) - min(x))})
      expr_ramp = log(expr$meta2d_rAMP)
      ind_bhq005 = expr$meta2d_BH.Q<0.05
    }
    ind_cir = 'ind_bhq005'
    aa0 = cbind(expr_mean,expr_range,ind_bhq005)
    aa = aa0[order(aa0[,'expr_range'],decreasing = T),];thresh_ind =
      which(aa[,ind_cir]==1)[round(sum(aa[,ind_cir])*0.95)];aa = aa[1:thresh_ind,]
    table1 = c(nrow(aa0),nrow(aa),round(nrow(aa)/nrow(aa0),2),
               sum(aa0[,ind_cir]),sum(aa[,ind_cir]),round(sum(aa[,ind_cir])/sum(aa0[,ind_cir]),2))
    aa = aa[order(aa[,'expr_mean']),]
    dd = exp(aa[,'expr_mean'])
    ee = replicate(1000,{
      # print(a <<- a+1)
      ind_simu = sample(1:nrow(aa),size = sum(aa[,ind_cir]),prob = dd)
      sum(aa[ind_simu,ind_cir])
    })
    ee1 = replicate(1000,{
      ind_simu = sample(1:nrow(aa0),size = sum(aa0[,ind_cir]))
      sum(aa0[ind_simu,ind_cir])
    })
    if (i==19) {
      ee1 = ee1+abs(rnorm(length(ee1),0,0.5))
    }
    m = sum(aa0[,ind_cir])
    n = nrow(aa0) - sum(aa0[,ind_cir])
    k = sum(aa0[,ind_cir])
    dhydata = dhyper(0:m,m,n,k) 
    plot(
      density(ee)
      ,xlim = c(min(ee1),max(ee))
      ,ylim = c(0,max(table(ee)/length(ee),dhydata[dhydata>0.001]))
      ,xlab = NA,main = vector_NameAllTissues[i],ylab = 'frequency'
      ,type = 'l',xaxt = 'n')
    abline(v=mean(ee),lty = 2)
    axis(side = 1)
    lines(density(ee1),col='red')
    abline(v=mean(ee1),lty = 2,col='red')
    table1 = c(table1,c(round(mean(ee),2),round(m^2/(m+n),2),t.test(ee,mu = round(m^2/(m+n),2))$p.value,max(ee)))
    table2 = rbind(table2,table1)
    print(vector_NameAllTissues[i])
    print(table1)
  }
  rownames(table2) = vector_NameAllTissues[id_ShowNumAllTissues]
  colnames(table2) = c('num_all','num_sele','ratio','cir_all','cir_sele','ratio_cir','cir_real','cir_simu','t-test.p','cir_RealMax')
  print(table2)
  #write.csv(table2,'table2.csv')
  par(mfcol = c(1,1))
} 
```

    ## [1] "Human_blood"
    ##  [1] 16996.00 12648.00     0.74  2011.00  1910.00     0.95   440.58   237.95
    ##  [9]     0.00   478.00

    ## [1] "Human_skin"
    ##  [1] 16996.00 10849.00     0.64  1086.00  1032.00     0.95   112.75    69.39
    ##  [9]     0.00   138.00

    ## [1] "Mouse_LIV"
    ##  [1] 19788.00 13848.00     0.70  5170.00  4912.00     0.95  2382.97  1350.76
    ##  [9]     0.00  2454.00

    ## [1] "Mouse_KID"
    ##  [1] 19788.00 14029.00     0.71  5018.00  4767.00     0.95  2020.59  1272.50
    ##  [9]     0.00  2087.00

    ## [1] "Fly_old"
    ##  [1] 10212.00  8571.00     0.84  1756.00  1668.00     0.95   367.48   301.95
    ##  [9]     0.00   404.00

    ## [1] "Fly_young"
    ##  [1] 10212.00  8476.00     0.83  1817.00  1726.00     0.95   424.42   323.30
    ##  [9]     0.00   471.00

    ## [1] "Yeast_high"
    ##  [1] 6138.00 5406.00    0.88 4568.00 4340.00    0.95 3520.66 3399.58    0.00
    ## [10] 3554.00

    ## [1] "Yeast_low"
    ##  [1] 6138.00 5456.00    0.89 3867.00 3674.00    0.95 2526.54 2436.25    0.00
    ## [10] 2575.00

    ## [1] "Neurospora"
    ##  [1]  9.733000e+03  7.398000e+03  7.600000e-01  1.730000e+02  1.640000e+02
    ##  [6]  9.500000e-01  5.670000e+00  3.080000e+00 1.503982e-205  1.300000e+01

    ## [1] "Arabidopsis_SD LEAF"
    ##  [1] 19914.00 13929.00     0.70  5947.00  5650.00     0.95  2997.99  1775.98
    ##  [9]     0.00  3061.00

    ## [1] "Arabidopsis_SD M"
    ##  [1] 19914.00 16137.00     0.81  2569.00  2441.00     0.95   692.68   331.41
    ##  [9]     0.00   734.00

    ## [1] "Arabidopsis_SD VA"
    ##  [1] 19914.00 15206.00     0.76  2937.00  2790.00     0.95   837.53   433.16
    ##  [9]     0.00   886.00

    ## [1] "Chlamydomounas"
    ##  [1] 16104.00 14576.00     0.91 10303.00  9788.00     0.95  6744.61  6591.64
    ##  [9]     0.00  6814.00

![](Figure4-markdown-code-and-plot_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

    ## [1] "Cyanobacteria"
    ##  [1] 2490.00 1872.00    0.75  572.00  543.00    0.95  164.54  131.40    0.00
    ## [10]  190.00
    ##                     num_all num_sele ratio cir_all cir_sele ratio_cir cir_real
    ## Human_blood           16996    12648  0.74    2011     1910      0.95   440.58
    ## Human_skin            16996    10849  0.64    1086     1032      0.95   112.75
    ## Mouse_LIV             19788    13848  0.70    5170     4912      0.95  2382.97
    ## Mouse_KID             19788    14029  0.71    5018     4767      0.95  2020.59
    ## Fly_old               10212     8571  0.84    1756     1668      0.95   367.48
    ## Fly_young             10212     8476  0.83    1817     1726      0.95   424.42
    ## Yeast_high             6138     5406  0.88    4568     4340      0.95  3520.66
    ## Yeast_low              6138     5456  0.89    3867     3674      0.95  2526.54
    ## Neurospora             9733     7398  0.76     173      164      0.95     5.67
    ## Arabidopsis_SD LEAF   19914    13929  0.70    5947     5650      0.95  2997.99
    ## Arabidopsis_SD M      19914    16137  0.81    2569     2441      0.95   692.68
    ## Arabidopsis_SD VA     19914    15206  0.76    2937     2790      0.95   837.53
    ## Chlamydomounas        16104    14576  0.91   10303     9788      0.95  6744.61
    ## Cyanobacteria          2490     1872  0.75     572      543      0.95   164.54
    ##                     cir_simu      t-test.p cir_RealMax
    ## Human_blood           237.95  0.000000e+00         478
    ## Human_skin             69.39  0.000000e+00         138
    ## Mouse_LIV            1350.76  0.000000e+00        2454
    ## Mouse_KID            1272.50  0.000000e+00        2087
    ## Fly_old               301.95  0.000000e+00         404
    ## Fly_young             323.30  0.000000e+00         471
    ## Yeast_high           3399.58  0.000000e+00        3554
    ## Yeast_low            2436.25  0.000000e+00        2575
    ## Neurospora              3.08 1.503982e-205          13
    ## Arabidopsis_SD LEAF  1775.98  0.000000e+00        3061
    ## Arabidopsis_SD M      331.41  0.000000e+00         734
    ## Arabidopsis_SD VA     433.16  0.000000e+00         886
    ## Chlamydomounas       6591.64  0.000000e+00        6814
    ## Cyanobacteria         131.40  0.000000e+00         190

``` r
library(plyr)
```

    ## ------------------------------------------------------------------------------

    ## You have loaded plyr after dplyr - this is likely to cause problems.
    ## If you need functions from both plyr and dplyr, please load plyr first, then dplyr:
    ## library(plyr); library(dplyr)

    ## ------------------------------------------------------------------------------

    ## 
    ## Attaching package: 'plyr'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     arrange, count, desc, failwith, id, mutate, rename, summarise,
    ##     summarize

``` r
library(ggpubr)
```

    ## 
    ## Attaching package: 'ggpubr'

    ## The following object is masked from 'package:plyr':
    ## 
    ##     mutate

``` r
set.seed(1)
path_DirEightSpieces = 'ContainedData/Meta2dResult/AllSpecies/'
  vector_FilenameAllTissues = c("Human_blood.csv","Human_skin.csv","Mouse_LIV.csv","Mouse_KID.csv","Fly_old.csv","Fly_young.csv","Yeast_high.csv","Yeast_low.csv","Neurospora.csv","Arabidopsis_SD LEAF.csv","Arabidopsis_SD M.csv",
                                "Arabidopsis_SD VA.csv","Chlamydomounas.csv","Cyanobacteria.csv")
vector_NameAllTissues = c("Human_blood","Human_skin","Mouse_LIV","Mouse_KID","Fly_old","Fly_young","Yeast_high","Yeast_low","Neurospora","Arabidopsis_SD LEAF","Arabidopsis_SD M",
                            "Arabidopsis_SD VA","Chlamydomounas","Cyanobacteria")
id_ShowNumAllTissues = c(1:14)
range1 <- c(19,19,24,24,24,24,24,24,24,24,24,24,19,19)
range2 <- c(28,22,35,35,35,35,43,47,35,35,35,35,26,28)
p = list()
table2 = c()
set.seed(1)
for(i in id_ShowNumAllTissues)
{
  {
    expr = read.csv(paste0(path_DirEightSpieces,vector_FilenameAllTissues[i]))
    ind_pc = expr$CycID != "#N/A"
    expr = expr[ind_pc,]    
    expr_mean = log(expr$meta2d_Base)
    expr_amp = log(expr$meta2d_AMP)
    expr_range <- log(apply(expr[,c(range1[i]:range2[i])],1, function(x){return(max(x) - min(x))}))
    #expr_ramp = log(expr$meta2d_rAMP)
    expr_ramp = log(expr$meta2d_rAMP)
    ind_bhq005 = expr$meta2d_BH.Q<0.05
  }
 
  ind_cir = 'ind_bhq005'
  aa0 = cbind(expr_mean,expr_range,ind_bhq005)

  aa = aa0[order(aa0[,'expr_range'],decreasing = T),];thresh_ind = 
    which(aa[,ind_cir]==1)[round(sum(aa[,ind_cir])*0.95)];aa = aa[1:thresh_ind,]
  
  table1 = c(nrow(aa0),nrow(aa),round(nrow(aa)/nrow(aa0),2),
             sum(aa0[,ind_cir]),sum(aa[,ind_cir]),round(sum(aa[,ind_cir])/sum(aa0[,ind_cir]),2))

  aa = aa[order(aa[,'expr_mean']),]
  dd = exp(aa[,'expr_mean'])
  ee = replicate(1000,{
    ind_simu = sample(1:nrow(aa),size = sum(aa[,ind_cir]),prob = dd)
    sum(aa[ind_simu,ind_cir])
  })
  ee1 = replicate(1000,{
    # print(a <<- a+1)
    ind_simu = sample(1:nrow(aa0),size = sum(aa0[,ind_cir]))
    sum(aa0[ind_simu,ind_cir])
  })
  # points(density(ee1),col='red',pch = '.')
  if (i==9) {
    ee1 = ee1+abs(rnorm(length(ee1),0,0.5))
    ee = ee+abs(rnorm(length(ee),0,0.5))
  }
  m = sum(aa0[,ind_cir])
  n = nrow(aa0) - sum(aa0[,ind_cir])
  k = sum(aa0[,ind_cir])
  dhydata = dhyper(0:m,m,n,k)
  df = data.frame(
    overlap = c(ee,ee1)
    ,method = rep(c('simulation','random'),c(length(ee),length(ee1))))
  mu <- ddply(df, "method", summarise, grp.mean=mean(overlap))
  p_i <- ggplot(df, aes(x=overlap, fill=method)) +
    geom_density(alpha=0.4) + 
    geom_vline(data=mu, aes(xintercept=grp.mean, color=method),
               linetype="dashed")+ 
    labs(title = vector_NameAllTissues[i])+
    theme(plot.title=element_text(hjust=0.5))
  p[[i]] = p_i
  
  table1 = c(table1,c(round(mean(ee),2),round(m^2/(m+n),2),t.test(ee,mu = round(m^2/(m+n),2))$p.value,max(ee)))
  table2 = rbind(table2,table1)
  print(vector_NameAllTissues[i])
  print(table1)
}
```

    ## [1] "Human_blood"
    ##  [1] 16996.00 12648.00     0.74  2011.00  1910.00     0.95   440.58   237.95
    ##  [9]     0.00   478.00
    ## [1] "Human_skin"
    ##  [1] 16996.00 10849.00     0.64  1086.00  1032.00     0.95   112.75    69.39
    ##  [9]     0.00   138.00
    ## [1] "Mouse_LIV"
    ##  [1] 19788.00 13851.00     0.70  5170.00  4912.00     0.95  2382.30  1350.76
    ##  [9]     0.00  2455.00
    ## [1] "Mouse_KID"
    ##  [1] 19788.00 14029.00     0.71  5018.00  4767.00     0.95  2020.59  1272.50
    ##  [9]     0.00  2087.00
    ## [1] "Fly_old"
    ##  [1] 10212.00  8571.00     0.84  1756.00  1668.00     0.95   367.48   301.95
    ##  [9]     0.00   404.00
    ## [1] "Fly_young"
    ##  [1] 10212.00  8476.00     0.83  1817.00  1726.00     0.95   424.42   323.30
    ##  [9]     0.00   471.00
    ## [1] "Yeast_high"
    ##  [1] 6138.00 5406.00    0.88 4568.00 4340.00    0.95 3520.66 3399.58    0.00
    ## [10] 3554.00
    ## [1] "Yeast_low"
    ##  [1] 6138.00 5456.00    0.89 3867.00 3674.00    0.95 2526.54 2436.25    0.00
    ## [10] 2575.00
    ## [1] "Neurospora"
    ##  [1]  9.733000e+03  7.398000e+03  7.600000e-01  1.730000e+02  1.640000e+02
    ##  [6]  9.500000e-01  6.060000e+00  3.080000e+00 2.309955e-240  1.361700e+01
    ## [1] "Arabidopsis_SD LEAF"
    ##  [1] 19914.00 13929.00     0.70  5947.00  5650.00     0.95  2997.18  1775.98
    ##  [9]     0.00  3056.00
    ## [1] "Arabidopsis_SD M"
    ##  [1] 19914.00 16137.00     0.81  2569.00  2441.00     0.95   693.56   331.41
    ##  [9]     0.00   743.00
    ## [1] "Arabidopsis_SD VA"
    ##  [1] 19914.00 15206.00     0.76  2937.00  2790.00     0.95   837.84   433.16
    ##  [9]     0.00   877.00
    ## [1] "Chlamydomounas"
    ##  [1] 16104.00 14576.00     0.91 10303.00  9788.00     0.95  6743.93  6591.64
    ##  [9]     0.00  6808.00
    ## [1] "Cyanobacteria"
    ##  [1] 2490.00 1872.00    0.75  572.00  543.00    0.95  164.62  131.40    0.00
    ## [10]  188.00

``` r
ggarrange(plotlist = p[id_ShowNumAllTissues]
          ,ncol = 4,nrow = 4
          # ,labels = vector_NameAllTissues[id_ShowNumAllTissues]
          ,labels = LETTERS[1:14]
          ,font.label = list(size = 12,face='bold',col='black')
          ,common.legend = T)
```

![](Figure4-markdown-code-and-plot_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->
