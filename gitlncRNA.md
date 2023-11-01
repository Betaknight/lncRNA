lncRNA
================
Bruno
2023-10-12

Librerias

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(tidyr)
library(forcats)
library(stringr)
library(ggplot2)
library(here)
```

    ## here() starts at C:/Users/mwaso/Documents/Proyecto_D_laeve/20-Transcriptomic_Bulk/24-lncRNA

``` r
directorio <- here()
source("scriptlncrna.R")
```

    ## New names:
    ## * `` -> `...1`
    ## * `` -> `...8`

    ## New names:
    ## Rows: 79530 Columns: 12
    ## -- Column specification
    ## -------------------------------------------------------- Delimiter: "\t" chr
    ## (3): X1, X2, X6 dbl (9): X3, X4, X5, X7, X8, X9, X10, X11, X12
    ## i Use `spec()` to retrieve the full column specification for this data. i
    ## Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## Rows: 224804 Columns: 5
    ## -- Column specification
    ## -------------------------------------------------------- Delimiter: "\t" chr
    ## (1): #transcript dbl (4): plus_strand_1stReads, minus_strand_1stReads,
    ## total_reads, diff_ratio
    ## i Use `spec()` to retrieve the full column specification for this data. i
    ## Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## Rows: 206854 Columns: 5
    ## -- Column specification
    ## -------------------------------------------------------- Delimiter: "\t" chr
    ## (1): #transcript dbl (4): plus_strand_1stReads, minus_strand_1stReads,
    ## total_reads, diff_ratio
    ## i Use `spec()` to retrieve the full column specification for this data. i
    ## Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## -- Column specification
    ## -------------------------------------------------------- cols( X1 =
    ## col_character(), X2 = col_double() )
    ## * `` -> `...1`
    ## * `` -> `...8`

You can check some preliminary analysis and summary statistics
[here](./SomeStatistics.md) Union de las columnas: plus, minus, total
strands y diff ratio de AMP e IRR con query y target Los histogramas son
respecto al dif-ratio de la suma de amp y irr para querry y target

``` r
ggplot(data = ampirr, aes(x = diff_ratio_ampirr.qer)) + geom_histogram() + labs(x = "diff ratio amp + irr", title = "Query") + theme(plot.title = element_text(hjust = 0.5)) + xlim(-1,1)
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

    ## Warning: Removed 2182 rows containing non-finite values (`stat_bin()`).

    ## Warning: Removed 2 rows containing missing values (`geom_bar()`).

![](gitlncRNA_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
ggplot(data = ampirr, aes(x = diff_ratio_ampirr.target)) + geom_histogram() + labs(x = "diff ratio amp + irr", title = "Target") + theme(plot.title = element_text(hjust = 0.5)) + xlim(-1,1)
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

    ## Warning: Removed 8510 rows containing non-finite values (`stat_bin()`).
    ## Removed 2 rows containing missing values (`geom_bar()`).

![](gitlncRNA_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

Expresion diferencial de

``` r
#amp
filter(qer_sub_diff, (log2FoldChange.qer.amp * log2FoldChange.sub.amp) < 0) %>% ggplot(aes(x = log2FoldChange.sub.amp, y = log2FoldChange.qer.amp)) + geom_point() + labs(title = "Amp") + theme(plot.title = element_text(hjust = 0.5))
```

![](gitlncRNA_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
#irr
filter(qer_sub_diff, (log2FoldChange.qer.irr * log2FoldChange.sub.irr) < 0) %>% ggplot(aes(x = log2FoldChange.sub.irr, y = log2FoldChange.qer.irr)) + geom_point() + labs(title = "Irr") + theme(plot.title = element_text(hjust = 0.5))
```

![](gitlncRNA_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->
