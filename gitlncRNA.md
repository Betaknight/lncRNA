lncRNA
================
Bruno
2023-10-12

Librerias

You can check some preliminary analysis and summary statistics
[here](./SomeStatistics.md)

Los histogramas son respecto al dif-ratio de la suma de amp y irr para
querry y target, la maypría de transcritos tienen

``` r
#Histogram of the ratio of reads on the positive and the negative strand (ideal is -1)
ggplot(data = ampirr, aes(x = diff_ratio_ampirr.qer)) + geom_histogram() + 
  labs(x = "Difference ratio of reads", title = "Query mRNAs show the expected predominance of minus strand reads", subtitle = "Difference ratio taking into account both experiments", y = "Number of transcripts") + theme(plot.title = element_text(hjust = 0.5)) + xlim(-1,1)
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

    ## Warning: Removed 2182 rows containing non-finite values (`stat_bin()`).

    ## Warning: Removed 2 rows containing missing values (`geom_bar()`).

![](gitlncRNA_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
ggplot(data = ampirr, aes(x = diff_ratio_ampirr.target)) + geom_histogram() + 
  labs(x = "Difference ratio of reads", title = "Antisense mRNAs show an even distribution of read strand", subtitle = "Difference ratio taking into account both experiments", y = "Number of transcripts") + theme(plot.title = element_text(hjust = 0.5)) + xlim(-1,1)
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

    ## Warning: Removed 8510 rows containing non-finite values (`stat_bin()`).
    ## Removed 2 rows containing missing values (`geom_bar()`).

![](gitlncRNA_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

Expresión diferencial para las condiciones de amputación y radiación

``` r
#Grafica de puntos de los cuadrantes +,- y -,+ #
#amp
filter(qer_sub_diff, (log2FoldChange.qer.amp * log2FoldChange.sub.amp) < 0) %>% ggplot(aes(x = log2FoldChange.sub.amp, y = log2FoldChange.qer.amp)) + geom_point() + labs(title = "DEG relationship between Query and Subject" , subtitle = "DEG Amp", x = "log2foldchange antisense transcript" , y = "log2foldchange RNA protein coding") + theme(plot.title = element_text(hjust = 0.5))
```

![](gitlncRNA_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
#irr
filter(qer_sub_diff, (log2FoldChange.qer.irr * log2FoldChange.sub.irr) < 0) %>% ggplot(aes(x = log2FoldChange.sub.irr, y = log2FoldChange.qer.irr)) + geom_point() + labs(title = "DEG relationship between Query and Subject" , subtitle = "DEG Irr", x = "log2foldchange antisense transcript" , y = "log2foldchange RNA protein coding") + theme(plot.title = element_text(hjust = 0.5))
```

![](gitlncRNA_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

Descripción de casos: caso 1 El alineamiento es solamente dentro del CDS
caso 2 El alineamiento es tan largo que cae en las 3 zonas (5UTR, CDS y
3UTR) caso 3 El alineamiento empieza en la zona 5UTR y termina dentro de
CDS caso 4 El alineamiento empieza en la zona de CDS y termina en 3UTR
caso 5 El alineamiento es solamente en la zona 5UTR caso 6 El
alineamiento es solamente en la zona 3UTR

``` r
prot %>% ggplot(aes(x = caso)) + geom_histogram(stat = "count") + labs(title = "La mayoría de alineamientos son solamente en la zona 3UTR" ,subtitle = "Distribución de casos")
```

![](gitlncRNA_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

\#Código para histogramas Histogramas Cuando el alineamiento es en la
zona de CDS, la mayoría de alineamientos en la zona de CDS se alinea con
el 50% del

``` r
prot %>% filter(CDSporcentaje>0) %>%   ggplot(aes(x = CDSporcentaje)) + geom_histogram()
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](gitlncRNA_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

Cuando el alineamiento es en la zona 3UTR, la gran mayoría de veces se
alinea con todo el largo del área 3UTR

``` r
prot %>% filter(UTR3porcentaje>0) %>%   ggplot(aes(x = UTR3porcentaje)) + geom_histogram()
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](gitlncRNA_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

Cuando el alineamiento es en la zona 5UTR, la mayoría de veces el
alineamiento es en todo el largo del área 5UTR

``` r
prot %>% filter(UTR5porcentaje>0) %>%   ggplot(aes(x = UTR5porcentaje)) + geom_histogram()
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](gitlncRNA_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

Correlación positiva para el caso 1 (el alineamiento es únicamente en la
zona CDS) Pearson’s product-moment correlation

data: x$log2FoldChange.sub.amp and x$log2FoldChange.qer.amp t = 22.313,
df = 21331, p-value \< 2.2e-16 alternative hypothesis: true correlation
is greater than 0 95 percent confidence interval: 0.1400006 1.0000000
sample estimates: cor 0.1510244

Explicación sobre los DF

Descripción del trabajo

DF de expresión diferencial: -results_full_amp_degs
-results_full_irr_degs

DF resultados del blast con transcritos antisentido:
-swissProt_self_blastAStranscripts

DF con anotaciones de BLASTX, prot ID, prot corrds, BLASTP, ontología,
nombre de los genes, organismos en el que se encuentra, length, etc
-dlaevis_assembly_uniprt

Data frames creadas:

-ampirr: 32 columnas, 79mil filas Tomando como base el DF swissProt, se
hizzo un left join de los DF de DEG. Además se hizo otro left join con
las longitudes

ampirr tiene 79mil entradas porque a veces un transcrito de query, hace
match con diferentes subjects. En total hay 9930 queries diferentes.

-Filtrados:32 collumnas, 23mil filas.

Tomando como base el DF ampirr, se hizo un filtro para seleccionar a los
transcritos mayores a 300. El df filtrado contiene 23639 entradas, es
decir, se quedarón fuera 50mil datos aprox

-qer_sub_diff:10 coulmnas, 23mil filas

tomando como base filtrados. Se hicieron left joins para agregar los
datos de log2FoldChange y padj para qseqid y sseqid de ambas condiciones
(amp y irr)

-olvidados: esta tabla tiene los 50mil datos que fueron excluidos en el
DF filtrados, es decir, aquí hay transcritos de entres 200 y 300
nucleotidos. -diff_rati_olvidados: al DF olvidados, se le agregó los
datos de log2FoldChange

-prot: Se seleccionó unicamente las columnas de trasncript_id y prot
coords del DF dlaevis_assembly_uniprt. Con código se separo prot coords
en dos columnas (prot_start y prot_end). Se eliminaron las filas con NA
en prot coord Se hizo un left join en prot con filtrados para agregar
los datos de “qseqid, sseqid, length, qstart, qend,sstart, send” Se
filtró para tener únicamente los datos en el que van en sentido Se
agregaron columnas con los datos de CDS_length y URT_3 Se hicieron las
columnas de los casos, porcentaje de zona cubierto, y porcentaje del lnc
en cada zona.
