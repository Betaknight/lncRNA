Some_statistics
================
Jeronimo Miranda
2024-01-29

## Transcript classification

We took the information of which mRNAs are protein coding from
Transdecoder. More than 400,000 out of the 464,338 transcripts do not
code for a detectable peptide. Afterwards, to have an idea of
annotation, we use TRinotate that aligns with blastn and tblastn against
the uniprot database to find putative orthologs. There is also a small
but significant population of transcripts that have a significant match
at the RNA level but for which transdecoder could not find a long enough
ORF. We assume for now that these are either artifactual or
non-functional isoforms and are marked as *“nonsense”*. Later, if a gene
has at least one nonsense transcript, and no transcripts with coding
potential, we will classify this gene as a pseudogene. The following
graph info comes from the dataframe `transcript_classification`.

![](SomeStatistics_files/figure-gfm/Coding%20and%20non%20coding%20transcripts-1.png)<!-- -->

### Annotation guided more granular classification of transcripts

We have already classified non-coding transcripts using the Transdecoder
annotation as support: transcripts that could be annotated just as
non-coding are classified as non-sense, because they do show some
similarity to uniprot transcripts, even when they do not seem to code
for a peptide. The previous classification was done at the transcript
level. The next classification, pertains to peptides and so, will only
apply to the 60,247 transcripts for which transdecoder found at least
one peptide. The relationship between transcripts and peptides is
one-to-many, though, because one transcript has six possible reading
frames and could even code for more than one peptide in the same reading
frame. This is reflected in the 70,572 unique prot_id’s from the
transdecoder file. \#ANOTACION “…The next classification, pertains to
peptides and so, will only apply to the 60,247 transcripts…” Yo usaría
“pertaining” en lugarde “petains” pero checando con Google bard ambas
son correctas. Google bard sugiere esto para hacerlo mas formal: The
next classification, which pertains to peptides, will only apply to the
60,247 transcripts.

![](SomeStatistics_files/figure-gfm/coding%20transcripts-1.png)<!-- -->

There are 28 thousand peptides with a protein match to uniprot. This is
somewhat on the ballpark of the human genome? On the other hand, 31% of
the peptides are predicted to be on the opposite strand (-). This is a
lot, specially given our dataset strandedness. These transcripts could
correspond to antisense transcripts, though. Obviously, even if
transdecoder predicted ORFs on the opposite strand, if the assembly for
these transcripts is correct, the peptides cannot be produced. We
classify as *“Truncated”* the 3714 RNAs that code for a peptide, and
show similarity at the level of the RNA but not enough to show
similarity at the tblastn level. These transcripts could be added to the
possible pseudogene category, although of course this is not necesarilly
the case.

Perhaps the most interesting subset are the 15737 peptides that show no
similarity to any uniprot protein or mRNA. These *“Unnanotated”* dataset
could correspond to novel mollusk or gasteropod genes, however, we will
wait for the gene-level integration before making this call.

### Gene level integration

First, because of the one-to-many relationship between transcripts and
peptides, we made a further classification of transcripts according to
the peptides they produce. A transcript can code for many peptides, but
we classify a transcript as *Swissprot* annotated if at least one of the
peptides has a valid top ortholog from the SWissprot database. Next, a
transcript is annotated as either Truncated, Unnanotated or Wrong strand
if any of the peptides they produce belong to one of these categories
and none of the peptides they produce belong to the others. Transcripts
that produce a combination of these are not shown.

Next, we group transcripts by their Trinity assembly “gene” and classify
them as: - “Swissprot” genes if **any** of the transcripts produce a
peptide with a Swissprot ortholog. - We classified them as *Unnanotated*
if **all** of the transcripts code for an unnanotated protein (this
might be an underestimation). - We classify a transcript as *Non coding*
if **all** of their transcripts are non coding - We classify a
transcript as pseudogene if **all** of their transcripts belong to
either non coding, truncated or

![](SomeStatistics_files/figure-gfm/gene%20level-1.png)<!-- -->

## Length distribution of transcripts

We will now go into length distribution of transcripts. As stated
before, there are ~464,000 transcripts. By looking at a summary, we
already see that the distribution is heavily right skewed. The minimum
size is 169, and the maximum is 30kb (Spoiler: this is the Titin gene,
the longest known protein in humans). Still, the third quartile is 687.
That is, 75% of transcripts are shorter than the mean.

``` r
summary(longitudes$total_length)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   169.0   270.0   385.0   696.2   687.0 30666.0

To visualize this distribution, we log transform the data (exp = e)

![](SomeStatistics_files/figure-gfm/histogram%20whole%20transcriptome-1.png)<!-- -->

Remember that our first classification was, transcript-wise by their
coding potential, differentiating the peptide producing from the
non-coding, and a third smaller group of transcripts that have a blastx
hit but do not contain an ORF long enough to match by tblastn. There is
a very clear difference in the length distribution of transcripts
producing a peptide in transdecoder: almost all of them are longer than
the transcriptome median (log(389) = 5.95). There is also a sort of
bimodal distribution in the log transformed length data for the ORF
containing transcripts that could be explored later. Overall, we can
conclude that a good portion of the non-coding transcriptome is shorter
than 400 nucleotides. \#ANOTACION Creo que el inicio de este párrafo hay
que checarlo.

![](SomeStatistics_files/figure-gfm/log%20length%20by%20coding%20potential-1.png)<!-- -->

By just plotting the length of the ORF containing transcripts we see an
interesting pattern: the swissprot annotated distribution now loses the
bimodality and looks very nicely bell shaped, if a bit skewed. We also
see that most transcripts classified as truncated have a typical length
(similar to the swissprot matching ones). These transcripts are the ones
that have a similarity at the RNA level but the protein does not. This
raises the confidence that “truncated” transcripts correspond to
pseudogenes (“truncated” referring to the predicted ORF). Although the
“Unnanotated” ORF containg transcripts are exciting, the length
distribution shows that they tend to be shorter, and so care should also
be taken to ensure that they are truly unannotated and not just to short
to be matched to a known protein/transcript.

``` r
filter(trinotate, !is.na(prot_id)) %>% left_join(longitudes) %>% ggplot() + geom_histogram(aes(x = log(total_length))) + theme_bw() + labs(x = "log(transcript length) [nucleotides]", title = "Almost all ORF containing transcripts are longer than the median", subtitle = "Histogram of log transformed transcript length") + geom_vline(aes(xintercept = log(median(total_length))), color = "red") + facet_wrap(~peptide_type, nrow = 2, scales = "free_y")
```

![](SomeStatistics_files/figure-gfm/log%20length%20of%20peptide%20containing%20orfs,%20with%20annotations-1.png)<!-- -->

## Length distribution of ORFs

Now the obvious step is to check the distribution of ORFs. First, we
look at the distribution of ORF length, dependent on its transcript.
When looking only at **Swissprot** annotated ORFs, we see again the
lognormal distribution of transcript size (Sommer & Cohen, 1980). The
ORF size is not quite lognormal, though, with several peaks in the
log-transformed density plot.

![](SomeStatistics_files/figure-gfm/scatter%20side%20density%20mrna%20vs%20orf-1.png)<!-- -->

Now we look at the same scatterplot comparison for the Truncated
category. To remember, these were classified as **Truncated** for having
similarity at the transcript level but lacking it at the level of the
translated protein blast. This suggested that they are pseudogenes, that
are still similar to a transcript, but the predicted ORF is not big
enough to show enough similarity.

When looking at the scatterplot below, we have more confidence that this
is the case. These transcripts still show a roughly lognormal density
distribution, but most points are located well below the diagonal and
their density distribution on the y-axis shows that the vast majority of
the predicted ORFs are too short (around $e^6 =~ 400$ nucleotides).

![](SomeStatistics_files/figure-gfm/truncated%20ORFs%20size%20comparison-1.png)<!-- -->

The distribution of ORF length for those on the WRONG transcript strand
is similar to the “Truncated” category. Most transcripts are below
$e^{6.5}$ or around 665 nucleotides. This also supports the idea that
these transcripts correspond to assembly artifacts. However, their
position on the opposite strand could also indicate that they are
antisense transcripts, as we will see later.

![](SomeStatistics_files/figure-gfm/Wrong%20strand%20ORFs-1.png)<!-- -->

This is a test for the ORFs for the *Unannotated* transcripts, will they
look more like the *Swissprot* class or like the other two classes? The
shape of the ORFs density distribution is similar to “bad” transcripts
classes, however, taking into account that the unannotated transcripts
are also generally shorter, this is not concerning. The scatterplot does
look more like the one from the *Swissprot* genes because they are more
correlated with their transcript length (the area right at or below the
diagonal is more populated). Many of these transcripts likely code for
true novel gasteropod or mollusc transcripts.

![](SomeStatistics_files/figure-gfm/unnanotated%20ORF%20size%20distribution-1.png)<!-- -->

## Untranslated regions

Now we would like to look into the untranslated regions and check their
length distribution. This will be important when we map where the
antisense transcripts are hybridizing to the ORF.

#### 5’ Unstranslated regions

![](SomeStatistics_files/figure-gfm/utr5%20size-1.png)<!-- -->

The log10 transformed data for 5’ UTR sizes looks roughly log normal.
The very small UTRs showing as columns to the left of the x-axis could
be due to assembly/sequencing artifacts. Interestingly, both the
*unnanotated* and *Swissprot* categories are within 100-1000
nucleotides, while the 5’UTR of the *Truncated* ORF mRNAs are longer
than 1000.

#### 3’UTR

The sizes of 3’ UTRs are somewhat similar. If maybe they tend to be
longer, towards 1kb. There does not seem to be a huge difference in the
number 3’UTRs that are too short. In fact, I leave the warning for this
one because the log10 transformation erased 3’ UTRs with 0 values.
Something that did not happen with the 5’UTR.

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning: Removed 3608 rows containing non-finite values (`stat_bin()`).

![](SomeStatistics_files/figure-gfm/size%20of%203%20utr-1.png)<!-- -->

As said, the log10 transformed distributions look similar. This makes me
wonder whether there is a correlation between 5’ and 3’ utr lengths.
There does not seem to be a clear correlation, a transcript can have 0
length 5’ and 3’ or just one or the other. In general, though, there is
a tendency for *Swissprot* annotated transcripts to have a longer 3’UTR.
This does not hold true for *truncated* or *unnanotated* transcripts.

``` r
trinotate_wORFsize %>% filter(strand == "+") %>% transmute(transcript_id, utr5 = as.integer(prot_begin), prot_length, utr3 = as.integer(total_length) - as.integer(prot_end), total_length, peptide_type) %>%  ggplot(aes(utr3, utr5)) + geom_density2d_filled() + geom_point(alpha = 0.1, size = 0.1) + scale_x_log10() + scale_y_log10() + facet_wrap(~peptide_type) + theme_bw() + labs(title = "No correlation between 5' and 3' UTR sizes", subtitle = "log10 transformed UTR size for different peptide classes", x = "3' UTR length (nucleotides)", y ="5' UTR length (nucleotides)")  + theme(legend.position = 'none')
```

![](SomeStatistics_files/figure-gfm/utr%20correlation-1.png)<!-- -->

# Antisense transcripts

These sequences were then aligned back to the whole transcriptome with
blastn. The resulting alignments were filtered to be plus/minus
alignments and to span longer than 300 bp. The Target transcripts thus
found are labelled as antisense transcripts, because they match the
protein coding mRNAs at least partially on their antisense strand. Out
of the 39,844 unique transcripts with a swissprot annotation, we find
only 9930 with at least one antisense transcript. On the other hand, the
total number of unique antisense transcripts found this way is 26,953.

The length distribution of the two differ, with the antisense
transcripts being closer to the whole transcriptome numbers (i. e.
skewing towards smaller).

``` r
select(ampirr, total_length.qer, total_length.target) %>% transmute(Protein_coding_mRNAs = total_length.qer, Antisense_Transcripts = total_length.target) %>% summary()
```

    ##  Protein_coding_mRNAs Antisense_Transcripts
    ##  Min.   :  201        Min.   :  197.0      
    ##  1st Qu.:  453        1st Qu.:  291.0      
    ##  Median : 1102        Median :  413.0      
    ##  Mean   : 1828        Mean   :  722.2      
    ##  3rd Qu.: 2603        3rd Qu.:  731.0      
    ##  Max.   :20424        Max.   :16451.0

We can see that through all quartiles, Protein coding mRNAs are longer
than antisense Transcripts, with the median protein coding mRNA being
1kb while half of antisense transcripts are shorter than 413 nucleotides
(this is after already filtering all non coding RNAs shorter than 200
bp). Compare to the total length distribution of the transcriptome: 75%
of all transcripts have a length of less than 687 nucleotides:

    ##  transcript_id       total_length    
    ##  Length:464338      Min.   :  169.0  
    ##  Class :character   1st Qu.:  270.0  
    ##  Mode  :character   Median :  385.0  
    ##                     Mean   :  696.2  
    ##                     3rd Qu.:  687.0  
    ##                     Max.   :30666.0

``` r
#pivot to longer dataframe
ampirr %>% transmute(qseqid, sseqid, Protein_coding_mRNAs = total_length.qer, Antisense_Transcripts = total_length.target) %>% 
  pivot_longer(cols = Protein_coding_mRNAs:Antisense_Transcripts, names_to = "mRNA_type", values_to = "mRNA_length") %>% 
  #This operation is to chose the transcript id from the two seqid columns
  transmute(transcript_id = if_else(mRNA_type == "Protein_coding_mRNAs", qseqid, sseqid), mRNA_type, mRNA_length) %>% unique() %>% 
  #now you can plot
  ggplot(aes(x = log(mRNA_length))) + geom_histogram() + scale_x_continuous(n.breaks = 10) + theme_bw() + facet_wrap(~mRNA_type, nrow = 2, scales = "free_y") + labs(x = "log of mRNA length [nucleotides]", title = "Antisense transcripts are generally shorter than their targets, but more numerous", subtitle = "Histogram distribution of transcript length")
```

![](SomeStatistics_files/figure-gfm/Grafica%20de%20longitud-1.png)<!-- -->

## Bibliography

Sommer, S. S., & Cohen, J. E. (1980). The size distributions of
proteins, mRNA, and nuclear RNA. Journal of molecular evolution, 15(1),
37–57. <https://doi.org/10.1007/BF01732582>
