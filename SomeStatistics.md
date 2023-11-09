Some_statistics
================
Jeronimo Miranda
2023-11-09

## Transcript classification

We took the information of which mRNAs are protein coding from
Transdecoder. More than 400,000 out of the 464,338 transcripts do not
code for a detectable peptide. Afterwards, to have an idea of
annotation, we use TRinotate that aligns with blastn and tblastn against
the uniprot database to find putative orthologs. There is also a small
but significant population of transcripts that have a significant match
at the RNA level but do not have a significant tblastn alignment (mostly
because transdecoder could not find the corresponding ORF). We assume
for now that these are either artifactual or non-functional isoforms and
are marked as *“nonsense”*. Later, if a gene has at least one nonsense
transcript, and no transcripts with coding potential, we will classify
this gene as a pseudogene.

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
transdecoder file.

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
the peptides they produce.

![](SomeStatistics_files/figure-gfm/gene%20level-1.png)<!-- -->

## Length distribution of transcripts

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
  ggplot(aes(x = mRNA_length)) + geom_histogram(binwidth = 200) + scale_x_continuous(n.breaks = 10) + xlim(c(0,16500)) + theme_bw() + facet_wrap(~mRNA_type, nrow = 2, scales = "free_y") + labs(x = "mRNA length [nucleotides]", title = "Antisense transcripts are generally shorter than their targets, but more numerous", subtitle = "Histogram distribution of transcript length")
```

![](SomeStatistics_files/figure-gfm/Grafica%20de%20longitud-1.png)<!-- -->
