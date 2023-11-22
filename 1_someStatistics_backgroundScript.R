###Must be run after scriptlncRNA.R
library(janitor)

# Load and clean ----------------------------------------------------------
trinotate <- read_delim(paste0(directorio, "/data/trinotate_annotation_report.xls"), 
                        delim = "\t", escape_double = FALSE, 
                        na = ".", trim_ws = TRUE)

trinotate <- clean_names(trinotate) %>% remove_empty()

#Separate the coordinates columns
trinotate <- trinotate %>% separate_wider_delim(prot_coords, delim = "[", names = c("prot_coordinates", "strand")) %>%
  mutate(strand = str_remove(strand, "]"))

# Classifications for donuts charts ----------------------------------------
#Classification of peptides
trinotate <- trinotate %>% mutate(peptide_type = case_when(
  strand == "-" ~ "Wrong strand",
  !is.na(sprot_top_blastp_hit) ~ "Swissprot",
  !is.na(sprot_top_blastx_hit) & !is.na(prot_id) & is.na(sprot_top_blastp_hit) ~ "Truncated",
  is.na(sprot_top_blastx_hit) & is.na(sprot_top_blastp_hit) & !is.na(prot_id) ~"Unnanotated"
)) 

#Because the transcripts are duplicated when more than one ortholog is found, we group by transcript

transcript_classification <- trinotate %>% group_by(transcript_id) %>% 
  summarise(number_gene_id = unique(number_gene_id),
            npeptides = sum(!is.na(prot_id)), nblasthits = sum(!is.na(sprot_top_blastx_hit)),
            at_least_one_swissprot = any(peptide_type == "Swissprot", na.rm = TRUE),
            at_least_one_wrongStrand = any(peptide_type == "Wrong strand", na.rm = TRUE),
            at_least_one_unnanotated = any(peptide_type == "Unnanotated", na.rm = TRUE),
            at_least_one_trunc       = any(peptide_type == "Truncated", na.rm = TRUE)
            ) %>%
  mutate(rna_type = case_when(
    npeptides > 0 ~ "Peptide",
    npeptides == 0 & nblasthits > 0 ~ "Nonsense",
    TRUE ~ "nonCoding"
            ),
    rna_type_2 = case_when(
      at_least_one_swissprot ~ "Swissprot",
      at_least_one_unnanotated & at_least_one_trunc ~ "2", #Doesnt exist
      at_least_one_wrongStrand & !at_least_one_unnanotated & at_least_one_trunc ~ "3",
      at_least_one_wrongStrand & at_least_one_unnanotated & !at_least_one_trunc ~ "Antisense and peptide",
      !at_least_one_wrongStrand & at_least_one_unnanotated & !at_least_one_trunc ~ "Unnanotated",
      at_least_one_wrongStrand & !at_least_one_unnanotated & !at_least_one_trunc ~ "Just wrong",
      !at_least_one_wrongStrand & !at_least_one_unnanotated & at_least_one_trunc ~ "True_ncated"
    ), rna_type_joint = coalesce(rna_type_2, rna_type)
    )


gene_classification <- transcript_classification %>% group_by(number_gene_id) %>% summarise(gene_type = case_when(
  any(rna_type_joint == "Swissprot") ~ "Swissprot",
  all(rna_type_joint == "Unnanotated") ~ "Unnanotated",
  all(rna_type_joint == "nonCoding") ~ "Non coding",
  all(rna_type_joint == "True_ncated" | rna_type_joint == "Nonsense" | rna_type_joint == "nonCoding") ~ "Pseudogene",
  TRUE ~ "OTHER"
)
)

# Length distribution of transcripts -------------------------------------
#Simple Joining dataframes to get transcript length

transcript_classification <- inner_join(transcript_classification, longitudes)

# Length distribution of ORFs

trinotate_wORFsize <- filter(trinotate, !is.na(prot_id)) %>% left_join(longitudes) %>% separate(prot_coordinates, sep = "-", into = c("prot_begin", "prot_end")) %>% mutate(prot_length = abs(as.integer(prot_begin) - as.integer(prot_end)))

## Differential expression of lncRNAs
diff_exp_amp <- filter(results_full_amp_degs, padj < 0.05)
diff_exp_amp <- left_join(diff_exp_amp, transcript_classification)
diff_exp_amp %>% group_by(rna_type_joint) %>% summarise(number = n())

diff_exp_irr <- filter(results_full_irr_degs, padj < 0.05)
diff_exp_irr <- left_join(diff_exp_irr, transcript_classification)
diff_exp_irr %>% group_by(rna_type_joint) %>% summarise(number = n())

#Which classification do antisense transcripts have?
antisense <- swissProt_self_blastAStranscripts %>% select(sseqid) %>% unique()
