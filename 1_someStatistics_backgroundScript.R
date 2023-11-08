###Must be run after scriptlncRNA.R
library(janitor)

trinotate <- read_delim(paste0(directorio, "/data/trinotate_annotation_report.xls"), 
                        delim = "\t", escape_double = FALSE, 
                        na = ".", trim_ws = TRUE)

trinotate <- clean_names(trinotate) %>% remove_empty()

#Because the transcripts are duplicated when more than one ortholog is found, we group by transcript

transcript_classification <- trinotate %>% group_by(transcript_id) %>% 
  summarise(number_gene_id = unique(number_gene_id),
            sprot_top_blastx_hit = unique(sprot_top_blastx_hit),
            npeptides = sum(!is.na(prot_id)), nblasthits = sum(!is.na(sprot_top_blastx_hit)),
            ntop_blastp_hit = sum(!is.na(sprot_top_blastp_hit))) %>%
  mutate(rna_type = case_when(
    npeptides > 0 ~ "Peptide",
    npeptides == 0 & nblasthits > 0 ~ "Nonsense",
    TRUE ~ "nonCoding"
            ))