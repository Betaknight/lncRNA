###Must be run after scriptlncRNA.R
library(janitor)

trinotate <- read_delim(paste0(directorio, "/data/trinotate_annotation_report.xls"), 
                        delim = "\t", escape_double = FALSE, 
                        na = ".", trim_ws = TRUE)

trinotate <- clean_names(trinotate) %>% remove_empty()

#Separate the coordinates columns
trinotate <- trinotate %>% separate_wider_delim(prot_coords, delim = "[", names = c("prot_coordinates", "strand")) %>%
  mutate(strand = str_remove(strand, "]"))

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
            npeptides = sum(!is.na(prot_id)), nblasthits = sum(!is.na(sprot_top_blastx_hit))) %>%
  mutate(rna_type = case_when(
    npeptides > 0 ~ "Peptide",
    npeptides == 0 & nblasthits > 0 ~ "Nonsense",
    TRUE ~ "nonCoding"
            ))

