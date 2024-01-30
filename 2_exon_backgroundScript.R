library(readr)
library(here)
library(dplyr)
library(ggplot2)
directorio <- here()


# load data ---------------------------------------------------------------
#exones intrones
gtf_trinity_exons <- read_delim(paste0(directorio,"/data/gtf_trinity_exons.csv"), 
                                delim = "\t", escape_double = FALSE, 
                                col_types = cols(exon_start = col_integer(), 
                                exon_end = col_integer()), trim_ws = TRUE)


Dinv_CNS_match_Dlaeve <- read_delim(paste0(directorio,"/data/Dinv_CNS_match_Dlaeve.txt"), 
                                    delim = "\t", escape_double = FALSE, 
                                    col_names = FALSE, col_types = cols(X3 = col_integer(), 
                                                                        X4 = col_integer(), X5 = col_integer(), 
                                                                        X6 = col_integer(), X7 = col_integer(), 
                                                                        X8 = col_integer(), X9 = col_integer()), 
                                    trim_ws = TRUE)

colnames(Dinv_CNS_match_Dlaeve) <- c("qseqid", "sseqid", "length", "qlen", "slen", "qstart", "qend", "sstart", "send", "evalue")

Dret_CNS_match_Dlaeve <- read_delim(paste0(directorio,"/data/Dret_CNS_match_Dlaeve.txt"), 
                                    delim = "\t", escape_double = FALSE, 
                                    col_names = FALSE, col_types = cols(X3 = col_integer(), 
                                                                        X4 = col_integer(), X5 = col_integer(), 
                                                                        X6 = col_integer(), X7 = col_integer(), 
                                                                        X8 = col_integer(), X9 = col_integer()), 
                                    trim_ws = TRUE)
colnames(Dret_CNS_match_Dlaeve) <- c("qseqid", "sseqid", "length", "qlen", "slen", "qstart", "qend", "sstart", "send", "evalue")

StrandedAssGCcontent <- read_table(paste0(directorio,"/data/StrandedAssGCcontent"), 
                                   col_types = cols(content = col_skip()))

#get the number of exons per transcript
gtf_trinity_exons <- gtf_trinity_exons %>% group_by(transcript_id) %>% mutate(number_of_exons = n(), exon_length = (1 + exon_end - exon_start), mean_exon_length = mean(exon_length)) %>%
  ungroup() %>% group_by(gene_id) %>% mutate(n_isoforms = n_distinct(transcript_id), gene_length = max(exon_end)) 

gtf_trinity_transcripts <- gtf_trinity_exons %>% ungroup() %>% group_by(transcript_id) %>% 
  summarise(gene_id = unique(gene_id),
            transcript_id = unique(transcript_id),
            n_isoforms = unique(n_isoforms),
            transcript_length = sum(exon_length),
            gene_length = unique(gene_length),
            mean_exon_length = unique(mean_exon_length),
            number_of_exons = unique(number_of_exons)
            ) %>% ungroup() %>% group_by(gene_id) %>% mutate(all_isoforms = sum(transcript_length), 
                                                             coverage = all_isoforms/(n_isoforms*gene_length))

gtf_trinity_transcripts <- inner_join(transcript_classification, gtf_trinity_transcripts) %>% 
  select(transcript_id, gene_id, rna_type, rna_type_2, rna_type_joint, total_length, transcript_length, gene_length, mean_exon_length, number_of_exons, coverage)

gtf_trinity_transcripts <- select(results_full_amp_degs, transcript_id, baseMean, log2FoldChange, padj) %>% inner_join(gtf_trinity_transcripts)

gtf_trinity_transcripts <- left_join(gtf_trinity_transcripts,ss_analysisAMP)

#############

data_donut <- gtf_trinity_transcripts %>% filter(baseMean != 0, abs(diff_ratio) > 0.75, mean_exon_length > 20, coverage > 0.4) %>% ungroup %>% group_by(rna_type) %>%
  summarise(rna_type = unique(rna_type), number = n()) %>%
  mutate(category = case_match(rna_type, "Nonsense" ~ "Nonsense", "nonCoding" ~ "Non coding", .default = "Peptide"), label = paste0(category, "\n n = ", number), ymax = cumsum(number), ymin = lag(ymax, default = 0))

gtf_trinity_transcripts <- left_join(gtf_trinity_transcripts,ss_analysisAMP)

ggtf_trinity_transcripts %>% ggplot(aes(n_exons)) + geom_histogram() + scale_x_log10()

gtf_trinity_transcripts %>% ggplot(aes(n_exons)) + geom_histogram(binwidth = 1)

inner_join(transcript_classification, gtf_trinity_transcripts)


