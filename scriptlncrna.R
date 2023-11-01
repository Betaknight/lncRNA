library(readr)
#amputated_df_crbb_annotated_df_len_crbb_trin <- read_csv("C:/Users/cbzom/Downloads/amputated_df_crbb_annotated - df_len_crbb_trin.csv", na = c("NA",".",""))
#irradiated_df_crbb_annotation_df_len_crbb_trin <- read_csv("C:/Users/cbzom/Downloads/irradiated_df_crbb_annotation - df_len_crbb_trin.csv", na = c("NA","","."))
library(dplyr)
library(tidyr)
library(forcats)
library(stringr)
library(ggplot2)
library(here)
directorio <- here()
results_full_amp_degs <- read_csv(paste0(directorio,"/data/results_full_amp_degs.csv"), 
                                  col_types = cols(...8 = col_skip()))

results_full_irr_degs <- read_csv(paste0(directorio,"/data/results_full_irr_degs.csv"), 
                                  col_types = cols(...8 = col_skip()))

#cargar resutados de expresion diferencial
colnames(results_full_amp_degs)[1] <- "transcript_id"
colnames(results_full_irr_degs)[1] <- "transcript_id"

#Cargar blast results transcritos de antisentidos
swissProt_self_blastAStranscripts <- read_delim(paste0(directorio, "/data/swissProt_self_blastAStranscripts.csv"), 
                                                delim = "\t", escape_double = FALSE, 
                                                col_names = FALSE, trim_ws = TRUE)

colnames(swissProt_self_blastAStranscripts) <- c("qseqid", "sseqid", "length", "mismatch", "pident", "frames", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

ss_analysisAMP <- read_delim("ss_analysisAMP.dat", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)
ss_analysisIrr <- read_delim("ss_analysisIrr.dat", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)
View(ss_analysisIrr)







colnames(ss_analysisAMP)[1] <- "transcript_id"
colnames(ss_analysisIrr)[1] <- "transcript_id"

#Union de los plus, minus, total strands y diff ratio de AMP e IRR con query y target

ampirr <- left_join(swissProt_self_blastAStranscripts, ss_analysisAMP, by = join_by("qseqid" == "transcript_id"))
ampirr <- left_join(ampirr, ss_analysisIrr, by = join_by("qseqid" == "transcript_id"), suffix = c(".amp.qer", ".irr.qer"))
ampirr <- ampirr %>% mutate(diff_ratio_ampirr.qer = ((plus_strand_1stReads.amp.qer + plus_strand_1stReads.irr.qer) - (minus_strand_1stReads.amp.qer + minus_strand_1stReads.irr.qer))/(total_reads.amp.qer + total_reads.irr.qer))
ampirr <- left_join(ampirr, ss_analysisAMP, by = join_by("sseqid" == "transcript_id"))
ampirr <- left_join(ampirr, ss_analysisIrr, by = join_by("sseqid" == "transcript_id"), suffix = c(".amp.target", ".irr.target"))
ampirr <- ampirr %>% mutate(diff_ratio_ampirr.target = ((plus_strand_1stReads.amp.target + plus_strand_1stReads.irr.target) - (minus_strand_1stReads.amp.target + minus_strand_1stReads.irr.target))/(total_reads.amp.target + total_reads.irr.target))

ggplot(data = ampirr, aes(x = diff_ratio_ampirr.qer)) + geom_histogram() + labs(x = "diff ratio amp + irr", title = "Query") + theme(plot.title = element_text(hjust = 0.5)) + xlim(-1,1)
ggplot(data = ampirr, aes(x = diff_ratio_ampirr.target)) + geom_histogram() + labs(x = "diff ratio amp + irr", title = "Target") + theme(plot.title = element_text(hjust = 0.5)) + xlim(-1,1)


library(readr)
longitudes <- read_table("longitudes.txt", 
                         col_names = FALSE)
View(longitudes)

colnames(longitudes) <- c("transcript_id", "total_length")
ampirr <- inner_join(ampirr, longitudes, by = join_by("qseqid" == "transcript_id"))
ampirr <- inner_join(ampirr, longitudes, by = join_by("sseqid" == "transcript_id"), suffix = c(".qer", ".target"))

select(ampirr, qseqid, total_length.qer) %>% unique() %>% ggplot(aes(x = total_length.qer)) + geom_histogram(binwidth = 200) + scale_x_continuous(n.breaks = 10)

# Filtrado de matches mayores a 300

filtrados <- filter(ampirr, length > 300)

#experimento
qer_sub_diff <- left_join(select(filtrados, qseqid, sseqid),select(results_full_amp_degs,transcript_id, log2FoldChange, padj), by = join_by("qseqid" == "transcript_id")) %>% 
left_join(select(results_full_amp_degs,transcript_id, log2FoldChange, padj), by =join_by("sseqid" == "transcript_id"), suffix = c(".qer.amp", ".sub.amp")) %>% 
left_join(select(results_full_irr_degs, transcript_id, log2FoldChange, padj), by = join_by("qseqid" == "transcript_id")) %>% 
left_join(select(results_full_irr_degs, transcript_id, log2FoldChange, padj), by = join_by("sseqid" == "transcript_id"), suffix = c(".qer.irr", ".sub.irr"))

#obtener los cuadrantes +,- y -,+
#filter(qer_sub_diff, (log2FoldChange.qer.amp * log2FoldChange.sub.amp) < 0) %>% ggplot(aes(x = log2FoldChange.sub.amp, y = log2FoldChange.qer.amp)) + geom_point()
# irr : filter(qer_sub_diff, (log2FoldChange.qer.irr * log2FoldChange.sub.irr) < 0) %>% ggplot(aes(x = log2FoldChange.sub.irr, y = log2FoldChange.qer.irr)) + geom_point()


# Amp, padj menor a 0.05 de qer.
#filter(qer_sub_diff, (log2FoldChange.qer.amp * log2FoldChange.sub.amp) < 0) %>% filter(padj.sub.amp < 0.05) %>% ggplot(aes(x = log2FoldChange.sub.amp, y = log2FoldChange.qer.amp)) + geom_point()

# Irr 
# filter(qer_sub_diff, (log2FoldChange.qer.irr * log2FoldChange.sub.irr) < 0) %>% filter(padj.sub.irr < 0.05) %>% ggplot(aes(x = log2FoldChange.sub.irr, y = log2FoldChange.qer.irr)) + geom_point()

antisense_sig <- qer_sub_diff %>% filter(padj.sub.amp < 0.05)
# qer_sub_diff %>% filter(padj.sub.irr < 0.05) AUN NO ESTA ASIGNADO A UN OBJEO

# hacer le left join con elq uery, avance aun no asignado a un objeto
#left_join(antisense_sig, dlaevis_assembly_uniprt, by =join_by("qseqid" == "transcript_id")) %>% select(qseqid, sseqid, protein_names) %>% unique() %>% View
