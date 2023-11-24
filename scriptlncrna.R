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
dlaevis_assembly_uniprt <- read_delim(paste0(directorio,"/data/dlaevis_assembly_uniprt.csv"), 
                                      delim = "\t", escape_double = FALSE, 
                                      col_types = cols(length = col_integer()), 
                                      trim_ws = TRUE)
                                      
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

ss_analysisAMP <- read_delim(paste0(directorio,"/data/ss_analysisAMP.dat"), 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)
ss_analysisIrr <- read_delim(paste0(directorio,"/data/ss_analysisIrr.dat"), 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)
View(ss_analysisIrr)

colnames(ss_analysisAMP)[1] <- "transcript_id"
colnames(ss_analysisIrr)[1] <- "transcript_id"

#Union de los reads plus, minus, total strands y diff ratio de AMP e IRR con query y target

ampirr <- left_join(swissProt_self_blastAStranscripts, ss_analysisAMP, by = join_by("qseqid" == "transcript_id"))
ampirr <- left_join(ampirr, ss_analysisIrr, by = join_by("qseqid" == "transcript_id"), suffix = c(".amp.qer", ".irr.qer"))
ampirr <- ampirr %>% mutate(diff_ratio_ampirr.qer = ((plus_strand_1stReads.amp.qer + plus_strand_1stReads.irr.qer) - (minus_strand_1stReads.amp.qer + minus_strand_1stReads.irr.qer))/(total_reads.amp.qer + total_reads.irr.qer))
ampirr <- left_join(ampirr, ss_analysisAMP, by = join_by("sseqid" == "transcript_id"))
ampirr <- left_join(ampirr, ss_analysisIrr, by = join_by("sseqid" == "transcript_id"), suffix = c(".amp.target", ".irr.target"))
ampirr <- ampirr %>% mutate(diff_ratio_ampirr.target = ((plus_strand_1stReads.amp.target + plus_strand_1stReads.irr.target) - (minus_strand_1stReads.amp.target + minus_strand_1stReads.irr.target))/(total_reads.amp.target + total_reads.irr.target))

#length of mRNAs
longitudes <- read_table(paste0(directorio,"/data/longitudes.txt"), 
                         col_names = FALSE)
View(longitudes)

colnames(longitudes) <- c("transcript_id", "total_length")
ampirr <- inner_join(ampirr, longitudes, by = join_by("qseqid" == "transcript_id"))
ampirr <- inner_join(ampirr, longitudes, by = join_by("sseqid" == "transcript_id"), suffix = c(".qer", ".target"))

#length distribution en escala log10 of mRNAs
select(ampirr, qseqid, total_length.qer) %>% unique() %>% ggplot() + geom_histogram(aes(x =total_length.qer)) + scale_x_log10()
# código para definir la escala del eje X y cúantos habrá- binwith= contenedor/ ancho bidireccional?* ggplot()) + geom_histogram(binwidth = 200) + scale_x_continuous(n.breaks = 10)

select(ampirr, qseqid, total_length.target) %>% unique() %>% ggplot() + geom_histogram(aes(x =total_length.target)) + scale_x_log10()

#olvidados son los datos entre 200 y 300 bases
olvidados <- filter(ampirr, length < 300)
diff_rati_olvidados <- left_join(select(olvidados, qseqid, sseqid),select(results_full_amp_degs,transcript_id, log2FoldChange, padj), by = join_by("qseqid" == "transcript_id")) %>% 
  left_join(select(results_full_amp_degs,transcript_id, log2FoldChange, padj), by =join_by("sseqid" == "transcript_id"), suffix = c(".qer.amp", ".sub.amp")) %>% 
  left_join(select(results_full_irr_degs, transcript_id, log2FoldChange, padj), by = join_by("qseqid" == "transcript_id")) %>% 
  left_join(select(results_full_irr_degs, transcript_id, log2FoldChange, padj), by = join_by("sseqid" == "transcript_id"), suffix = c(".qer.irr", ".sub.irr"))
#amp 318 datos
filter(diff_rati_olvidados, padj.qer.amp < 0.05 | padj.sub.amp < 0.05) #%>% left_join(select(dlaevis_assembly_uniprt, transcript_id,gene_names),by = join_by("qseqid" == "transcript_id")) %>%  View # hacer un left join con dlaevis para saber el nombre del gen
# Gráfica de puntos con valores de padj significativos.
filter(diff_rati_olvidados, padj.qer.amp < 0.05 | padj.sub.amp < 0.05) %>% ggplot(aes(x = log2FoldChange.sub.amp, y = log2FoldChange.qer.amp)) + geom_point()

#irr 205 datos
filter(diff_rati_olvidados, padj.qer.irr < 0.05 | padj.sub.irr < 0.05) #%>% View
# Gráfica de puntos con valores de padj significativos.
filter(diff_rati_olvidados, padj.qer.irr < 0.05 | padj.sub.irr < 0.05) %>% ggplot(aes(x = log2FoldChange.sub.amp, y = log2FoldChange.qer.amp)) + geom_point()


# Filtrado de matches mayores a 300

filtrados <- filter(ampirr, length > 300)

#Joining the Antisense blast table with log2foldchange and padj
qer_sub_diff <- left_join(select(filtrados, qseqid, sseqid),select(results_full_amp_degs,transcript_id, log2FoldChange, padj), by = join_by("qseqid" == "transcript_id")) %>% 
left_join(select(results_full_amp_degs,transcript_id, log2FoldChange, padj), by =join_by("sseqid" == "transcript_id"), suffix = c(".qer.amp", ".sub.amp")) %>% 
left_join(select(results_full_irr_degs, transcript_id, log2FoldChange, padj), by = join_by("qseqid" == "transcript_id")) %>% 
left_join(select(results_full_irr_degs, transcript_id, log2FoldChange, padj), by = join_by("sseqid" == "transcript_id"), suffix = c(".qer.irr", ".sub.irr"))

#query = RNA protein coding
#target antisense transcript

#amp
filter(qer_sub_diff, padj.qer.amp < 0.05 | padj.sub.amp < 0.05) # %>% left_join(select(dlaevis_assembly_uniprt, transcript_id,gene_names),by = join_by("qseqid" == "transcript_id")) %>%  View # hacer un left join con dlaevis para saber el nombre del gen
#Gráfica de puntos con valores de padj significativos.
filter(qer_sub_diff, padj.qer.irr < 0.05 | padj.sub.irr < 0.05) %>% ggplot(aes(x = log2FoldChange.sub.amp, y = log2FoldChange.qer.amp)) + geom_point()
#irr
filter(qer_sub_diff, padj.qer.irr < 0.05 | padj.sub.irr < 0.05) %>% View
#Gráfica de puntos con valores de padj significativos.
filter(qer_sub_diff, padj.qer.amp < 0.05 | padj.sub.amp < 0.05) %>% ggplot(aes(x = log2FoldChange.sub.amp, y = log2FoldChange.qer.amp)) + geom_point()

###
#scripts para poder 
###

#Creacion de un DF con el inicio y fin de la prot
prot <- select(dlaevis_assembly_uniprt, transcript_id, prot_coords) %>% separate(prot_coords, c("prot_start", "prot_end"), sep = "-") %>% mutate(prot_end = str_remove(prot_end, "[+]"))
prot$prot_end <- gsub("[\\[\\]", "", prot$prot_end)
prot$prot_end <- gsub("]", "", prot$prot_end)
# convertir las strings a numeros
prot$prot_start <- as.numeric(prot$prot_start)
prot$prot_end <- as.numeric(prot$prot_end)
prot<- unique(prot)
prot<- filter(prot, prot_start != 0)
prot<- filter(prot, prot_end != 0)

# Strings antisentido generalmente tienen NA en BlastP
# juntar prot_coords (end y start), con las coordenadas de matches de query y en sentido positivo
prot <- left_join(select(filtrados, qseqid, sseqid, length, qstart, qend,sstart, send), prot, by = join_by("qseqid" == "transcript_id")) %>%  filter(prot_start < prot_end)

# el caso 3, incluye tambien a los del caso 5, y el 4 incluye tambien a los del caso 6, por eso el caso 5 y 6 van primero
prot$caso <- case_when(
  prot$qstart >= prot$prot_end & prot$qend > prot$prot_end ~ "caso 6",
  prot$qstart < prot$prot_start & prot$qend <= prot$prot_start ~ "caso 5",
  prot$qstart >= prot$prot_start & prot$qend <= prot$prot_end ~ "caso 1",
  prot$qstart <= prot$prot_start & prot$qend >= prot$prot_end & prot$length > (prot$prot_end - prot$prot_start) ~ "caso 2",
  prot$qstart <= prot$prot_start & prot$qend <= prot$prot_end ~ "caso 3",
  prot$qstart >= prot$prot_start & prot$qend >= prot$prot_end ~ "caso 4",
)

#descripción de casos
#caso 1 El match cae dentro del CDS
#caso 2 El match es tan largo que cae tanto en 5UTR, CDS y 3UTR
#caso 3 El match empieza en la zona 5UTR y termina en cds
#caso 4 El match empieza en la zona de CDS y termina en 3UTR
#caso 5 El match es en la zona 5UTR
#caso 6 El match es en la zona 3UTR

prot <- mutate(prot, prot_end - prot_start) # + 1 para que coincida
colnames(prot)[11] <- "CDS_length"
prot <- left_join(prot, longitudes, by = join_by("qseqid" == "transcript_id"))
prot <- mutate(prot, total_length - prot_end)
colnames(prot)[13] <- "UTR_3"

# mutate(CDS_length = prot_end - prot_start)
#largo del 5UTR = inicio -1
#select(prot, qseqid, total_length) %>% unique() %>% ggplot(aes(x = total_length)) + geom_histogram() + labs(x = "bp", title = "Total_length") + theme(plot.title = element_text(hjust = 0.5))

#script para practicar las tablas pivotales
#prot %>% group_by(qseqid) %>% summarise(no = n()) %>%  left_join(select(prot,qseqid, total_length), by = "qseqid") %>%  View
#prot %>% group_by(qseqid) %>% summarise(No = n(), suma = sum(length), total = unique(total_length))

prot %>% group_by(caso) %>% summarise(No = n())

#Cambiar el tipo de dato de double a integer
#columnas <- c("length",qstart","qend","prot_start","prot_end","CDS_length","total_length","UTR_3")
#prot[columnas] <- lapply(prot[columnas], as.integer)

#Lógica
# Porcentaje de 5UTR,CDS y 3UTR cubierto con el transcrito antisentido
#caso 1 Dado que el alineamiento es en la zona CDS el porcentaje cubierto de CDS debe ser >=100% y las zonas de 5 y 3 UTR deben tener 0%
#caso 2 El alineamiento se da en las zonas 5utr y cds por lo tanto el porcentaje de la zona 3utr debe ser 0
#caso 3 prot end - qstart = x% en CDS, y% en 3UTR total lenght- prot end = 100%, prot end - qend = y%, 
#caso 4 
#caso 5 100% del match en zona 5utr, 0 en las demás
#caso 6 100% del match en la zona 3utr, 0 eb las demás
#Script

prot$UTR5porcentaje <- case_when(
  prot$caso == "caso 1" ~ 0,
  prot$caso == "caso 2" ~ ((prot$prot_start-prot$qstart)*100/(prot$prot_start-1)),
  prot$caso == "caso 3" ~ ((prot$prot_start-prot$qstart)*100/(prot$prot_start-1)),
  prot$caso == "caso 4" ~ 0,
  prot$caso == "caso 5" ~ 100,
  prot$caso == "caso 6" ~ 0,
)

prot$CDSporcentaje <- case_when(
  prot$caso == "caso 1" ~ (prot$length*100/prot$CDS_length),
  prot$caso == "caso 2" ~ 100,
  prot$caso == "caso 3" ~ ((prot$qend-prot$prot_start)*100/prot$CDS_length),
  prot$caso == "caso 4" ~ ((prot$prot_end-prot$qstart)*100/prot$CDS_length),
  prot$caso == "caso 5" ~ 0,
  prot$caso == "caso 6" ~ 0,
)

prot$UTR3porcentaje <- case_when(
  prot$caso == "caso 1" ~ 0,
  prot$caso == "caso 2" ~ ((prot$qend-prot$prot_end)*100/prot$UTR_3),
  prot$caso == "caso 3" ~ 0,
  prot$caso == "caso 4" ~ ((prot$qend-prot$prot_end)*100/prot$UTR_3),
  prot$caso == "caso 5" ~ 0,
  prot$caso == "caso 6" ~ 100,
)

left_join(select(prot,qseqid,sseqid,caso), qer_sub_diff, by = "qseqid") %>% filter(caso == "caso 6") %>% 
ggplot(aes(x = log2FoldChange.sub.amp, y = log2FoldChange.qer.amp)) + geom_point() + labs(title = "DEG relationship between Query and Subject" , subtitle = "DEG Amp", x = "log2foldchange antisense transcript" , y = "log2foldchange RNA protein coding") + theme(plot.title = element_text(hjust = 0.5))


#cor(x$log2FoldChange.sub.amp, x$log2FoldChange.qer.amp, use = "pairwise.complete.obs")
#cor.test(x$log2FoldChange.sub.amp, x$log2FoldChange.qer.amp,alternative = "greater", use = "pairwise.complete.obs")

prot %>% ggplot(aes(x = caso)) + geom_histogram()

#agregar el dato de la mediana
#grafica de log, transcritos VS UTR 3 y 5
#ggplot(prot, aes(x = ?????)) + geom_histogram() + scale_x_log()

#Logica para saber que porcentaje del transcrito está en cada zona.Debe sumar 100%

prot$porcentaje_en5UTR <- case_when(
  prot$caso == "caso 1" ~ 0,
  prot$caso == "caso 2" ~ ((prot$prot_start-prot$qstart)*100/(prot$qend-prot$qstart)),
  prot$caso == "caso 3" ~ ((prot$prot_start-prot$qstart)*100/(prot$qend-prot$qstart)),
  prot$caso == "caso 4" ~ 0,
  prot$caso == "caso 5" ~ 100,
  prot$caso == "caso 6" ~ 0,
)

prot$porcentaje_enCDS <- case_when(
  prot$caso == "caso 1" ~ 100,
  prot$caso == "caso 2" ~ 100,
  prot$caso == "caso 3" ~ ((prot$qend-prot$prot_start)*100/(prot$qend-prot$qstart)),
  prot$caso == "caso 4" ~ ((prot$prot_end-prot$qstart)*100/(prot$qend-prot$qstart)),
  prot$caso == "caso 5" ~ 0,
  prot$caso == "caso 6" ~ 0,
)

prot$porcentaje_en3UTR <- case_when(
  prot$caso == "caso 1" ~ 0,
  prot$caso == "caso 2" ~ ((prot$qend-prot$prot_end)*100/(prot$qend-prot$qstart)), #, scientific = FALSE, digits = 2)), 
  prot$caso == "caso 3" ~ 0,
  prot$caso == "caso 4" ~ ((prot$qend-prot$prot_end)*100/(prot$qend-prot$qstart)), #, scientific = FALSE, digits = 2)),
  prot$caso == "caso 5" ~ 0,
  prot$caso == "caso 6" ~ 100,
)
#sugiero menor a 4.8% porque en general estas son 40 bases o menos. hay un transcrito largo (5mil) que tiene un 4.8% en 3UTR que representa aprox 250 bases
#ilter(prot, porcentaje_en5UTR < 5 & porcentaje_en5UTR > 0) %>%  View
#o que por caso, la diferencia sea de 10 bases o menos