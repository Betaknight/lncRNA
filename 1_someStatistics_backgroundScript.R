###Must be run after scriptlncRNA.R
library(janitor)

trinotate <- read_delim(paste0(directorio, "/data/trinotate_annotation_report.xls"), 
                        delim = "\t", escape_double = FALSE, 
                        na = ".", trim_ws = TRUE)
make_clean_names(trinotate
                 )
