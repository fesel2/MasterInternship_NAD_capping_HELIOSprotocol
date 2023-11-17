library(biomaRt)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
#args = c("tmp/query_list.csv", "tmp/mart_query_afterR.csv")

library(biomaRt)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least two arguments must be supplied (input file, output file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  stop("At least two arguments must be supplied (input file, output file).n", call.=FALSE)
}

query_list <- read.csv(args[1], sep=",")

ensemble = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
mart_query <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'description',
                                 'gene_biotype', 'entrezgene_id', 'start_position', 
                                 'end_position', 'chromosome_name'),
                    filters = 'ensembl_gene_id', values = query_list$Ensemble_ID, mart = ensemble)
write.csv(mart_query, args[2], row.names = FALSE)
