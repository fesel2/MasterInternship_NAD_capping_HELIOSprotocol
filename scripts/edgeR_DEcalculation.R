library(edgeR)

# each step is performed for all treatments

  # import data
norm_counts <- read.csv(
  "tmp/norm_counts.csv", row.names=1)

treatments = c("Control", "FK866", "NRH", "Rotenone")
names = c("Control", "FK", "NRH", "Rot")
i = 1

for (condition in treatments) {
  norm_counts_subset <- norm_counts[c(
    paste("bc03", condition, sep = "_"),
    paste("bc05", condition, sep = "_"),
    paste("bc07", condition, sep = "_"),
    paste("bc04", condition, sep = "_"),
    paste("bc06", condition, sep = "_"),
    paste("bc08", condition, sep = "_")
  )]
  
  # this is required to tell only look on genes which actually have reads
  smallest_row_sum = min(rowSums(norm_counts_subset))
  
  # create group binary
  group_binary <- factor(c(0, 0, 0, 1, 1, 1))
  
  # create a matrix from DataFrame
  hek_matrix <- data.matrix(norm_counts_subset)
  
  # Create a design matrix
  design_exact <- model.matrix(~ group_binary)
  
  # Create the DGEList object
  dgehek <- DGEList(counts=hek_matrix, group = group_binary)
  
  # Filter genes with low expression levels
  keep <- filterByExpr(min.total.count = smallest_row_sum + 1, dgehek)
  dgehek <- dgehek[keep, keep.lib.sizes = FALSE]
  
  # Calculate normalization factors and estimate dispersion
  #dgehek <- calcNormFactors(dgehek, method = "TMM")
  dgehek <- estimateDisp(dgehek)
  
  # Run the exact test
  exact_result <- exactTest(dgehek)
  
  # Extract top tags (genes)
  top_tags <- topTags(exact_result, n= 20)
  results_edgeR <- topTags(exact_result, n = nrow(dgehek), sort.by = "PValue")
  
  #write for export
  df <- results_edgeR$table
  write.csv(df, file = paste("tmp/", names[i], ".csv", sep = ""), row.names = TRUE)
  
  # clear variables for next iteration
  rm(dgehek)
  i = i+1
}
  
