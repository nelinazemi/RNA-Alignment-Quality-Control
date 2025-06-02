library(tibble)
library(DESeq2)
library(pheatmap)
library(dplyr)
library(tidyr)
library(ggplot2)

raw_counts = read.delim("counts.txt", 
                      header = T, 
                      stringsAsFactors = F)

xp_design <- read.delim("experimental_design_modified.txt",
                        header = T,
                        stringsAsFactors = F)

xp_design$sample
colnames(raw_counts)


p_mean_sd_scaled <- 
  raw_counts %>% 
  as.data.frame() %>% 
  pivot_longer(cols = - gene, names_to = "sample", values_to = "counts") %>% 
  group_by(gene) %>%
  summarise(gene_average = mean(counts), gene_stdev = sd(counts)) %>% 
  ungroup() %>% 
  ggplot(., aes(x = log10(gene_average), y = log10(gene_stdev))) +
  geom_point(alpha = 0.5, fill = "grey", colour = "black") +
  labs(x = "Gene count average (log10 scale)",
       y = "Gene count standard deviation (log10 scale)") +
  ggtitle("Mean - Standard deviation relationship\n(no variance stabilisation ")

p_mean_sd_scaled

raw_counts <- raw_counts %>%
  column_to_rownames("gene")

head(raw_counts)

all(colnames(raw_counts) %in% xp_design$sample)

dds <- DESeqDataSetFromMatrix(countData = raw_counts, 
                              colData = xp_design, 
                              design = ~ seed + infected + dpi)

dds

dds = estimateSizeFactors(dds)

dds = estimateDispersions(object = dds, 
                          fitType = "parametric", 
                          quiet = TRUE)

sizeFactors(dds)

size_factors_df <- tibble(sample = names(dds$sizeFactor), correction_factor = dds$sizeFactor)

p <- ggplot(size_factors_df, aes(x = sample, y = correction_factor, group = 1)) +
  geom_point() + 
  geom_line() +
  theme(axis.text.x = element_text(angle = 90, size = 5)) +
  scale_y_continuous(limits = c(0.5, 2))

vsd = varianceStabilizingTransformation(object = dds, 
                                        blind = TRUE,
                                        fitType = "parametric")

variance_stabilised_counts <- assay(vsd)

head(variance_stabilised_counts)

p_mean_sd_vst <- 
  variance_stabilised_counts %>% 
  as.data.frame() %>% 
  rownames_to_column("gene") %>% 
  pivot_longer(cols = - gene, names_to = "sample", values_to = "counts") %>% 
  group_by(gene) %>% 
  summarise(gene_average = mean(counts), gene_stdev = sd(counts)) %>% 
  ungroup() %>% 
  ggplot(., aes(x = gene_average, y = gene_stdev)) +
  geom_point(alpha = 0.5, fill = "grey", colour = "black") +
  labs(x = "Gene count average (variance stabilised)", 
       y = "Gene count standard deviation (variance stabilised)") +
  ggtitle("Mean - Standard deviation relationship\n(after variance stabilisation ")

p_mean_sd_vst

# define a custom R function called "mypca()"
mypca <- function(x, center = TRUE, scale = TRUE){
  # remove constant variables
  constant_val = apply(x, 2, 'sd')
  x_reduced = x[, constant_val > 0]
  
  # perform SVD
  SVD <- svd(scale(x_reduced, center = center, scale = scale))
  
  # create scores data frame
  scores <- as.data.frame(SVD$u %*% diag(SVD$d))
  rownames(scores) <- rownames(x)
  colnames(scores) <- paste0("PC", c(1:dim(scores)[2]))
  
  # create loadings data frame
  loadings <- data.frame(SVD$v)
  rownames(loadings) <- colnames(x_reduced)
  colnames(loadings) <- paste0("PC", c(1:dim(loadings)[2]))
  
  # create data frame for explained variances
  explained_var <- as.data.frame(round((SVD$d^2) / sum(SVD$d^2)*100, digits = 1))
  rownames(explained_var) <- paste0("PC", c(1:dim(loadings)[2]))
  colnames(explained_var) <- "exp_var"
  
  # return result
  return (list("scores" = scores, "loadings" = loadings, "explained_var" = explained_var))
}

t_variance_stabilised_counts <- t(variance_stabilised_counts)

# before computing the PCA, check that samples are in rows and genes in columns
pca_results <- mypca(t_variance_stabilised_counts, 
                     center = TRUE, 
                     scale = TRUE)

ggplot(pca_results$explained_var, 
       aes(x = seq(from = 1, to = nrow(pca_results$explained_var)), 
           y = exp_var)) +
  ylab('explained variance (%)') + 
  ggtitle('Explained variance per component') + 
  geom_bar(stat = "identity") +
  labs(x = "Principal Component number") +
  scale_x_continuous(breaks = seq(
    from = 1,
    to = nrow(pca_results$explained_var)))

cumsum(pca_results$explained_var)

scores <- pca_results$scores

scores[1:5,1:5]

scores_with_conditions <- 
  scores %>% 
  rownames_to_column("sample") %>% # to prepare to join on the "sample" column
  left_join(x = .,                 # this means that we are passing the 'scores' dataframe 
            y = xp_design,         # this dataframe contains the sample to condition correspondence
            by = "sample")

head(scores_with_conditions)

explained_variance <- 
  pca_results$explained_var %>% 
  pull("exp_var")

ggplot(scores_with_conditions,
       aes(PC1, PC2, color = seed)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ",explained_variance[1],"% variance")) +
  ylab(paste0("PC2: ",explained_variance[2],"% variance")) + 
  coord_fixed(ratio = 1) +
  ggtitle("PCA score plot with the infection condition overlaid")