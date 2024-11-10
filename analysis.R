#---- Task CSF - analysis --------------------------------------------------####
# Date:         07.11.2024
# Author:       Katarzyna Chyzynska
# Filename:     analysis.R
# Description:  Differential expression analysis
#               > read in gene counts
#               > FOXP3 expression in stimulated vs non-treated cells
#               > most differentially expressed genes
# Project:      task_CSF
#----------------------------------------------------------------------------###

#---- Libraries ------------------------------------------------------------####

library(data.table)
library(magrittr)
library(edgeR)
library(ggplot2)

#----------------------------------------------------------------------------###

#---- Read and Process Data ------------------------------------------------####

foxp3_ensid <- "ENSG00000049768"

# sample info
samples_info  <- data.table(
  sample_name = c("SRR5223500", "SRR5223522", "SRR5223543",
                  "SRR5223505", "SRR5223547", "SRR5223570"),
  condition = c(rep("control", 3), rep("treated", 3))
) %>% .[order(sample_name)]

# RNA-seq count data
count_data <- read.table("./alignments/counts.txt", header = TRUE) %>%
  .[, c(1, 7:12)] %>%
  set_colnames(c("ensembl_id", samples_info$sample_name)) %>%
  DGEList()

# design for differential expression analysis
design <- model.matrix(~ samples_info$condition, data = count_data$samples)

# filter data
keep <- count_data %>% filterByExpr(design = design, group = design[,2])
count_data %<>% .[keep, , keep.lib.sizes=FALSE]

# normalize data (limma-voom)
norm_data <- count_data %>% voom(design, plot = TRUE)

# replicate correlations
norm_data$E[,!as.logical(design[,2])] %>% cor() # control
norm_data$E[, as.logical(design[,2])] %>% cor() # treated

#----------------------------------------------------------------------------###

#---- Differential Expression Analysis -------------------------------------####

# fit linear model and apply empirical Bayes smoothing
fit <- lmFit(norm_data, design)
fit <- eBayes(fit)

# get top differentially expressed genes
degs <- topTable(fit, coef = 2, number = Inf) %>% setDT()
top_degs <- subset(degs, (adj.P.Val < 0.05 & abs(logFC) > 1))

# visualize FOXP3
ggplot(degs, aes(logFC, -log10(adj.P.Val))) + geom_point() + theme_bw() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_vline(xintercept = c(-1,1), linetype = "dashed") +
  geom_point(data = degs[ensembl_id == foxp3_ensid,],
             aes(logFC, -log10(adj.P.Val)), colour="red", shape = 1, size = 3) +
  geom_text(data = degs[ensembl_id == foxp3_ensid,],
            label = "FOXP3", vjust = 2, color = "red")

#----------------------------------------------------------------------------###

