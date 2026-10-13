title: "Nociceptive interoception drives a self-reinforcing neurovascular circuit in aortic dissection"


# Pseudotime --------------------------------------------------------------
## CytoTRACE2
library(CytoTRACE2)

AD_SMC <- cytotrace2(AD_SMC, is_seurat = TRUE, slot_type = "counts", species = "human", seed = 1234)
head(AD_SMC@meta.data)

annotation <- data.frame(phenotype=AD_SMC@meta.data$major_cell_type) %>% set_rownames(., colnames(AD_SMC))
plots <- plotData(cytotrace2_result=AD_SMC, annotation=annotation, is_seurat=TRUE)
p1 <- plots$CytoTRACE2_UMAP
p2 <- plots$CytoTRACE2_Potency_UMAP
p3 <- plots$CytoTRACE2_Relative_UMAP
p4 <- plots$CytoTRACE2_Boxplot_byPheno
(p1+p2+p3+p4) + plot_layout(ncol = 2)

update_geom_defaults("point", list(size = 0.1))
plots <- plotData(cytotrace2_result=AD_SMC, annotation=annotation, is_seurat=TRUE)

(plots$CytoTRACE2_UMAP +
    plots$CytoTRACE2_Potency_UMAP +
    plots$CytoTRACE2_Relative_UMAP +
    plots$CytoTRACE2_Boxplot_byPheno) +
  plot_layout(ncol = 2)

## Slingshot
library(slingshot)
AD_SMC_sce <- as.SingleCellExperiment(AD_SMC)
umap = AD_SMC@reductions$umap@cell.embeddings
colnames(umap) = c('UMAP-1', 'UMAP-2')
reducedDims(AD_SMC_sce) = SimpleList(UMAP = umap)
rd = umap
plot(rd,pch=16, asp = 1)
meta = AD_SMC@meta.data
colData(AD_SMC_sce)$orig.ident = meta$orig.ident
colData(AD_SMC_sce)$Group = meta$Group
colData(AD_SMC_sce)$Subgroup = meta$Subgroup
colData(AD_SMC_sce)$major_cell_type = meta$major_cell_type
sds <- slingshot(AD_SMC_sce,
                 clusterLabels = AD_SMC_sce@colData@listData[["major_cell_type"]],
                 reducedDim = 'UMAP',
                 reweight = FALSE)
SlingshotDataSet(sds)

library(RColorBrewer)
plot(reducedDims(sds)$UMAP, pch=16, asp = 1)
lines(SlingshotDataSet(sds), lwd=2, col=brewer.pal(9,"Set1"))

pseudotime_mat <- slingPseudotime(sds)
pseudotime_median <- apply(pseudotime_mat, 1, function(x) median(x, na.rm = TRUE))
AD_SMC[["slingshot"]] <- pseudotime_median
head(AD_SMC@meta.data)

mydata<- FetchData(AD_SMC,vars = c("UMAP_1","UMAP_2","slingshot"))
my_colors <- c("#334a55","#3a6f6f","#4d9a8e","#70a68b","#92af81","#a5b47e","#bcbc80","#e4c57b","#e6b471") #,"#e18d64","#e76e50"

ggplot(mydata, aes(x = UMAP_1, y = UMAP_2, colour = slingshot)) +
  geom_point(size = 0.1) +
  scale_color_gradientn(colors = my_colors) +
  theme_bw() +
  theme(
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black")
  ) + NoLegend()

## Monocle3
library(monocle3)

expression_matrix <- as(as.matrix(AD_SMC@assays$RNA@counts), 'sparseMatrix')
cell_metadata <- AD_SMC@meta.data

gene_metadata <- data.frame(gene_short_name = rownames(expression_matrix))
rownames(gene_metadata) <- rownames(expression_matrix)

cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_metadata)

cds <- preprocess_cds(cds, num_dim = 100)
cds <- reduce_dimension(cds)
umap_coords <- AD_SMC@reductions$umap@cell.embeddings
reducedDims(cds)$UMAP <- umap_coords

cds <- cluster_cells(cds)
cds <- learn_graph(cds)
cds <- order_cells(cds)
plot_cells(cds)
AD_SMC@meta.data$Monocle_pseudotime <- pseudotime(cds)
head(AD_SMC@meta.data)

mydata<- FetchData(AD_SMC,vars = c("UMAP_1","UMAP_2","Monocle_pseudotime"))
ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = Monocle_pseudotime))+
  geom_point(size = 0.1)+
  scale_color_viridis(option = "magma") +
  theme_bw() + theme(panel.border = element_blank(), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black")) + NoLegend()


# Pseudotime Quant --------------------------------------------------------
head(AD_SMC@meta.data)
meta <- AD_SMC@meta.data
meta <- meta[!is.na(meta$Monocle_pseudotime), ]
meta$major_cell_type <- as.factor(meta$major_cell_type)

df_prop <- meta %>%
  mutate(Monocle_pseudotime_bin = cut(Monocle_pseudotime, breaks = 10, include.lowest = TRUE)) %>%
  count(Monocle_pseudotime_bin, major_cell_type) %>%
  complete(Monocle_pseudotime_bin, major_cell_type, fill = list(n = 0)) %>%
  group_by(Monocle_pseudotime_bin) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# bin center
df_prop$bin_center <- sapply(as.character(df_prop$Monocle_pseudotime_bin), function(x) {
  x <- gsub("\\[|\\]|\\(|\\)", "", x)
  vals <- as.numeric(strsplit(x, ",")[[1]])
  mean(vals)
})

color3 = c("#5586B3", "#45958e", "#B6C796", "#9F77CF", "#ead4d4",  "#C98C9A", "#fad4a6","#72A77A", "#cee0e6", "#9ac1f0", "#f5c0bf", "#D66F62", "#C8A2C8")
ggplot(df_prop, aes(x = bin_center, y = prop, fill = major_cell_type)) +
  geom_area(position = "stack", alpha = 0.8) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_fill_manual(values = color3) +
  theme_classic(base_size = 14) +
  xlab("Monocle pseudotime") +
  ylab("Proportion") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    legend.title = element_blank()
  )

## get median pseudotime
table(AD_SMC@meta.data$major_cell_type)

avg_pseudotime <- AD_SMC@meta.data %>%
  group_by(major_cell_type) %>%
  summarise(mean_pseudotime = mean(Monocle_pseudotime, na.rm = TRUE)) %>%
  ungroup()

target_order <- c(
  "Synthetic SMC", "Mesenchymal-like SMC", "Immune-regulatory SMC", "IFN-enriched SMC",
  "Proinflammatory SMC", "Neuroresponsive SMC", "Contractile SMC", "Stressed SMC",
  "Secretory SMC", "Fibroblast-like SMC", "Macrophage-like SMC"
)

avg_pseudotime_ordered <- avg_pseudotime %>%
  mutate(major_cell_type = factor(major_cell_type, levels = target_order)) %>%
  arrange(major_cell_type)

print(avg_pseudotime_ordered)
avg_pseudotime_ordered$OR <- c(-0.55, 1.77, -0.59, -0.81, -0.7, 1.16, -0.47, -0.49, 1.62, -0.78, -0.15)
avg_pseudotime_filtered <- avg_pseudotime_ordered %>%
  mutate(major_cell_type = str_remove(major_cell_type, " SMC"))
print(avg_pseudotime_filtered)
avg_pseudotime_filtered <- as.data.frame(avg_pseudotime_filtered)
rownames(avg_pseudotime_filtered) <- avg_pseudotime_filtered[,1]
avg_pseudotime_filtered <- avg_pseudotime_filtered[,-1]

size <- c("8630","8090","6902","6359","4845","4530","3851","3384","2506","1988","1343")
avg_pseudotime_filtered$size.log <- log(as.numeric(size))


avg_pseudotime_filtered$pseudotime_scaled <- (avg_pseudotime_filtered$mean_pseudotime - min(avg_pseudotime_filtered$mean_pseudotime)) /
  (max(avg_pseudotime_filtered$mean_pseudotime) - min(avg_pseudotime_filtered$mean_pseudotime))

avg_pseudotime_filtered$OR_scaled <- (avg_pseudotime_filtered$OR - min(avg_pseudotime_filtered$OR)) / 
  (max(avg_pseudotime_filtered$OR) - min(avg_pseudotime_filtered$OR))


color3 = c("#5586B3", "#45958e", "#B6C796", "#9F77CF", "#ead4d4",  "#C98C9A", "#fad4a6","#72A77A", "#cee0e6", "#9ac1f0", "#f5c0bf")

ggplot(avg_pseudotime_filtered, aes(x = pseudotime_scaled, y = OR_scaled)) +
  geom_point(aes(size = size.log, color = color3)) +
  scale_color_identity() +
  scale_x_continuous(limits = c(0, 1)) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )


head(AD_SMC@meta.data)

library(dplyr)
library(tidyr)
library(ggplot2)


meta <- AD_SMC@meta.data


meta <- meta[!is.na(meta$Monocle_pseudotime), ]
meta$major_cell_type <- as.factor(meta$major_cell_type)
meta$Monocle_pseudotime <- as.numeric(meta$Monocle_pseudotime)


meta <- meta %>%
  mutate(Monocle_pseudotime_bin = cut(Monocle_pseudotime,
                                      breaks = 60,
                                      include.lowest = TRUE,
                                      labels = FALSE))


df_prop <- meta %>%
  group_by(Monocle_pseudotime_bin, major_cell_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Monocle_pseudotime_bin) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

ggplot(df_prop, aes(x = Monocle_pseudotime_bin, y = prop, fill = major_cell_type)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(breaks = seq(1, 30, by = 2)) +
  labs(x = "Monocle pseudotime bin", y = "Cell proportion") +
  scale_fill_manual(values = color3) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) #+ NoLegend()


# Suppl Fig
colnames(AD_SMC@meta.data)

library(dplyr)
library(ggplot2)

df <- AD_SMC@meta.data


df <- df %>%
  group_by(major_cell_type) %>%
  mutate(median_score = median(slingshot, na.rm = TRUE)) %>%
  ungroup()

df$major_cell_type <- factor(df$major_cell_type, 
                             levels = df %>% 
                               group_by(major_cell_type) %>% 
                               summarise(med = median(slingshot, na.rm = TRUE)) %>% 
                               arrange(med) %>% 
                               pull(major_cell_type))

ggplot(df, aes(x = major_cell_type, y = slingshot, fill = major_cell_type)) +
  geom_boxplot(outlier.shape = 21, outlier.color = "black", outlier.size = 2) +
  theme_classic() +
  xlab("Major Cell Type") +
  ylab("CytoTRACE2 value") +
  ggtitle("Monocle scores per major cell type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# SMC CytoSig --------------------------------------------------
install.packages("scaper")
library(scaper)
AD_SMC <- scapeForSeurat(AD_SMC, database = "cytosig", cytokine = "all", normalize = F)
head(AD_SMC@meta.data)
(as.data.frame(AD_SMC@assays$scape@data))[1:10,1:10]


score <- as.data.frame(AD_SMC@assays$scape@data)
score[1:4,1:4]
score_change <- score %>% 
  select_all(~str_replace_all(., "\\.", "-"))
identical(colnames(score_change) , rownames(AD_SMC@meta.data))
AD_SMC@meta.data <- cbind(AD_SMC@meta.data,t(score_change) )


df <- AD_SMC@meta.data
colnames(AD_SMC@meta.data)
avg_df = aggregate(df[,30:ncol(df)],list(df$major_cell_type),mean)

df <- AD_SMC@meta.data
colnames(AD_SMC@meta.data)
avg_df = aggregate(df[,30:ncol(df)],list(df$Group),mean)


cytokines <- c("ActivinA", "BDNF", "BMP2", "BMP4", "BMP6", "CD40LG", "CXCL12", 
               "EGF", "FGF2", "GCSF", "GDF11", "GMCSF", "HGF", "IFNG", "IFNL", 
               "IL10", "IL12", "IL13", "IL15", "IL17A", "IL1A", "IL1B", "IL2", 
               "IL21", "IL22", "IL27", "IL3", "IL4", "IL6", "LIF", "LTA", "MCSF", 
               "NO", "OSM", "TGFB1", "TGFB3", "TNFA", "TNFSF12", "TRAIL", 
               "VEGFA", "WNT3A")

plot_cytokine_vs_pseudotime <- function(cytokine_name) {
  ggplot(AD_SMC@meta.data, aes_string(x = "pseudotime", y = cytokine_name)) +
    geom_point() +
    labs(x = "Pseudotime", y = cytokine_name, title = paste(cytokine_name, "vs. Pseudotime")) +
    theme_minimal()
}
plots <- lapply(cytokines, plot_cytokine_vs_pseudotime)

compute_correlation <- function(cytokine_name) {
  correlation <- cor(AD_SMC@meta.data$Monocle_pseudotime, AD_SMC@meta.data[[cytokine_name]], use = "complete.obs")
  print(paste(cytokine_name, "correlation:", correlation))
}

# Compute correlations and store in a named vector
correlations <- sapply(cytokines, function(cytokine_name) {
  cor(AD_SMC@meta.data$Monocle_pseudotime, AD_SMC@meta.data[[cytokine_name]], use = "complete.obs")
})

# Convert to a data frame
correlation_df <- data.frame(Cytokine = names(correlations),Correlation = correlations)


cytokines <- c("ActivinA", "BDNF", "BMP2", "BMP4", "BMP6", "CD40LG", "CXCL12", 
               "EGF", "FGF2", "GCSF", "GDF11", "GMCSF", "HGF", "IFNG", "IFNL", 
               "IL10", "IL12", "IL13", "IL15", "IL17A", "IL1A", "IL1B", "IL2", 
               "IL21", "IL22", "IL27", "IL3", "IL4", "IL6", "LIF", "LTA", "MCSF", 
               "NO", "OSM", "TGFB1", "TGFB3", "TNFA", "TNFSF12", "TRAIL", 
               "VEGFA", "WNT3A")


compute_correlation <- function(cytokine) {
  # Extract data
  cytokine_data <- as.numeric(AD_SMC@meta.data[[cytokines]])
  
  # Check if data is numeric
  if (!is.numeric(cytokine_data)) {
    stop(paste("Data for", cytokine, "is not numeric."))
  }
  
  # Compute correlations
  AD_SMC_correlation <- cor(cytokine_data, as.numeric(AD_SMC@meta.data$major_cell_type == "IFN-enriched SMC"), use = "complete.obs")
  AD_correlation <- cor(cytokine_data, as.numeric(AD_SMC@meta.data$Group == "AD"), use = "complete.obs")
  pseudotime_correlation <- cor(cytokine_data, as.numeric(AD_SMC@meta.data$Monocle_pseudotime), use = "complete.obs")
  
  # Create a data frame for the cytokine
  data.frame(
    Cytokine = cytokine,
    AD_SMC_correlation = AD_SMC_correlation,
    AD_correlation = AD_correlation,
    pseudotime_correlation = pseudotime_correlation
  )
}


compute_correlation_p_values <- function(cytokine) {
  # Extract data
  cytokine_data <- as.numeric(AD_SMC@meta.data[[cytokine]])
  pseudotime_data <- as.numeric(AD_SMC@meta.data$Monocle_pseudotime)
  Neuroresponsive_data <- as.numeric(AD_SMC@meta.data$major_cell_type == "IFN-enriched SMC") #Neuroresponsive SMC;IFN-enriched SMC
  AD_data <- as.numeric(AD_SMC@meta.data$Group == "AD")
  
  # Check if data is numeric
  if (!is.numeric(cytokine_data)) {
    stop(paste("Data for", cytokine, "is not numeric."))
  }
  
  # Compute correlations and p-values
  cytokine_pseudotime_test <- cor.test(cytokine_data, pseudotime_data, use = "complete.obs")
  Neuroresponsive_test <- cor.test(cytokine_data, Neuroresponsive_data, use = "complete.obs")
  AD_test <- cor.test(cytokine_data, AD_data, use = "complete.obs")
  
  # Create a data frame for the cytokine
  data.frame(
    Cytokine = cytokine,
    Cytokine_Pseudotime_Correlation_r = cytokine_pseudotime_test$estimate,
    Cytokine_Pseudotime_Correlation_P_Value = cytokine_pseudotime_test$p.value,
    Neuroresponsive_Correlation_P_Value = Neuroresponsive_test$p.value,
    AD_Correlation_P_Value = AD_test$p.value
  )
}
correlation_p_values <- lapply(cytokines, compute_correlation_p_values)
correlation_p_values_df <- do.call(rbind, correlation_p_values)
print(correlation_p_values_df)

correlation_p_values_df <- correlation_p_values_df %>%
  mutate(
    Cytokine_Pseudotime_Correlation_P_Value = -log2(Cytokine_Pseudotime_Correlation_P_Value + 1e-10),
    AD_Correlation_P_Value = -log2(AD_Correlation_P_Value + 1e-10),
    Neuroresponsive_Correlation_P_Value = -log2(Neuroresponsive_Correlation_P_Value + 1e-10)
  )
correlation_p_values_df$scaled_Cytokine_Pseudotime_Correlation_r <- scale(correlation_p_values_df$Cytokine_Pseudotime_Correlation_r)
print(correlation_p_values_df)


ggplot(correlation_p_values_df, aes(
  x = Neuroresponsive_Correlation_P_Value,
  y = AD_Correlation_P_Value,
  color = scaled_Cytokine_Pseudotime_Correlation_r,
  size = Cytokine_Pseudotime_Correlation_P_Value)) +
  geom_point() +
  scale_color_gradientn(colors = rev(paletteer::paletteer_c("grDevices::RdYlBu", 20))) +

  labs(
    x = "Neuroresponsive program correlation -log(p_value)",
    y = "AD correlation -log(p_value)",
    color = "Scaled avg.correlation with Pseudotime ",
    size = "Pseudotime Correlation -log(p_value)"
  ) +
  geom_text_repel(aes(label = Cytokine), 
                  size = 3, color = "black",
                  box.padding = 0.5, point.padding = 0.5, segment.color = 'grey50', max.overlaps = 50) +
  geom_vline(xintercept = 1.5, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 1.5, linetype = "dashed", color = "grey") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.background = element_rect(fill = "white", color = "black", size = 1),
    plot.background = element_blank(),  
    legend.position = "right", 
    axis.title = element_text(color = "black"), 
    axis.text = element_text(color = "black"),  
    legend.title = element_text(color = "black"),  
    legend.text = element_text(color = "black"), 
    plot.title = element_text(color = "black") 
  )


pdf("cytokine_plot_mean_expression.pdf", width = 5.5, height = 4, onefile = TRUE)

p1 <- ggplot(correlation_df_labeled, aes(
  x = neglogP_neuro,
  y = neglogP_AD,
  color = scaled_r_pseudotime,
  size = mean_expression
)) +
  geom_point(alpha = 0.7) +
  geom_text_repel(
    aes(label = label_plot),
    size = 3,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey50",
    min.segment.length = 0,
    max.overlaps = Inf,
    na.rm = TRUE
  ) +
  scale_color_gradientn(
    colors = rev(paletteer_c("grDevices::RdYlBu", 20)),
    name = "Scaled r\n(pseudotime)"
  ) +
  scale_size_continuous(
    name = "Mean Expression",
    range = c(1, 8)  
  ) +
  labs(
    x = "-log2(p-value) Neuroresponsive",
    y = "-log2(p-value) AD",
    title = "Cytokine Correlations with Mean Expression"
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

print(p1)
dev.off()


pdf("cytokine_plot_log_expression.pdf", width = 10, height = 8, onefile = TRUE)

p2 <- ggplot(correlation_df_labeled, aes(
  x = neglogP_neuro,
  y = neglogP_AD,
  color = scaled_r_pseudotime,
  size = log_mean_expression
)) +
  geom_point(alpha = 0.7) +
  geom_text_repel(
    aes(label = label_plot),
    size = 3,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey50",
    min.segment.length = 0,
    max.overlaps = Inf,
    na.rm = TRUE
  ) +
  scale_color_gradientn(
    colors = rev(paletteer_c("grDevices::RdYlBu", 20)),
    name = "Scaled r\n(pseudotime)"
  ) +
  scale_size_continuous(
    name = "Log2(Mean Expression+1)",
    range = c(1, 8)
  ) +
  labs(
    x = "-log2(p-value) Neuroresponsive",
    y = "-log2(p-value) AD",
    title = "Cytokine Correlations with Log2(Mean Expression)"
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

print(p2)
dev.off()


pdf("cytokine_plot_scaled_expression.pdf", width = 10, height = 8, onefile = TRUE)

p3 <- ggplot(correlation_df_labeled, aes(
  x = neglogP_neuro,
  y = neglogP_AD,
  color = scaled_r_pseudotime,
  size = scaled_mean_expression
)) +
  geom_point(alpha = 0.7) +
  geom_text_repel(
    aes(label = label_plot),
    size = 3,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey50",
    min.segment.length = 0,
    max.overlaps = Inf,
    na.rm = TRUE
  ) +
  scale_color_gradientn(
    colors = rev(paletteer_c("grDevices::RdYlBu", 20)),
    name = "Scaled r\n(pseudotime)"
  ) +
  scale_size_continuous(
    name = "Scaled Mean\nExpression",
    range = c(1, 8)
  ) +
  labs(
    x = "-log2(p-value) Neuroresponsive",
    y = "-log2(p-value) AD",
    title = "Cytokine Correlations with Scaled Mean Expression"
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

print(p3)
dev.off()


file.exists("cytokine_plot_mean_expression.pdf")
file.exists("cytokine_plot_log_expression.pdf")
file.exists("cytokine_plot_scaled_expression.pdf")




# scMetabolism ------------------------------------------------------------

library(scMetabolism)
Idents(AD_SMC) <- AD_SMC@meta.data$major_cell_type
countexp.Seurat <- sc.metabolism.Seurat(obj = AD_SMC, method = "AUCell", imputation = F, ncores = 2, metabolism.type = "KEGG")

input.pathway <- rownames(countexp.Seurat@assays$METABOLISM$score)[1:85]
DotPlot.metabolism(obj = countexp.Seurat, 
                   pathway = input.pathway, 
                   phenotype = "major_cell_type",
                   norm = "y")

score <- countexp.Seurat@assays$METABOLISM$score
score_change <- score %>% 
  select_all(~str_replace_all(., "\\.", "-"))
identical(colnames(score_change) , rownames(countexp.Seurat@meta.data))
countexp.Seurat@meta.data <- cbind(countexp.Seurat@meta.data,t(score_change) )

countexp.Seurat@meta.data$major_cell_type <- as.character(countexp.Seurat@meta.data$major_cell_type)
df <- countexp.Seurat@meta.data
colnames(countexp.Seurat@meta.data)
head(df)

table(AD_SMC@meta.data$major_cell_type)
selected_celltypes <- c("Contractile SMC","Fibroblast-like SMC","IFN-enriched SMC","Macrophage-like SMC",
                        "Mesenchymal-like SMC", "Neuroresponsive SMC","Secretory SMC",
                        "Stressed SMC","Synthetic SMC")

input.pathway
selected_pathway <- c(
  "Citrate cycle (TCA cycle)",
  "Butanoate metabolism",
  "Valine, leucine and isoleucine degradation",
  "Alanine, aspartate and glutamate metabolism",
  
  "alpha-Linolenic acid metabolism",
  "Ether lipid metabolism",
  "Metabolism of xenobiotics by cytochrome P450",
  
  "Oxidative phosphorylation",
  "Glycine, serine and threonine metabolism",
  "Inositol phosphate metabolism",
  
  "Arginine biosynthesis",
  "Fatty acid biosynthesis",
  "Nitrogen metabolism",
  "Steroid biosynthesis",
  "Vitamin B6 metabolism",
  
  "Lysine degradation",
  "Tyrosine metabolism",
  "Histidine metabolism",
  
  "Fructose and mannose metabolism",
  "Glycolysis / Gluconeogenesis",
  "Nicotinate and nicotinamide metabolism",
  "Pentose phosphate pathway",
  "Porphyrin and chlorophyll metabolism",
  "Thiamine metabolism",
  
  "Glycosaminoglycan biosynthesis - keratan sulfate",
  "Mannose type O-glycan biosynthesis",
  
  "Primary bile acid biosynthesis",
  "Glyoxylate and dicarboxylate metabolism",
  
  "Glycosaminoglycan biosynthesis - chondroitin sulfate / dermatan sulfate",
  "Glycosaminoglycan biosynthesis - heparan sulfate / heparin",
  "Sulfur metabolism",
  "Steroid hormone biosynthesis",
  "Ubiquinone and other terpenoid-quinone biosynthesis"
)



# SMC SASP ----------------------------------------------------------------

senescence_genes_pos <- c("IL6","CXCL8","SERPINE1","IGFBP3","MMP3","MMP1","CXCL1","GDF15","IL1A","IL1B",
                          "CDKN1A","TP53","GADD45A","GADD45B")
senescence_genes_neg <- c("LMNB1","HMGB1")


AD_SMC <- AddModuleScore(AD_SMC, features = list(senescence_genes_pos), name = "Sen_pos")

AD_SMC <- AddModuleScore(AD_SMC, features = list(senescence_genes_neg), name = "Sen_neg")
AD_SMC$sen_score <- AD_SMC$Sen_pos1 - AD_SMC$Sen_neg1

VlnPlot(AD_SMC, "sen_score", group.by="seurat_clusters", pt.size=0.0000000000000000001)

cluster_levels <- levels(Idents(AD_SMC))
mapping <- data.frame(
  cluster = cluster_levels,
  color   = color3[seq_along(cluster_levels)],
  stringsAsFactors = FALSE
)

cluster_to_type <- AD_SMC@meta.data %>%
  dplyr::select(seurat_clusters, major_cell_type) %>%
  dplyr::distinct(seurat_clusters, major_cell_type)
mapping2 <- dplyr::left_join(mapping, cluster_to_type,
                             by = c("cluster" = "seurat_clusters"))
mapping2$major_cell_type <- factor(
  mapping2$major_cell_type,
  levels = levels(AD_SMC$major_cell_type)
)
mapping2 <- mapping2[order(mapping2$major_cell_type), ]
print(mapping2)
color_major <- mapping2$color

color_major <- c("#5586B3", "#45958e", "#B6C796", "#9F77CF", 
                 "#ead4d4", "#C98C9A", "#fad4a6", "#72A77A", 
                 "#cee0e6", "#9ac1f0", "#f5c0bf")
VlnPlot(AD_SMC, "sen_score", pt.size=0, col = color_major)

library(dplyr)


df <- FetchData(AD_SMC, vars = c("sen_score", "major_cell_type"))
target_group <- "Neuroresponsive SMC"
other_groups <- setdiff(unique(df$major_cell_type), target_group)
pvals <- sapply(other_groups, function(g) {
  df_sub <- df %>% filter(major_cell_type %in% c(target_group, g))
  df_sub$group <- ifelse(df_sub$major_cell_type == target_group, target_group, g)
  
  test <- wilcox.test(sen_score ~ group, data = df_sub)
  test$p.value
})


pval_df <- data.frame(
  Comparison = paste(target_group, "vs", other_groups),
  P_value = pvals
)

print(pval_df)


# SMC decoupleR ---------------------------------------------------------------
library(decoupleR)
library(decoupleRdata)
library(dorothea)
library(tibble)
library(OmnipathR)
library(tidyr)
library(dplyr)

expr_mat <- GetAssayData(AD_SMC, slot = "data", assay = "RNA")  # genes  cells
data(dorothea_hs)
dorothea <- dorothea_hs %>% filter(confidence %in% c("A", "B", "C"))


regulon <- dorothea %>%
  select(source = tf, target, likelihood = confidence) %>%
  mutate(likelihood = case_when(
    likelihood == "A" ~ 1,
    likelihood == "B" ~ 0.8,
    likelihood == "C" ~ 0.6,
    TRUE ~ 0
  ))

tf_activity <- run_wmean(mat = expr_mat, network = regulon)
dim(tf_activity) 
head(tf_activity)
head(rownames(tf_activity), 5)

cells_neuro <- colnames(AD_SMC)[AD_SMC$major_cell_type == "Neuroresponsive SMC"]
tf_activity_neuro <- tf_activity %>%
  filter(condition %in% cells_neuro)


tf_mean_neuro <- tf_activity_neuro %>%
  group_by(source) %>%
  summarise(mean_activity = mean(score, na.rm = TRUE)) %>%
  arrange(desc(mean_activity))
dim(tf_mean_neuro)


top_tfs <- tf_mean_neuro %>% head(200)


net <- decoupleR::get_progeny(organism='human', top=500)
net
head(net)
table(net$source)
mat <- as.matrix(AD_SMC@assays$RNA@data)

acts <- run_mlm(
  mat = mat,
  net = net,
  .source = 'source',
  .target = 'target',
  .mor = 'weight',
  minsize = 5
)

tf_matrix <- acts %>%
  tidyr::pivot_wider(id_cols = source, names_from = condition, values_from = score) %>%
  tibble::column_to_rownames('source') %>%
  as.matrix()

identical(colnames(tf_matrix), colnames(AD_SMC))
tf_assay <- Seurat::CreateAssayObject(data = tf_matrix)
AD_SMC[["pathwaysmlm"]] <- tf_assay
DefaultAssay(AD_SMC) <- "pathwaysmlm"
AD_SMC <- Seurat::ScaleData(object = AD_SMC,assay = "pathwaysmlm",features = rownames(tf_matrix))

unique(acts$source)
length(unique(acts$source))


df <- t(as.matrix(AD_SMC@assays$pathwaysmlm@data)) %>%
  as.data.frame() %>%
  dplyr::mutate(major_cell_type = AD_SMC$major_cell_type) %>% 
  tidyr::pivot_longer(cols = -major_cell_type, names_to = "source", values_to = "score") %>%
  dplyr::group_by(major_cell_type, source) %>%
  dplyr::summarise(mean = mean(score), .groups = "drop")


top_acts_mat <- df %>%
  tidyr::pivot_wider(id_cols = major_cell_type, names_from = source, values_from = mean) %>%
  tibble::column_to_rownames(var = "major_cell_type") %>%
  as.matrix()


top_acts_mat_scaled <- scale(top_acts_mat)
top_acts_mat_scaled <- t(scale(t(top_acts_mat)))


desired_order <- c("Contractile SMC", "Fibroblast-like SMC", "IFN-enriched SMC", 
                   "Immune-regulatory SMC", "Macrophage-like SMC", "Mesenchymal-like SMC",
                   "Neuroresponsive SMC", "Proinflammatory SMC", "Secretory SMC",
                   "Stressed SMC", "Synthetic SMC")
top_acts_mat_scaled <- top_acts_mat_scaled[desired_order, ]

colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
colors.use <- colorRampPalette(colors = colors)(100)
colors.use <- colors <- colorRampPalette(c("lightgrey", "#fef0d9", "#fdbb84", "#e34a33"))(100)
colors.use <- colorRampPalette(c("#f7f4f9", "#fcd5ce", "#f08080"))(100)
my_breaks <- c(seq(-1.25, 0, length.out = ceiling(100 / 2) + 1),
               seq(0.05, 1.25, length.out = floor(100 / 2)))



pheatmap::pheatmap(
  mat = top_acts_mat_scaled,
  color = colors,
  breaks = my_breaks,
  border_color = "white",
  cellwidth = 20,
  cellheight = 20,
  treeheight_row = 20,
  treeheight_col = 20,
  main = "TF activity (mean per cell state)",
  cluster_rows = FALSE
)

library(Seurat)
library(decoupleR)
library(dorothea)
library(dplyr)
library(tidyr)


expr_mat <- GetAssayData(AD_SMC, layer = "data", assay = "RNA")
expr_mat <- as.matrix(expr_mat)

cell_types <- AD_SMC$major_cell_type
names(cell_types) <- colnames(AD_SMC)


data(dorothea_hs)
regulon <- dorothea_hs %>%
  filter(confidence %in% c("A", "B")) %>%
  transmute(source = tf, target = target, mor = mor, weight = 1)


cell_type_list <- split(colnames(expr_mat), cell_types)
regulon_by_type <- list()

for (ct in names(cell_type_list)) {
  message(sprintf("Running VIPER on %s (%d cells)", ct, length(cell_type_list[[ct]])))
  cells <- cell_type_list[[ct]]
  sub_expr <- expr_mat[, cells, drop = FALSE]
  
  res <- decoupleR::run_viper(mat = sub_expr, network = regulon, minsize = 4)
  res$group <- ct
  regulon_by_type[[ct]] <- res
}


regulon_all <- bind_rows(regulon_by_type)


rss_df <- regulon_all %>%
  group_by(group, source) %>%
  summarise(score = mean(score, na.rm = TRUE), .groups = "drop")

rss_sub <- rss_df %>%
  filter(group == "Neuroresponsive SMC") %>%
  arrange(desc(score)) %>%
  mutate(rank = row_number(),
         label = ifelse(rank <= 5, source, ""))


reb1_score <- rss_df %>%
  filter(source == "CREB1", group == "Neuroresponsive SMC") %>%
  pull(score)

reb1_rank <- rss_sub %>%
  filter(source == "CREB1") %>%
  pull(rank)

reb1_point <- data.frame(
  rank = ifelse(length(reb1_rank) > 0, reb1_rank, NA),
  score = reb1_score,
  label = "CREB1"
)


ggplot(rss_sub, aes(x = rank, y = score)) +
  geom_point(color = "steelblue", size = 2) +
  geom_text_repel(aes(label = label), max.overlaps = Inf, size = 4) +
  

  geom_point(data = reb1_point, aes(x = rank, y = score), color = "black", size = 3) +
  geom_text_repel(data = reb1_point, aes(x = rank, y = score, label = label),
                  color = "black", size = 4, fontface = "bold") +
  
  labs(
    title = "",
    x = "Rank",
    y = "Regulon Specificity Score"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.6),
    axis.ticks = element_line(color = "black")
  )


# SMC SCENIC --------------------------------------------------------------
library(Seurat)
library(dplyr)


data(dorothea_hs)
regulon <- dorothea_hs %>%
  filter(confidence %in% c("A", "B")) %>%
  transmute(source = tf, target = target, mor = mor, weight = 1)

regulons_list <- split(regulon$target, regulon$source)


expr_mat <- GetAssayData(AD_SMC, layer = "data", assay = "RNA")
expr_mat <- as.matrix(expr_mat)


cellInfo <- data.frame(
  celltype = AD_SMC$major_cell_type,
  row.names = colnames(AD_SMC)
)

library(AUCell)
cells_rankings <- AUCell_buildRankings(expr_mat, nCores=4, plotStats=FALSE)

regulonAUC <- AUCell_calcAUC(regulons_list, cells_rankings)
AUCmat <- getAUC(regulonAUC)
rownames(AUCmat) <- gsub("[(+)]", "", rownames(AUCmat))

library(SCENIC)
cellInfo <- cellInfo[colnames(AUCmat), , drop = FALSE]
rss <- calcRSS(AUC = AUCmat, cellAnnotation = cellInfo$celltype)
rss <- rss[, sort(colnames(rss))]


library(tidyverse)

rss_df <- as.data.frame(rss)
rss_df$regulon <- rownames(rss_df)

rss_long <- rss_df %>%
  pivot_longer(cols = -regulon, names_to = "group", values_to = "rss_score")


rss_sub <- rss_long %>%
  filter(group == "Neuroresponsive SMC") %>%
  arrange(desc(rss_score)) %>%
  mutate(
    rank = row_number(),
    label = ifelse(rank <= 10, regulon, "")
  )


highlight_point <- rss_sub %>%
  filter(regulon == "CREB1")


library(ggrepel)
ggplot(rss_sub, aes(x = rank, y = rss_score)) +
  geom_point(color = "steelblue", size = 2) +
  geom_text_repel(aes(label = label), max.overlaps = Inf, size = 4) +
  geom_point(data = highlight_point, aes(x = rank, y = rss_score),
             color = "black", size = 3) +
  geom_text_repel(data = highlight_point, aes(x = rank, y = rss_score, label = regulon),
                  color = "black", size = 4, fontface = "bold") +
  labs(
    title = "",
    x = "Rank",
    y = "Regulon Specificity Score"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.6),
    axis.ticks = element_line(color = "black")
  )

