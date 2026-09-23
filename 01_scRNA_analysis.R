# Load required packages for the analysis
library(Seurat)
library(dplyr)
library(ggplot2)

# Set working directory to the project root
setwd("/Home/Documents/Conferences/SGBCC africa")

# Load cell-level metadata (sample, disease status, cell type annotations)
metadata <- read.csv("Integrated_Dataset_metadata.csv")

# Quick structural check of the metadata
dim(metadata)
colnames(metadata)
table(metadata$Disease)


# Compute per-sample proportions of each major cell type
major_composition <- metadata %>%
  count(Sample, Disease, major_celltype) %>%
  group_by(Sample) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

# Plot: major cell type composition, primary vs metastatic
ggplot(major_composition, aes(x = Disease, y = proportion, fill = Disease)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1) +
  facet_wrap(~major_celltype, scales = "free_y") +
  scale_fill_manual(values = c(Primary = "blue", Metastasis = "red")) +
  theme_classic() +
  labs(y = "Proportion of sample", x = NULL,
       title = "Major cell type composition by disease status")

ggsave("outputs/figures/celltype_composition_major.png",
       width = 9, height = 6, dpi = 300)

# Inspect minor cell type labels to identify stromal subtypes (CAFs, pericytes)
sort(unique(metadata$minor_celltype))


# Load the integrated Seurat object (Ozmen et al., Zenodo 10.5281/zenodo.13743374)
seurat_obj <- readRDS("Integrated_Dataset_normalized_seurobj.rds")

# Confirm cell identifiers match between the object and the metadata table
mean(metadata$CellID[metadata$major_celltype == "Malignant"] %in% colnames(seurat_obj))

# Subset to malignant cells only
malignant_ids <- metadata$CellID[metadata$major_celltype == "Malignant"]
malignant <- subset(seurat_obj, cells = malignant_ids)

# Load Hallmark gene sets and extract TNFA and interferon-alpha signatures
library(msigdbr)
hallmark <- msigdbr(species = "Homo sapiens", category = "H")
tnfa_genes <- hallmark$gene_symbol[hallmark$gs_name == "HALLMARK_TNFA_SIGNALING_VIA_NFKB"]
ifna_genes <- hallmark$gene_symbol[hallmark$gs_name == "HALLMARK_INTERFERON_ALPHA_RESPONSE"]

# Score each malignant cell for both pathways
malignant <- AddModuleScore(malignant, features = list(tnfa_genes, ifna_genes), name = "HM_")

# Attach disease status to the scored object
malignant$Disease <- metadata$Disease[match(colnames(malignant), metadata$CellID)]

# Plot: TNFA vs interferon-alpha signaling, primary vs metastatic malignant cells
ggplot(malignant@meta.data, aes(x = HM_2, y = HM_1, color = Disease)) +
  geom_point(size = 0.3, alpha = 0.3) +
  scale_color_manual(values = c(Primary = "blue", Metastasis = "red")) +
  theme_classic() +
  labs(x = "Interferon-alpha response score", y = "TNFα/NF-κB signaling score",
       title = "Malignant cell signaling states by disease status")

ggsave("outputs/figures/malignant_TNFa_vs_IFNa.png",
       width = 7, height = 6, dpi = 300)


# Load the integrated Seurat object (Ozmen et al., Zenodo 10.5281/zenodo.13743374)
seurat_obj <- readRDS("Integrated_Dataset_normalized_seurobj.rds")

# Confirm cell identifiers match between the object and the metadata table
mean(metadata$CellID[metadata$major_celltype == "Malignant"] %in% colnames(seurat_obj))

# Subset to malignant cells only
malignant_ids <- metadata$CellID[metadata$major_celltype == "Malignant"]
malignant <- subset(seurat_obj, cells = malignant_ids)

# Load Hallmark gene sets and extract TNFA and interferon-alpha signatures
library(msigdbr)
hallmark <- msigdbr(species = "Homo sapiens", category = "H")
tnfa_genes <- hallmark$gene_symbol[hallmark$gs_name == "HALLMARK_TNFA_SIGNALING_VIA_NFKB"]
ifna_genes <- hallmark$gene_symbol[hallmark$gs_name == "HALLMARK_INTERFERON_ALPHA_RESPONSE"]

# Score each malignant cell for both pathways
malignant <- AddModuleScore(malignant, features = list(tnfa_genes, ifna_genes), name = "HM_")

# Attach disease status to the scored object
malignant$Disease <- metadata$Disease[match(colnames(malignant), metadata$CellID)]

# Plot: TNFA vs interferon-alpha signaling, primary vs metastatic malignant cells
ggplot(malignant@meta.data, aes(x = HM_2, y = HM_1, color = Disease)) +
  geom_point(size = 0.3, alpha = 0.3) +
  scale_color_manual(values = c(Primary = "blue", Metastasis = "red")) +
  theme_classic() +
  labs(x = "Interferon-alpha response score", y = "TNFα/NF-κB signaling score",
       title = "Malignant cell signaling states by disease status")

ggsave("outputs/figures/malignant_TNFa_vs_IFNa.png",
       width = 7, height = 6, dpi = 300)


# Define CAF and pericyte subtypes present in the metadata
stromal_subtypes <- c("F01 iCAF", "F02 mCAF", "F03 apCAF", "F04 Pericyte")

# Compute per-sample proportion of each stromal subtype among all stromal cells
stromal_composition <- metadata %>%
  filter(minor_celltype %in% stromal_subtypes | grepl("^F0", minor_celltype)) %>%
  count(Sample, Disease, minor_celltype) %>%
  group_by(Sample) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup() %>%
  filter(minor_celltype %in% stromal_subtypes)

# Plot: stromal subtype composition, primary vs metastatic
ggplot(stromal_composition, aes(x = Disease, y = proportion, fill = Disease)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1) +
  facet_wrap(~minor_celltype, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c(Primary = "blue", Metastasis = "red")) +
  theme_classic() +
  labs(y = "Proportion of stromal cells", x = NULL,
       title = "Fibroblast subtype composition by disease status")

ggsave("outputs/figures/stromal_composition.png",
       width = 10, height = 5, dpi = 300)


# Define macrophage/monocyte subtypes of interest (matching key myeloid states)
myeloid_subtypes <- c("M18 Macrophage-FOLR2", "M09 Macrophage-CX3CR",
                      "M08 Macrophage-CCL2", "M07 Macrophage-SPP1",
                      "M21 Macrophage-MMP9", "M04 Monocyte Intermediate-FCN1")

# Compute per-sample proportion of each myeloid subtype among all myeloid cells
myeloid_composition <- metadata %>%
  filter(grepl("^M", minor_celltype)) %>%
  count(Sample, Disease, minor_celltype) %>%
  group_by(Sample) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup() %>%
  filter(minor_celltype %in% myeloid_subtypes)

# Plot: myeloid subtype composition, primary vs metastatic
ggplot(myeloid_composition, aes(x = Disease, y = proportion, fill = Disease)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1) +
  facet_wrap(~minor_celltype, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c(Primary = "blue", Metastasis = "red")) +
  theme_classic() +
  labs(y = "Proportion of myeloid cells", x = NULL,
       title = "Macrophage/monocyte subtype composition by disease status")

ggsave("outputs/figures/myeloid_composition.png",
       width = 10, height = 6, dpi = 300)


# Define T cell subtypes of interest (matching key immune states)
tcell_subtypes <- c("T05 Tregulatory-FOXP3", "T09 CD8 Texhausted-CXCL13",
                    "T01 CD4 Tcentralmemory-LMNA", "T07 Tactivated-IFI6")

# Compute per-sample proportion of each T cell subtype among all T cells
tcell_composition <- metadata %>%
  filter(grepl("^T", minor_celltype)) %>%
  count(Sample, Disease, minor_celltype) %>%
  group_by(Sample) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup() %>%
  filter(minor_celltype %in% tcell_subtypes)

# Plot: T cell subtype composition, primary vs metastatic
ggplot(tcell_composition, aes(x = Disease, y = proportion, fill = Disease)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1) +
  facet_wrap(~minor_celltype, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c(Primary = "blue", Metastasis = "red")) +
  theme_classic() +
  labs(y = "Proportion of T cells", x = NULL,
       title = "T cell subtype composition by disease status")

ggsave("outputs/figures/tcell_composition.png",
       width = 10, height = 5, dpi = 300)

save.image("outputs/scRNA_analysis_workspace.RData")
saveRDS(malignant, "outputs/malignant_scored.rds")