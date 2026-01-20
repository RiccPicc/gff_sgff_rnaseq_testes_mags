rm(list=ls(all=TRUE))
library(openxlsx)
library(taxonomizr)
library(tidyverse)
library(ggplot2)
library(phyloseq)
library(decontam)
library(ggpubr)
library(gridExtra)
library("pheatmap")
library(DESeq2)
library(writexl)

load("tax.ranks.RData") # row 19 performed in this RData
load("12.01.2026.RData")

#### phyloseq object creation output ####
# This function returns a phyloseq object given the output from hgtseq
ps_from_hgtseq <- function(hits, sampledata, taxa="Bacteria", decontamination=T, min_reads=10, tax.ranks=NULL){
  print(paste(Sys.time(), "- Starting creating phyloseq object from hgtseq output"))
  
  if (is.null(tax.ranks)){
    prepareDatabase('accessionTaxa.sql',types='prot',protocol='http')
    tax.ranks <<- getTaxonomy(hits$taxID,'accessionTaxa.sql') # global
  } 
  
  # filtering taxa ids from taxonomic ranks database
  taxa.tab <- filter_taxa_from_db(taxa, tax.ranks)
  taxmat <- tax_table_to_mat(taxa.tab)
  
  # make otu table from hgtseq table
  otu_tab <- otu_from_hgt(hits)
  otumat <- otu_table_to_mat(otu_tab, taxa.tab)
  
  # store samples information into a dataframe
  sd <- as.data.frame(sampledata)
  sd$quant_reading <- colSums(otu_tab)[-1] # add columns containing info about the quantity of reads per sample belonging to Bacteria according to Kraken2
  rownames(sd) <- sd$sample
  sd <- sample_data(sd)
  
  # phyloseq object
  ps <- phyloseq(otumat, taxmat, sd)
  
  # decontamination
  if (decontamination){
    ps <- decontam_data(ps)
    print(paste(Sys.time(), "- Decontamination step finished"))
  }
  
  ps <- prune_taxa(taxa_sums(ps) > min_reads, ps)
  
  print(paste(Sys.time(), "- phyloseq object created with success!"))
  return(ps)
}


# This fuction filters from the database a given superkingdom taxa (i.e. Bacteria, Archaea, Virus, Eukaryota)
filter_taxa_from_db <- function(taxa, tax.ranks){
  tax.ID <- as.data.frame(rownames(tax.ranks))
  colnames(tax.ID) <- c("taxID")
  tax.copy <- as_tibble(cbind(tax.ID, tax.ranks)) %>% 
    mutate(taxID = as.numeric(taxID)) %>%
    unique()
  taxa.tab <- tax.copy %>%
    filter(superkingdom == taxa, !is.na(phylum)) 
  
  return(taxa.tab)
}

# This function creates the otu table of available species found with the hgtseq run
otu_from_hgt <- function(hits){
  otu_tab <- hits %>%
    as_tibble() %>%
    filter(readtype != "single_mapped") %>%
    select(sample, taxID, read_count) %>%
    group_by(sample, taxID) %>%
    summarise(read_count = sum(read_count)) %>%
    pivot_wider(names_from = sample, values_from = read_count, values_fill = 0)
  
  return(otu_tab)
}

# This function converts the taxa table into a taxa matrix (for phyloseq object)
tax_table_to_mat <- function(taxa.tab){
  taxmat <- taxa.tab %>% as.data.frame()
  rownames(taxmat) <- paste0("OTU", as.numeric(taxa.tab$taxID))
  taxmat <- taxmat[,-1]
  colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  taxmat <- tax_table(as.matrix(taxmat))
  
  return(taxmat)
}

# This function converts the otu table into a otu matrix (for phyloseq object)
otu_table_to_mat <- function(otu_tab, taxa.tab){
  otumat <- otu_tab %>% filter(taxID %in% taxa.tab$taxID) %>% as.data.frame()
  rownames(otumat) <- paste0("OTU", otumat$taxID) # rownames = OTU + taxID NCBI
  otumat <- as.matrix(otumat[,-1])
  otumat <- otu_table(otumat, taxa_are_rows = TRUE)
  
  return(otumat)
}

# This function create a ps object removing contaminants by frequency (method usable when no negative controls are available)
decontam_data <- function(ps){
  print("Decontamination step may take a while")
  print(paste(Sys.time(), "- Starting decontamination by frequency"))
  contamdf.freq <- isContaminant(ps, method="frequency", conc="quant_reading")
  ps <- prune_taxa(!contamdf.freq$contaminant, ps)
  
  return(ps)
}

get_cutoff_bacteria <- function(ps, threshold=10){
  return(as.data.frame(otu_table(ps)) %>%
           rownames_to_column("Genus") %>%
           as_tibble() %>%
           rowwise() %>%
           mutate(mean_abundance = mean(c_across(-Genus))) %>%
           ungroup() %>%
           filter(mean_abundance > threshold))
}

get_number_bacteria <- function(ps, threshold=10){
  return( get_cutoff_bacteria(ps, threshold) %>%
            summarise(num_genera = n_distinct(Genus)) %>%
            pull(num_genera))
}

get_tissue_bacteria <- function(ps, threshold=10){
  cutoff_bacteria <- get_cutoff_bacteria(ps, threshold) %>% select(Genus, mean_abundance)
  genera_gram <- gram_annot(ps)
  joined_genera <- left_join(cutoff_bacteria, genera_gram, by = join_by(Genus == tax.ID)) %>% 
    dplyr::rename(OTU.ID = Genus, Genus = Genus.y) %>%
    relocate(Genus, .after = OTU.ID)
  joined_genera <- joined_genera %>% relocate(Gram , .after = Genus)
  return(joined_genera)   
}

differential_analysis <- function(ps){
  diagdds = phyloseq_to_deseq2(ps, ~ infection_status)
  diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
  res = results(diagdds, cooksCutoff = FALSE)
  alpha = 0.05
  sigtab = res[which(res$padj < alpha), ]
  if (dim(sigtab)[1] == 0){
    warning("No significant differences found")
    return(sigtab)
    }
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
  head(sigtab)
  
  df <- as.data.frame(colData(diagdds)[,c("infection_status")])
  counter <- counts(diagdds, normalized = T)[rownames(sigtab),]
  rownames(counter) <- sigtab$Genus
  counts_perc <- t(apply(counter, 1, function(x) x / max(x)))
  counts_perc[counts_perc == -Inf] <- 0
  counts_perc <- counts_perc[rowSums(counter > 0.1) >= 3, ]
  up_presence <- sigtab %>% filter(Genus %in% rownames(counts_perc), log2FoldChange>0)
  low_presence <- sigtab %>% filter(Genus %in% rownames(counts_perc), log2FoldChange<0)
  
  return(list(df, counter, counts_perc, up_presence, low_presence,diagdds))
}

#### plotting functions ####

# This function creates a table where OTUs are annotated by Gram coloration
gram_annot <- function(ps, is.tab=F){
  if (!is.tab){
    full.taxa.tab <- inner_join(tax_table(ps) %>% as.data.frame() %>% mutate(tax.ID = rownames(.)), 
                                otu_table(ps) %>% as.data.frame() %>% mutate(tax.ID = rownames(.)), 
                                by = "tax.ID")
  } else {
    full.taxa.tab <- ps
  }
  
  gram.tab <- full.taxa.tab %>% 
    mutate(
      Gram = case_when(
        Domain != "Bacteria" ~ "No bacteria",
        Phylum %in% c("Mycoplasmatota") ~ "Gram-neutral",
        Phylum %in% c("Bacillota", "Cyanobacteriota", "Chloroflexota", "Actynomicetota", "Deinococcus") ~ "Gram-positive",
        TRUE ~ "Gram-negative"
      )
    ) %>% # added gram coloration based on phylum
    select(-Domain, -Phylum, -Class, -Order, -Family, -Species)
  
  return(gram.tab)
}

# This function return the plot based on Gram coloration
plot_gram <- function(counts, meta, gram_tab,
                      gram = "Gram-negative",
                      stat.m = "wilcox.test",
                      div.m = "Shannon",
                      min.taxa = 10) {
  
  genera_keep <- gram_tab %>% filter(Gram == gram) %>% pull(taxID)
  counts_gram <- counts[rownames(counts) %in% genera_keep, ]
  
  counts_gram <- counts_gram[rowSums(counts_gram) > min.taxa, ]
  
  shannon <- diversity(t(counts_gram), index = "shannon")
  
  df_div <- meta %>%
    mutate(sampleID = rownames(meta)) %>%
    mutate(Shannon = shannon[match(sampleID, names(shannon))])
  
  p <- ggplot(df_div, aes(x = infection_status, y = Shannon)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.6) +
    facet_wrap(~ tissue, scales = "free_y") +
    theme_classic() +
    ggpubr::stat_compare_means(
      comparisons = list(c("neg", "pos")),
      method = stat.m,
      label = "p.signif"
    ) +
    labs(title = paste(div.m, "diversity of", gram, "bacteria"),
         x = "Infection status",
         y = paste(div.m, "index"))
  
  return(p)
}

# This function return the plot showing bacteria differences
plot_genus <- function(counts, meta, genus = "Escherichia",
                       stat.m = "wilcox.test", new.labs = NULL) {
  
  dat_long <- as.data.frame(counts) %>%
    rownames_to_column("Genus") %>%
    pivot_longer(-Genus, names_to = "SampleID", values_to = "Abundance")
  
  dat_merged <- dat_long %>%
    left_join(meta, by = c("SampleID" = "sample"))
  
  dat_filtered <- dat_merged %>% filter(Genus == genus)
  
  if (!is.null(new.labs)) {
    if (length(new.labs) == length(unique(dat_filtered$tissue))) {
      old_levels <- levels(as.factor(dat_filtered$tissue))
      message("Renaming ", paste(old_levels, collapse = ", "),
              " into ", paste(new.labs, collapse = ", "))
      dat_filtered$tissue <- factor(dat_filtered$tissue, labels = new.labs)
    } else {
      warning("Number of new labels does not match number of tissue levels.")
    }
  }
  
  p <- ggplot(dat_filtered, aes(x = infection_status, y = Abundance, fill = infection_status)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.6) +
    facet_wrap(~ tissue, scales = "free_y") +
    theme_classic() +
    scale_fill_manual(values = c("lightblue", "salmon")) +
    ggpubr::stat_compare_means(
      comparisons = list(c("neg", "pos")),
      method = stat.m,
      label = "p.signif"
    ) +
    labs(title = paste("Normalized abundance of", genus),
         x = "Infection status",
         y = "Normalized counts")
  
  return(p)
}

#### main ####
hits <- read.xlsx("glossina_classified_unmapped_read_counts_by-sample.xlsx") # mac

# make sample data table
sampledata <- hits %>%
  as_tibble() %>%
  dplyr::select(sample) %>%
  unique() %>% 
  mutate(
    tissue = if_else(str_detect(sample, "GLANDS"), "mags", "testes"),
    infection_status = if_else(str_detect(sample, "POS"), "pos", "neg"),
    sample_id = str_extract(sample, "S\\d+")
  )

# create ps object from hits
# hits <- read.xlsx("glossina_classified_unmapped_read_counts_by-sample.xlsx") # windows
ps <- ps_from_hgtseq(hits, sampledata, tax.ranks=tax.ranks)

# numbers of bacteria
ps_testes <- prune_samples(sample_data(ps)$tissue == "testes", ps)
ps_testes_genus <- tax_glom(ps_testes, taxrank = "Genus")
n_genera_testes <- get_number_bacteria(ps_testes_genus)
testes_bacteria <- get_tissue_bacteria(ps_testes_genus)
write_xlsx(testes_bacteria %>% select(-Gram), path="testes_genera.xlsx")

ps_mags <- prune_samples(sample_data(ps)$tissue == "mags", ps)
ps_mags_genus <- tax_glom(ps_mags, taxrank = "Genus")
n_genera_mags <- get_number_bacteria(ps_mags_genus)
mags_bacteria <- get_tissue_bacteria(ps_mags_genus)
write_xlsx(mags_bacteria %>% select(-Gram), path="mags_genera.xlsx")

testes_bacteria$tissue <- "testes"
mags_bacteria$tissue <- "MAGs"
testes_bacteria_full_tab <- testes_bacteria %>% select(OTU.ID, Genus, Gram, mean_abundance, tissue)
mags_bacteria_full_tab <- mags_bacteria %>% select(OTU.ID, Genus, Gram, mean_abundance, tissue)
full_tab <- rbind(testes_bacteria_full_tab, mags_bacteria_full_tab)

wb <- createWorkbook()
addWorksheet(wb, "Genera_in_tissues")
writeData(wb, "Genera_in_tissues", full_tab)
addWorksheet(wb, "Testes")
writeData(wb, "Testes", testes_bacteria)
addWorksheet(wb, "MAGs")
writeData(wb, "MAGs", mags_bacteria)
saveWorkbook(wb, "H:/Il mio Drive/Articoli/G.fuscipes_Spiroplasma_RNAseq/Tables/Table_S3.xlsx", overwrite = F)

# deseq2 part testes
results_testes <- differential_analysis(ps_testes_genus)
df <- results_testes[[1]]
testes_counts <- results_testes[[2]]
counts_perc <- results_testes[[3]]
up_presence <- gram_annot(results_testes[[4]], is.tab=T)
low_presence <- gram_annot(results_testes[[5]], is.tab=T)
write.xlsx(rbind(up_presence, low_presence) %>% arrange(desc(baseMean)),
           "H:/Il mio Drive/Articoli/G.fuscipes_Spiroplasma_RNAseq/Tables/Table_2.xlsx",
          overwrite = F)
dds <- results_testes[[6]]
rld <- rlog(dds, blind=FALSE)
pcaData <- plotPCA(rld, intgroup=c("infection_status"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=infection_status, shape=infection_status)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

outliers <- c("TESTISPOS_S17_Sp_Testes", "TESTISPOS_S18_Sp_Testes")
pcaData_filtered <- pcaData[!(rownames(pcaData) %in% outliers), ]
print(paste("Original sample count:", nrow(pcaData)))
print(paste("Filtered sample count:", nrow(pcaData_filtered)))
percentVar <- round(100 * attr(pcaData_filtered, "percentVar"))
ggplot(pcaData_filtered, aes(PC1, PC2, color=infection_status, shape=infection_status)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

results_testes <- differential_analysis(ps_testes_genus)
df <- results_testes[[1]]
testes_counts <- results_testes[[2]]
counts_perc <- results_testes[[3]]
up_presence <- gram_annot(results_testes[[4]], is.tab=T)
low_presence <- gram_annot(results_testes[[5]], is.tab=T)
dds <- results_testes[[6]]
rld <- rlog(dds, blind=FALSE)
pcaData <- plotPCA(rld, intgroup=c("infection_status"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=infection_status, shape=infection_status)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()


#remake deseq2 after removing outliers
results_testes <- differential_analysis(subset_samples(ps_testes_genus, !(sample_names(ps_testes_genus) %in% outliers)))
df <- results_testes[[1]]
testes_counts <- results_testes[[2]]
counts_perc <- results_testes[[3]]
up_presence <- gram_annot(results_testes[[4]], is.tab=T)
low_presence <- gram_annot(results_testes[[5]], is.tab=T)
write.xlsx(rbind(up_presence, low_presence) %>% arrange(desc(baseMean)),
           'table_2.xlsx',
           overwrite = F)


# deseq2 part mags
results_mags <- differential_analysis(ps_mags_genus)

png("heatmaps.png", width = 6, height = 6, units = "in", res = 300)
pheatmap(counts, cluster_rows=FALSE, show_rownames=T,
         cluster_cols=FALSE, show_colnames=F, annotation_col=df)
dev.off()

wol <- plot_genus(testes_counts, sampledata, genus="Wolbachia", new.labs=c("MAGs", "Testes")) +
  labs(fill='Infection Status') +
  theme(legend.position = 'right',
        text = element_text(size = 16),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.spacing.y = unit(.5, 'cm')) +
  ylab("Reads count") +
  scale_fill_discrete(labels=c('Spiroplasma\nNegative', 'Spiroplasma\nPositive')) +
  guides(fill = guide_legend(byrow = TRUE))
rom <- plot_genus(ps, genus="Romboutsia")
# print(plot_genus(ps, genus="Wigglesworthia", new.labs=c("MAGs", "Testes")))
# plot_genus(ps, genus="Sodalis")
# plot_genus(ps, genus="Spiroplasma")
ggsave("spiro_wolbachia.png", plot = wol, width = 6, height = 4, units = "in", dpi = 300)
ggsave("spiro_romboutsia.png", plot = rom, width = 6, height = 6, units = "in", dpi = 300)
ggsave("spiro_wiggle.png", plot = wol, width = 6, height = 6, units = "in", dpi = 300)
ggsave("spiro_sodalis.png", plot = wol, width = 6, height = 6, units = "in", dpi = 300)

#### analysis performed on Son 2021 ####
sample_data <- read.csv("H:/Il mio Drive/Spiroplasma/microbiota/son21/map_sample_srr.tsv", sep="\t")
kraken2_res <- read.csv("H:/Il mio Drive/Spiroplasma/microbiota/son21/results/taxprofiler/kraken2_db1.tsv", sep="\t")
colnames(kraken2_res) <- gsub("_db1\\.kraken2\\.kraken2\\.report", "", colnames(kraken2_res))
wolbachia_data <- kraken2_res %>%
  mutate(across(-taxonomy_id, ~ . / sum(.) * 100)) %>% 
  filter(taxonomy_id == 953)
matched_values <- sapply(sample_data$SRR, function(srr) {
  if (srr %in% colnames(wolbachia_data)[-1]) {
    wolbachia_data[1, which(colnames(wolbachia_data) == srr)]
  } else {
    NA  # Assign NA if no match is found
  }
})
sample_data$Wolbachia_Count <- as.numeric(matched_values)
wilcox.test(
  Wolbachia_Count ~ Status, 
  data = sample_data,
  exact = FALSE # Use approximate p-value for larger samples
)


#### plots from ps ####
otus <- otu_table(ps) %>% as.data.frame()
otus$otu <- rownames(otus)
taxatab <-tax_table(ps) %>% as.data.frame()
taxatab$otu <- rownames(taxatab)
fulldata <- otus  %>% as_tibble() %>%
  left_join(taxatab  %>% as_tibble(),by = join_by(otu))

wiggle <- fulldata %>% filter(Genus=="Wigglesworthia") %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE))
wolb <- fulldata %>% filter(Genus=="Wolbachia") %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE))
sodalis <- fulldata %>% filter(Genus=="Sodalis") %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) 
spiro <- fulldata %>% filter(Genus=="Spiroplasma") %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) 
  
dat_long <- as.data.frame(wiggle) %>%
  rownames_to_column("Genus") %>%
  pivot_longer(-Genus, names_to = "SampleID", values_to = "Abundance")

dat_merged <- dat_long %>%
  left_join(sampledata, by = c("SampleID" = "sample"))

new.labs=c("MAGs", "Testes")

if (!is.null(new.labs)) {
  if (length(new.labs) == length(unique(dat_merged$tissue))) {
    old_levels <- levels(as.factor(dat_merged$tissue))
    message("Renaming ", paste(old_levels, collapse = ", "),
            " into ", paste(new.labs, collapse = ", "))
    dat_merged$tissue <- factor(dat_merged$tissue, labels = new.labs)
  } else {
    warning("Number of new labels does not match number of tissue levels.")
  }
}
genus = "Wigglesworthia"
p <- ggplot(dat_merged %>% filter(SampleID != "TESTISPOS_S17_Sp_Testes"), aes(x = infection_status, y = Abundance, fill = infection_status)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  facet_wrap(~ tissue, scales = "free_y") +
  theme_classic() +
  scale_fill_manual(values = c("lightblue", "salmon")) +
  ggpubr::stat_compare_means(
    comparisons = list(c("neg", "pos")),
    method = "wilcox.test",
    label = "p.signif"
  ) +
  labs(title = paste("Normalized abundance of", genus),
       x = "Infection status",
       y = "Normalized counts")

return(p)

### heatmap of diff bacteria
dds <- phyloseq_to_deseq2(ps_testes_genus, ~ infection_status)
keep <- rowSums(counts(dds)) >= 10
ddsMF <- dds[keep,]
levels(ddsMF$infection_status)
design(ddsMF) <- formula(~ infection_status)
ddsMF <- DESeq(ddsMF, test="Wald", fitType="parametric")
resMF <- results(ddsMF, cooksCutoff = TRUE)
resultsNames(ddsMF)
vsd <- varianceStabilizingTransformation(ddsMF, blind=FALSE) # vst normalization

resMFType <- resMF[order(resMF$pvalue),] # order by adjusted p-value
alpha = 0.05
sigtab = resMFType[which(resMFType$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps_testes_genus)[rownames(sigtab), ], "matrix"))
dim(sigtab)
sigtab <- sigtab[-13] # removing doubled species column

tax.sig <- as.data.frame(cbind(rownames(sigtab),tax_table(ps_testes_genus)[rownames(sigtab),]))[,-8] # get tax info and removing subspecies

# log2fold change plot
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))

df <- as.data.frame(colData(ddsMF)[,c("infection_status")])
df$infection_status <- factor(df$`colData(ddsMF)[, c("infection_status")]`, levels=c('neg','pos'))
df <- df[order(factor(df$infection_status, levels=c('neg','pos'))),]
df <- df %>% dplyr::select(infection_status)

sig_rows <- sigtab %>% 
  mutate(genus = ifelse(!is.na(Genus), as.character(Genus), as.character(Family))) %>% 
  mutate(genus = ifelse(!is.na(genus), genus, as.character(Order))) %>%
  mutate(genus = ifelse(!is.na(genus), genus, as.character(Class))) %>%
  select(genus)
rownames(sig_rows) <- rownames(assay(vsd)[c(rownames(sigtab)),])
my_colors <- c("neg" = "steelblue", "pos" = "#D95F039A")
plot_deseq <- pheatmap(assay(vsd)[c(rownames(sigtab)),], cluster_rows=FALSE, show_rownames=T,
                       cluster_cols=FALSE, annotation_col=df, annotation_colors=list(infection_status=my_colors),
                       labels_row=sig_rows$genus)
png(paste('pheatmap_DESeq2_sig_Cooks.png',sep=''), width=6.5,units="in", height=12, res=1200)
print(plot_deseq)
dev.off()
