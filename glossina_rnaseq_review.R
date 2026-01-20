rm(list=ls(all=T))
library(DESeq2)
library(tidyverse)
library(writexl)
library(Biostrings)
library("pheatmap")
library(BaseSet)
library(GOfuncR)
library(topGO)
library(enrichplot)
library("RColorBrewer")
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library("grid")
library("gridExtra")
library(cowplot)
library(readxl)
library(biomaRt)
library(ggplotify)
library(patchwork)

load("negative_analysis.RData")
quartzFonts(helvetica = c("Helvetica Neue Light", 
                          "Helvetica Neue Bold", "Helvetica Neue Light Italic", 
                          "Helvetica Neue Bold Italic"))

# reads counts info
reads<- read_csv('fastqc_sequence_counts_plot.csv')
reads <- reads %>% 
  mutate(tissue = ifelse(grepl("Sp_acc_gl", Category), "MAGs", "Testes"), 
         sampleID = as.numeric(str_extract(Category, "(?<=^S)\\d+"))) %>%
  mutate(infection = ifelse(sampleID<=12, "Neg", "Pos")) %>%
  group_by(sampleID, tissue, infection) %>% 
  summarise(unique_reads = sum(`Unique Reads`, na.rm = TRUE),
            duplicate_reads = sum(`Duplicate Reads`, na.rm = TRUE)) %>%
  mutate(total_reads = sum(unique_reads, duplicate_reads),
         tissue = as.factor(tissue), 
         infection = as.factor(infection)) 
wilcox.test(reads$total_reads ~ reads$infection)
wilcox.test(reads$total_reads ~ reads$tissue)
write_xlsx(reads, path = 'table_s1.xlsx')

# data preparation
                                                                                                                                            par(family="helvetica")
sampledata <- tibble(
  sample = as_tibble(colData(dds))$sample,
  condition = sub("(^.+)(_S.+_Sp_.*)$", "\\1", as_tibble(colData(dds))$sample)
) %>%
  mutate(
    spiro = case_when(
      grepl("POS", condition) ~ "positive",
      grepl("NEG", condition) ~ "negative",
      TRUE ~ NA
    ),
    tissue = case_when(
      grepl("GLAND", condition) ~ "glands",
      grepl("TEST", condition) ~ "testis",
      TRUE ~ NA
    )
  ) %>%
  mutate(
    condition = factor(condition, levels = c("GLANDSNEG","TESTISNEG", "GLANDSPOS", "TESTISPOS")),
    spiro = factor(spiro, levels = c("negative", "positive")),
    tissue = factor(tissue, levels = c("glands", "testis"))
  )

####### JUST NEGATIVES 
negdata <- sampledata %>% 
  filter(spiro=="negative")

idx <- c(7:12,19:24) # campioni da eliminare 
dds_neg <- dds[,-idx] # mantengo i campioni dei neg
dds_neg <- DESeqDataSetFromMatrix(countData=assay(dds_neg), colData=negdata, design = ~ tissue)


######## filtering on number of samples

keep_neg <- rowSums(counts(dds_neg) >= 10) >= 6
dds_neg <- dds_neg[keep_neg,]

#### DE analysis for both

dds_neg <- DESeq(dds_neg)

################# results

res_neg = results(dds_neg)
res_neg = res_neg[order(res_neg$padj),]
res_neg_export = tibble(
  bind_cols(
    gene = row.names(res_neg),
    as_tibble(res_neg)
  )
)

dim(res_neg_export%>%
      filter(padj < 0.05 & abs(log2FoldChange) > 2))

testes.biased <- res_neg_export %>% filter(padj < 0.05 & log2FoldChange > 2) 
mags.biased <- res_neg_export %>% filter(padj < 0.05 & -log2FoldChange > 2) 
mean(testes.biased$baseMean)
median(testes.biased$baseMean)
dim(testes.biased)[1]
mean(mags.biased$baseMean)
median(mags.biased$baseMean)
dim(mags.biased)[1]

write_tsv(res_neg_export, file = "glossina_neg-design_de-results.tsv")
library(writexl)
write_xlsx(res_neg_export, path = "glossina_neg-design_de-results.xlsx")

#### shrink for visualisation

resLFC_neg <- lfcShrink(dds_neg, coef="tissue_testis_vs_glands", type="apeglm")
resLFC_neg <- na.omit(resLFC_neg[order(abs(resLFC_neg$log2FoldChange), decreasing=T),]) 
resLFC_neg_export <- tibble(
  bind_cols(
    gene = row.names(resLFC_neg),
    as_tibble(resLFC_neg)
  )
)

dim(resLFC_neg_export%>%
      filter(padj < 0.05 & abs(log2FoldChange) > 2))

resLFC_neg_plt <- resLFC_neg_export %>% 
  dplyr::mutate(l2FC = cut(log2FoldChange, 
                breaks=c(-Inf,-2,2,Inf), 
                right = FALSE, 
                labels = c("MAGs-biased genes", 
                           "Neutral genes",
                           "Testes-biased genes"))) %>%
  dplyr::mutate(l2FC = as.factor(ifelse(padj>0.05,"Neutral genes", as.character(l2FC))))

resLFC_neg_plt$l2FC <- factor(resLFC_neg_plt$l2FC, levels=c("MAGs-biased genes", 
                                                            "Neutral genes",
                                                            "Testes-biased genes"))

##### volcano plot
png('volcano-plot_neg-testis-glands.png', width=13,units="in", height=10, res=300)
volcano <- ggplot(resLFC_neg_plt, aes(x=log2FoldChange, y=-log10(padj))) + 
  geom_point(aes(color = l2FC, size=abs(log2FoldChange)),  alpha = 0.3) + 
  geom_hline(yintercept=-log10(0.05), linetype = "dashed", alpha = 0.8) +
  geom_vline(xintercept=c(2, -2),linetype = "dashed", alpha = 0.8) +
  # geom_hline(yintercept=350, linetype = "dashed", alpha = 0.8, color = NA) +
  scale_color_manual(values = c("forestgreen","grey10","goldenrod"))  + 
  labs(color='Tissue bias', size='|l2FC|') +
  theme_classic() +
  theme(legend.position = 'right',
        text = element_text(size = 16)) + 
  grids(linetype = "dashed")  + 
  guides(colour = guide_legend(override.aes = list(size=10))) +
  ylab(bquote("-" ~ log[10] ~ italic(P)))+
  xlab(bquote(log[2] ~ "fold change (l2FC)"))
print(volcano)
dev.off()

## the top gene is the same in the 2 designs
plotCounts(dds_neg, gene=which.min(res_neg$padj), intgroup="tissue")

ntd <- normTransform(dds_neg)
df <- as.data.frame(colData(dds_neg)[,"tissue"])
rownames(df) <- colnames(ntd)
colnames(df) <- "Tissue"
df$Tissue <- as.factor(ifelse(df$Tissue == "glands", "MAGs", "Testes"))
select <- order(rowMeans(counts(dds_neg,normalized=TRUE)),
                decreasing=TRUE)[1:50]
selected_genes <- assay(ntd)[select,]
anno_tab <- read_xlsx("./PostDoc/Spiroplasma/degs_tables/glossina_annotation-full_table.xlsx")
anno_tab_2 <- read.csv("./PostDoc/Spiroplasma/Gff2Dmel_proteins/Gff2Dmel.txt", sep="\t",header=F)
anno_tab <- anno_tab %>% dplyr::filter(grepl("LOC",Name)) %>% dplyr::select(gene_symbol, Name)
gff2dmel <- merge(anno_tab_2, anno_tab, by.x = "V1", by.y = "gene_symbol") %>% 
  dplyr::rename(gff_symbol = V1, gff_id = Name, dmel_symbol = V2) %>% 
  dplyr::select(gff_symbol, gff_id, dmel_symbol)
ensembl <- useEnsembl(biomart = "ensembl", 
                      dataset = "dmelanogaster_gene_ensembl", 
                      mirror = "www")
gene_list <- gsub("^gene-", "", rownames(selected_genes))
gene_list_tbl <- tibble(Name = gene_list, order = seq_along(gene_list))
gene_list_dmel <- gene_list_tbl %>%
  left_join(anno_tab, by = "Name") %>%
  dplyr::mutate(gene.name = ifelse(is.na(gene_symbol), Name, gene_symbol)) %>%
  dplyr::select(-gene_symbol, -order) %>%
  as_tibble()
gene_list_dmel <- gene_list_dmel %>%
  dplyr::left_join(gff2dmel, by = c("gene.name" = "gff_symbol")) %>%
  dplyr::distinct(gene.name, .keep_all = TRUE)
sig.tab.anno <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","entrezgene_id"),
                         filters = "ensembl_gene_id",
                         values = gene_list_dmel$dmel_symbol,
                         mart = ensembl)

gene_list_dmel <- left_join(gene_list_dmel, sig.tab.anno, by = join_by("dmel_symbol" == "ensembl_gene_id")) %>%
  dplyr::distinct(gene.name, .keep_all = TRUE) %>%
  dplyr::mutate(Name = ifelse(is.na(external_gene_name), Name, external_gene_name)) %>%
  dplyr::mutate(gff_id = ifelse(is.na(gff_id), Name, gff_id))

rownames(selected_genes) <- gene_list_dmel$Name

tiff('./Articoli/G.fuscipes_Spiroplasma_RNAseq/Gff-Spiro_2025-01-31/Figures/Fig_2.tiff', width=10.5,units="in", height=9, res=300)
pheatmap::pheatmap(selected_genes,
                   cluster_rows=F, 
                   show_rownames=T,
                   show_colnames=F,
                   cluster_cols=T, 
                   annotation_col=df,
                   annotation_colors = list(Tissue =c(MAGs="#24872195", Testes="#daa52095")),
                   colorRampPalette(c("steelblue","grey95","salmon"))(100),
                   fontsize = 12, cellheight=11,cellwidth=40)
dev.off()

## expression boxplots
expr_df <- assay(ntd)

resLFC_neg_vio <- bind_rows(testes.biased %>% dplyr::mutate(Tissue = "Testes"),
                            mags.biased %>% dplyr::mutate(Tissue = "MAGs")) %>%
  dplyr::mutate(logbase = log(baseMean))
resLFC_neg_vio$Tissue <- factor(resLFC_neg_vio$Tissue, levels = c("MAGs", "Testes"))
summary(resLFC_neg_vio$logbase ~ resLFC_neg_vio$Tissue)
resLFC_neg_vio %>%
  group_by(Tissue) %>%
  summarise(SD = sd(logbase, na.rm = TRUE))

png('boxplot_neg-testis-glands.png', width=13,units="in", height=10, res=300)
boxplot_neg <- ggplot(resLFC_neg_vio, aes(x=Tissue, y=logbase)) + 
  geom_boxplot(aes(color=Tissue, fill = Tissue), outliers = F) +
  geom_jitter(aes(color=Tissue), alpha = 0.4) +
  scale_fill_manual(values = c("#7bb779", "#e8c979")) +
  scale_color_manual(values = c("#248723", "#daa520")) +
  stat_compare_means(method = "t.test", label =  "p.signif", label.x = 1.5, label.y = 12.3, size = 8) +
  guides(colour = guide_legend(override.aes = list(size=16))) + 
  stat_summary(fun.y=mean, geom="point", size=2, color="black") +
  theme_classic() +
  theme(legend.position = 'right',
        axis.title.x = element_blank(),
        text = element_text(size = 16)) + 
  grids(linetype = "dashed") 
print(boxplot_neg)
dev.off()

tiff('./Articoli/G.fuscipes_Spiroplasma_RNAseq/Gff-Spiro_2025-01-31/Figures2/Fig_3.tiff', width=13,units="in", height=7, res=300)
plot_grid(volcano,boxplot_neg,nrow=1, ncol=2,labels=c("a", "b"),label_size = 15,rel_widths = c(1.5,1))
dev.off()


summary(resLFC_neg_vio$logbase ~ resLFC_neg_vio$Tissue)

#### PCA

vsd <- vst(dds_neg, blind=FALSE)
tiff('./Articoli/G.fuscipes_Spiroplasma_RNAseq/Gff-Spiro_2025-01-31/Figures2/Fig_1.tiff', width=10,units="in", height=8, res=300)
pcaData <- DESeq2::plotPCA(vsd, intgroup="tissue", returnData=TRUE)
pcaData$tissue <- as.factor(ifelse(pcaData$tissue  == "glands", "MAGs", "Testes"))
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2)) +
  geom_point(size=10,aes(color = tissue),  alpha = 0.5) +
  theme_classic() +
  labs(title = "") +
  theme(text=element_text(size = 30))+
  scale_color_manual(values = c("forestgreen","goldenrod"), name = "Tissue")  + 
  grids(linetype = "dashed") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(legend.position = 'right',
        text = element_text(size = 16))
dev.off()

#### most expresed genes
gene.counts <- assay(vsd) %>% as.data.frame()
gene.counts$gene <- gsub("^gene-", "", rownames(gene.counts)) 
gene.counts <- gene.counts %>% dplyr::select(gene, everything())
write_tsv(gene.counts, file = "glossina_neg-counts.tsv")

# Analisi con annotazioni de geni

anno_df <- read_tsv("./degs_tables/glossina_annotation_table.tsv")
sig_genes <- resLFC_neg_export%>%
  filter(padj < 0.05 & abs(log2FoldChange) > 2)
sig_genes$gene <- gsub("^gene-", "", sig_genes$gene)

# filtraggio per nome del gene - esempio: filter(anno_df, gene=="LOC119643313")

longest_substr <- function(x){substr(x[1], start = 1, stop = lcprefix(x[1], x[length(x)]))}
anno_sig_genes <- sig_genes %>% 
  mutate(product = map_chr(gene, ~ {
    df <- filter(anno_df, gene==.x, gbkey == "mRNA")
    if (nrow(df) == 0) {
      NA_character_
    } else {
      longest_substr(df$product)
    }
  }))

df <- anno_sig_genes
write_tsv(df, file="./glossina-negatives_annotated-significant-genes.tsv")
write_xlsx(df, path="./glossina-negatives_annotated-significant-genes.xlsx")

up.df <- df[df$log2FoldChange > 2,]
down.df <- df[df$log2FoldChange < -2,]
write_tsv(up.df, file="./glossina-negatives_annotated-upregulated-testes-genes.tsv")
write_xlsx(up.df, path="./glossina-negatives_annotated-upregulated-testes-genes.xlsx")
write_tsv(down.df, file="./glossina-negatives_annotated-downregulated-testes-genes.tsv")
write_xlsx(down.df, path="./glossina-negatives_annotated-downregulated-testes-genes.xlsx")



# check dell'espressione nei tessuti dei DEG significativi tra testis positivi e negativi 

# GO annotation of most significant expressed genes

gaf_file <- getGAF("./GCF_014805625.2_Yale_Gfus_2_genomic.gaf")
go_df <- as.data.frame(cbind(gaf_file@relations$elements,gaf_file@relations$sets))
colnames(go_df) <- c("gene","go_id")
go_df <- as_tibble(go_df)
DEGs <- sig_genes$gene
go_DEGs <- go_df[go_df$gene %in% DEGs,]
go_info <- as_tibble(get_names(go_DEGs$go_id))
go_DEGs_df <- unique(left_join(go_DEGs,go_info, by='go_id'))
# fare file con i termini di gene ontology

# enrichment analysis dei DEGs
### preparing the ranked list of genes using the shrinked lfc results for clusterProfiler
gsea_lfc <- sort(na.omit(resLFC_neg[resLFC_neg$padj<0.05,]$log2FoldChange), decreasing = TRUE)
names(gsea_lfc) <- gsub("^gene-", "", names(gsea_lfc))
upregulated <- names(gsea_lfc[gsea_lfc > 2])
downregulated <- names(gsea_lfc[gsea_lfc < -2])

go_df_new <- go_df %>%
  group_by(gene) %>%
  summarize(go_id = paste(go_id, collapse = ",")) %>%
  as.data.frame()
rownames(go_df_new) <- go_df_new$gene
go_df_new <- go_df_new[2]
# write.table(go_df_new, "gene2GO.txt")

get_gene_list <- function(genes, go_df){
  geneNames <- rownames(go_df)
  geneList <- factor(as.integer(geneNames %in% genes))
  names(geneList) <- geneNames
  geneList
}

up_list <- get_gene_list(upregulated, go_df_new)
down_list <- get_gene_list(downregulated, go_df_new)

gsea_topGO <- function(geneList, geneID2GO, kind){
  
  GOdata <- new("topGOdata", ontology = kind, allGenes = geneList,
                annot = annFUN.gene2GO, gene2GO = geneID2GO)
  
  resultTopGO.elim <- runTest(GOdata, algorithm = "elim", statistic = "Fisher" )
  Res <- GenTable(GOdata, elimKS = resultTopGO.elim,
                  orderBy = "elimKS", 
                  topNodes = 356)
  Res <- Res[Res$elimKS<0.05,]
  Res <- Res[,c("GO.ID","Term","elimKS")]
  Res$Term <- gsub(" [a-z]*\\.\\.\\.$", "", Res$Term)
  Res$Term <- gsub("\\.\\.\\.$", "", Res$Term)
  Res$Term <- paste(Res$GO.ID, Res$Term, sep=", ")
  Res$Term <- factor(Res$Term, levels=rev(Res$Term))
  Res
}

complete_gsea <- function(geneList, geneID2GO){
  bp_res <- gsea_topGO(geneList, geneID2GO, "BP") 
  cc_res <- gsea_topGO(geneList, geneID2GO, "CC")
  mf_res <- gsea_topGO(geneList, geneID2GO, "MF")
  bp_res$domain <- "BP"
  cc_res$domain <- "CC"
  mf_res$domain <- "MF"
  res <- rbind(bp_res, cc_res, mf_res)
  res
}

geneID2GO <- readMappings("gene2GO.txt", IDsep = ",", sep = " ")
up_res <- complete_gsea(up_list, geneID2GO)
down_res <- complete_gsea(down_list, geneID2GO)
up_res$regulation <- "upregulated"
down_res$regulation <- "downregulated"
final_res <- as_tibble(rbind(up_res,down_res))
final_res$enrichment <- ifelse(final_res$regulation == "upregulated", -log10(as.numeric(final_res$elimKS)), log10(as.numeric(final_res$elimKS)))
final_res <- final_res[order(final_res$enrichment, decreasing = T),]
BP <- final_res[final_res$domain == "BP",] %>%
  group_by(GO.ID) %>%
  filter(n() == 1) %>%
  ungroup()
CC <- final_res[final_res$domain == "CC",] %>%
  group_by(GO.ID) %>%
  filter(n() == 1) %>%
  ungroup()
MF <- final_res[final_res$domain == "MF",] %>%
  group_by(GO.ID) %>%
  filter(n() == 1) %>%
  ungroup()

make_gsea_plot <- function(gsea_res, domain){
  # gsea_res has to have Term (GO term), enrichment (log(pval) of the topGO) and regulation (upregulated or downregulated)
  gsea_res <- na.omit(gsea_res)
  p <- ggplot(gsea_res, aes(x=reorder(Term, enrichment), y=enrichment, fill=regulation)) +
    stat_summary(geom = "bar", fun = mean, position = "dodge") +
    xlab(domain) +
    ylab("Enrichment") +
    # ggtitle(substitute(paste(italic(Spiroplasma), domain," influenced by DEGs in testis versus MAGs"))) +
    scale_y_continuous(breaks = round(seq(min(gsea_res$enrichment), max(gsea_res$enrichment), by = 2), 1)) +
    theme_classic(base_size=26) +
    theme(
      legend.position='right',
      legend.background=element_rect(),
      plot.title=element_text(angle=0, size=28, vjust=1, hjust=0.5, face="bold"),
      axis.text.x=element_text(angle=0, size=20, hjust=1.10),
      axis.text.y=element_text(angle=0, size=18,  vjust=0.5),
      axis.title=element_text(size=24),
      legend.key=element_blank(),     #removes the border
      legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
      legend.text=element_text(size=20),  #Text size
      title=element_text(size=20)) +
    guides(fill=guide_legend(title="Tissue",override.aes=list(size=2.5))) +
    coord_flip()+
    scale_fill_manual(values=alpha(c("forestgreen","goldenrod"),.3), labels=c("MAGs", "Testes")) + 
    grids(linetype = "dashed") 
  print(p)
}

png('./Articoli/G.fuscipes_Spiroplasma_RNAseq/Gff-Spiro_2025-11-19/Supplementary_Material/Fig_S1.png', width=15,units="in", height=15, res=300)
mfp <- make_gsea_plot(MF, "Molecular Functions")
mfp
dev.off()

png('./Articoli/G.fuscipes_Spiroplasma_RNAseq/Gff-Spiro_2025-11-19/Supplementary_Material/Fig_S2.png', width=15,units="in", height=15, res=300)
bpp <- make_gsea_plot(BP, "Biological Processes")
bpp
dev.off()

tiff('glossina_negatives_gsea_CC.tiff', width=20,units="in", height=25, res=300)
ccp <- make_gsea_plot(CC, "Cellular Components")
ccp
dev.off()



tiff('glossina_negatives_gsea.tiff', width=18,units="in", height=15, res=300)
ggplot(final_res, aes(x=reorder(Term, enrichment), y=enrichment, fill=domain)) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("GO terms") +
  ylab("Enrichment") +
  ggtitle(substitute(paste(italic(Glossina)," domains influenced by DEGs in testis versus MAGs"))) +
  scale_y_continuous(breaks = round(seq(min(final_res$enrichment), max(final_res$enrichment), by = 2), 1)) +
  theme_bw(base_size=26) +
  theme(
    legend.position='right',
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=28, vjust=1, hjust=0.5, face="bold"),
    axis.text.x=element_text(angle=0, size=20, hjust=1.10),
    axis.text.y=element_text(angle=0, size=20,  vjust=0.5),
    axis.title=element_text(size=24),
    legend.key=element_blank(),
    legend.key.size=unit(1, "cm"),
    legend.text=element_text(size=20),
    title=element_text(size=20)) +
  guides(fill=guide_legend(title="Regulation",override.aes=list(size=2.5))) +
  coord_flip() +
  scale_fill_manual(values=alpha(brewer.pal(n = 3, name = "Set2"),.7))
dev.off()

write_tsv(final_res%>%
            group_by(GO.ID) %>%
            filter(n() == 1) %>%
            ungroup(), file="./glossina-negatives_GSEA.tsv")
write_xlsx(final_res%>%
             group_by(GO.ID) %>%
             filter(n() == 1) %>%
             ungroup(), path="./glossina-negatives_GSEA.xlsx")

# uploading table for droso translation for table S1
orthologs <- read.table('gff2Dmel_gene_names.txt', header=T)
res_neg_export$gene <- gsub("^gene-", "", res_neg_export$gene) 
write_xlsx(res_neg_export %>% left_join(orthologs, by=c("gene"="gff_id")), 
           path='table_s1.xlsx')

# figure 1 - update 20/11/25
p1 <- ggplot(pcaData, aes(PC1, PC2)) +
  geom_point(size=3,aes(color = tissue),  alpha = 0.5) +
  theme_classic() +
  labs(title = "") +
  theme(text=element_text(size = 30))+
  scale_color_manual(values = c("forestgreen","goldenrod"), name = "Tissue")  + 
  grids(linetype = "dashed") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(legend.position = 'right',
        text = element_text(size = 16))

p2 <- pheatmap::pheatmap(selected_genes,
                         cluster_rows=F, 
                         show_rownames=T,
                         show_colnames=F,
                         cluster_cols=T, 
                         annotation_col=df,
                         annotation_colors = list(Tissue =c(MAGs="#24872195", Testes="#daa52095")),
                         colorRampPalette(c("steelblue","grey95","salmon"))(100),
                         fontsize = 12, silent=T)

p3 <- volcano 

p4 <- boxplot_neg

blank <- ggplot() + theme_void()

left_col <- plot_grid(
  p1, p4,
  ncol = 1,
  labels = c("a", "d"),
  label_size = 12,
  rel_heights = c(1, 1)
)

middle_col <- plot_grid(
  ggplotify::as.ggplot(p2$gtable),
  labels = "b",
  label_size = 12
)

right_col <- plot_grid(
  blank, 
  p3, 
  blank, 
  ncol = 1,
  labels = c("","c",""),
  rel_heights = c(0.2, 0.6, 0.2),  
  label_size = 12
)

fig1 <- plot_grid(
  left_col, middle_col, right_col,
  ncol = 3,
  rel_widths = c(1, 1.3, 1.3)
)

png("Figure_1.png", 
    units = "in", width=15, height = 8, res=300)
print(fig1)
dev.off()

save.image('negative_analysis.RData')

#########################################



############ ANALISI LESCAI ############


dds_stratified <- DESeqDataSetFromMatrix(countData=assay(dds), colData=sampledata, design = ~ tissue + spiro)
dds_multi <- DESeqDataSetFromMatrix(countData=assay(dds), colData=sampledata, design = ~ condition)

######## filtering on number of samples

keep_stratified <- rowSums(counts(dds_stratified) >= 10) >= 6
dds_stratified <- dds_stratified[keep_stratified,]

keep_multi <- rowSums(counts(dds_multi) >= 10) >= 6
dds_multi <- dds_multi[keep_multi,]

####### DE analysis for both

dds_stratified <- DESeq(dds_stratified)

dds_multi <- DESeq(dds_multi)

################# results

res_stratified = results(dds_stratified)
res_stratified = res_stratified[order(res_stratified$padj),]
res_stratified_export = tibble(
  bind_cols(
    gene = row.names(res_stratified),
    as_tibble(res_stratified)
  )
)
  
write_tsv(res_stratified_export, file = "glossina_stratified-design_de-results.tsv")

res_multi = results(dds_multi)
res_multi = res_multi[order(res_multi$padj),]

resultsNames(dds_stratified)

resultsNames(dds_multi)

#### shrink for visualisation

resLFC_stratified <- lfcShrink(dds_stratified, coef="spiro_positive_vs_negative", type="apeglm")

resLFC_multi <- lfcShrink(dds_multi, coef="condition_TESTISPOS_vs_GLANDSNEG", type="apeglm")

resLFC_neg <- lfcShrink(dds_stratified, coef="tissue_testis_vs_glands", type="apeglm")

###### checks

plotMA(resLFC_neg)

## the top gene is the same in the 2 designs
plotCounts(dds_stratified, gene=which.min(res_stratified$padj), intgroup=c("tissue","spiro"))
plotCounts(dds_multi, gene=which.min(res_multi$padj), intgroup=c("tissue","spiro"))

library("pheatmap")
ntd <- normTransform(dds_stratified)
select <- order(rowMeans(counts(dds_stratified,normalized=TRUE)),
                decreasing=TRUE)[1:40]
df <- as.data.frame(colData(dds_stratified)[,c("spiro","tissue")])
out <- pheatmap(assay(ntd)[select,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df)
interesting.genes <- rownames(assay(ntd)[out$tree_row[["order"]],])[2:21] # estraggo i geni che sembrano differenti all'interno dei testicoli


#### PCA

vsd <- vst(dds_stratified, blind=FALSE)
plotPCA(vsd, intgroup=c("tissue", "spiro"))

##########################################
## extract only testis from multi design
########################################

results_testis_pos_neg = results(dds_multi, contrast = c("condition", "TESTISPOS", "TESTISNEG"))
results_testis_pos_neg = results_testis_pos_neg[order(results_testis_pos_neg$padj),]
results_testis_pos_neg_export = tibble(
  bind_cols(
    gene = row.names(results_testis_pos_neg),
    as_tibble(results_testis_pos_neg)
  )
)
write_tsv(results_testis_pos_neg_export, file = "glossina_results_contrast_testis-only.tsv")

vsd_testisonly = vsd[, vsd$condition %in% c("TESTISPOS", "TESTISNEG")]
plotPCA(vsd_testisonly, intgroup=c("spiro"))

plotCounts(dds_multi, "gene-LOC119644472")

save.image("glossina_analysis_20230306.RData")
# load("glossina_analysis_20230306.RData")