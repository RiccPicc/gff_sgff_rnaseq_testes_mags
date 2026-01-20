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
library(biomaRt)
library(ggrepel)

quartzFonts(helvetica = c("Helvetica Neue Light", 
                          "Helvetica Neue Bold", "Helvetica Neue Light Italic", 
                          "Helvetica Neue Bold Italic"))
par(family="helvetica")

load("testes_mag_analysis_data.RData") 

# View(sampledata)
# la tabella mostra tutti i campioni a disposizione
# seleziono solo i campioni testicoli per fare lanciare l'analisi differenziale 
# solo su questi campioni
glandsdata <- sampledata %>% 
  filter(tissue=="glands")

idx <- c(13:24) # campioni da eliminare (i MAG sono i primi 12 campioni)
dds_glands <- dds[,-idx] # rimuovo i campioni dei MAG
dds_glands <- DESeqDataSetFromMatrix(countData=assay(dds_glands), colData=glandsdata, design = ~ spiro)


######## filtering on number of samples

keep_glands <- rowSums(counts(dds_glands) >= 10) >= 6
dds_glands <- dds_glands[keep_glands,]

#### DE analysis for both

dds_glands <- DESeq(dds_glands)

################# results

res_glands = results(dds_glands)
res_glands = res_glands[order(abs(res_glands$log2FoldChange),decreasing=T),]
res_glands_export = tibble(
  bind_cols(
    gene = row.names(res_glands),
    as_tibble(res_glands)
  )
)

#### shrink for visualisation

resLFC_glands <- lfcShrink(dds_glands, coef="spiro_positive_vs_negative", type="apeglm")
resLFC_glands <- na.omit(resLFC_glands[order(abs(resLFC_glands$log2FoldChange), decreasing=T),]) 
resLFC_glands_export <- tibble(
  bind_cols(
    gene = row.names(resLFC_glands),
    as_tibble(resLFC_glands)
  )
)

write_tsv(resLFC_glands_export, file = "glossina_glands-design_de-results.tsv")
write_xlsx(resLFC_glands_export, path = "glossina_glands-design_de-results.xlsx")

###### checks

plotMA(resLFC_glands)

##### adding Dmel gene names
gff2dmel <- read.csv("./PostDoc/Spiroplasma/Gff2Dmel_proteins/Gff2Dmel.txt", sep="\t",header=F)
gff.anno <- read.csv("./PostDoc/Spiroplasma/degs_tables/glossina_annotation-full_table.tsv", sep="\t")
gff.anno <- gff.anno %>% dplyr::filter(grepl("LOC",Name)) %>% dplyr::select(gene_symbol, Name)
gff2dmel <- merge(gff2dmel, gff.anno, by.x = "V1", by.y = "gene_symbol") %>% 
  dplyr::rename(gff_symbol = V1, gff_id = Name, dmel_symbol = V2) %>% 
  dplyr::select(gff_symbol, gff_id, dmel_symbol)
resLFC_glands_dmel <- resLFC_glands_export %>% dplyr::mutate(gene = gsub("^gene-", "", gene))
resLFC_glands_dmel <- merge(gff2dmel, resLFC_glands_dmel, by.x = "gff_id", by.y = "gene")
ensembl <- useEnsembl(biomart = "ensembl", 
                      dataset = "dmelanogaster_gene_ensembl",
                      mirror = "www")
sig.tab.updated <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","entrezgene_id"),
                         filters = "ensembl_gene_id",
                         values = resLFC_glands_dmel$dmel_symbol,
                         mart = ensembl)
resLFC_glands_dmel <- merge(resLFC_glands_dmel, sig.tab.updated, by.x = "dmel_symbol", by.y = "ensembl_gene_id") %>%
  dplyr::arrange(dplyr::desc(abs(log2FoldChange))) %>%
  dplyr::distinct(gff_symbol, .keep_all = TRUE) %>%
  dplyr::mutate(labs = ifelse(padj<0.05 & abs(log2FoldChange)>0.5, external_gene_name, NA)) %>%
  dplyr::mutate(l2FC = cut(log2FoldChange, 
                           breaks=c(-Inf,-1.5, -0.5, 0.5, 1.5,Inf), 
                           right = FALSE, 
                           labels = c("l2FC < -1.5", 
                                      "-1.5 > l2FC > -0.5", 
                                      "-0.5 > l2FC > 0.5\nor p < 0.05",
                                      "0.5 > l2FC > 1.5", 
                                      "l2FC > 1.5"))) %>%
  dplyr::mutate(l2FC = as.factor(ifelse(padj>0.05,"-0.5 > l2FC > 0.5\nor p < 0.05", as.character(l2FC))))


resLFC_glands_dmel$l2FC <- factor(resLFC_glands_dmel$l2FC, levels=c("l2FC < -1.5", 
                                                                    "-1.5 > l2FC > -0.5", 
                                                                    "-0.5 > l2FC > 0.5\nor p < 0.05",
                                                                    "0.5 > l2FC > 1.5", 
                                                                    "l2FC > 1.5"))
tiff('./Articoli/G.fuscipes_Spiroplasma_RNAseq/Gff-Spiro_2025-01-31/Figures2/Fig_6.tiff', width=10,units="in", height=8, res=300)
volcano_mags <- ggplot(resLFC_glands_dmel, aes(x=log2FoldChange, y=-log10(padj))) + 
  geom_point(aes(color = l2FC, size=abs(log2FoldChange)),  alpha = 0.3) + 
  geom_label_repel(
    aes(x=log2FoldChange, y=-log10(padj), label = labs),
    fontface = 'bold', 
    color = 'grey20',
    box.padding = 0.35, point.padding = 0.5,
    segment.color = 'grey50',
    max.overlaps=nrow(iris)
  ) +
  geom_hline(yintercept=-log10(0.05), linetype = "dashed", alpha = 0.8) +
  geom_vline(xintercept=c(0.5, -0.5),linetype = "dotted", alpha = 0.8) +
  geom_vline(xintercept=c(1.5, -1.5),linetype = "dashed", alpha = 0.8) +
  geom_hline(yintercept=15, linetype = "dashed", alpha = 0.8, color = NA) + # to harmonize with testes plot
  scale_color_manual(values = c("grey10","darksalmon"))  + 
  labs(color='l2FC', size='|l2FC|')+
  theme_classic() +
  theme(legend.position = 'right',
        text = element_text(size = 16)) + 
  grids(linetype = "dashed") + 
  guides(colour = guide_legend(override.aes = list(size=10))) +
  ylab(bquote("-" ~ log[10] ~ italic(P)))+
  xlab(bquote(log[2] ~ "fold change (l2FC)"))
print(volcano_mags)
dev.off()

tiff('H:/Il mio Drive/Articoli/G.fuscipes_Spiroplasma_RNAseq/Figures/Fig_3.tiff', width=20,units="in", height=10, res=300)
plot_grid(volcano_mags+theme(legend.position = 'none'), 
          volcano_testes+theme(legend.position = 'left'),
          nrow=1, 
          ncol=2,
          labels=c("a", "b"),
          label_size = 12,
          rel_widths = c(1.25,1.5))
dev.off()

## the top gene is the same in the 2 designs
plotCounts(dds_glands, gene=which.min(res_glands$padj), intgroup="spiro")


# analizzo i geni più espressi nei testicoli overall
expr_cutoff <- 40
ntd <- normTransform(dds_glands)
select <- order(rowMeans(counts(dds_glands,normalized=TRUE)),
                decreasing=TRUE)[1:expr_cutoff]
df <- as.data.frame(colData(dds_glands)[,c("spiro")])
rownames(df) <- colnames(assay(ntd))
out <- pheatmap(assay(ntd)[select,], cluster_rows=TRUE, show_rownames=TRUE,
                cluster_cols=TRUE, annotation_col=df)
df <- as.data.frame(assay(ntd)[select,])
df$counts_mean <- rowMeans(df)
most_expr_genes <- as_tibble(cbind(rownames(df),df$counts_mean)) %>% dplyr::rename(gene=V1, counts_mean=V2)
most_expr_genes$gene <- gsub("^gene-", "", most_expr_genes$gene)

anno2 <- read_tsv("glossina_annotation2_table.tsv")
anno2$ID <- gsub("^exon_", "", anno2$ID)
anno2$ID <- gsub("_t\\d+-E\\d+$", "", anno2$ID)
anno1 <- read_tsv("glossina_annotation_table.tsv")
pseudo.df <- read_tsv("glossina_pseudo-annotation_table.tsv")
anno2$ID <- gsub("\\..*$", "", anno2$ID)
anno2$symbol <- ifelse(grepl("^LOC", anno2$ID), anno2$ID, NA)
anno2$group_ID <- ifelse(grepl("^LOC", anno2$ID), anno2$gene_id, anno2$ID)
anno2 <- anno2 %>%
  group_by(group_ID) %>%
  summarize(gene_id = dplyr::first(na.omit(gene_id)),
            gene_symbol = dplyr::first(na.omit(symbol)),
            description = dplyr::first(na.omit(description))) %>%
  dplyr::select(-group_ID)
anno2 <- dplyr::rename(anno2, gene = gene_symbol)
pseudo.df$ebi_biotype <- paste("nbis",pseudo.df$ebi_biotype,row_number(pseudo.df),sep="-")
pseudo.df$description <- "unspecified psudogene product"
pseudo.df <- dplyr::rename(pseudo.df, gene=ebi_biotype)
pseudo.df <- dplyr::rename(pseudo.df, gene_id=ID)

add_pseudo <- function(df, pseudo.df = NULL){
  merged_df <- base::merge(df, pseudo.df, by = "gene", all.x = TRUE)
  merged_df$gene_id.x[is.na(merged_df$gene_id.x)] <- merged_df$gene_id.y[is.na(merged_df$gene_id.x)]
  merged_df$description.x[is.na(merged_df$description.x)] <- merged_df$description.y[is.na(merged_df$description.x)]
  merged_df <- merged_df[, !(names(merged_df) %in% c("gene_id.y", "description.y"))]
  colnames(merged_df) <- colnames(df)
  df <- merged_df
  df
}
most_expr_genes <- unique(dplyr::left_join(most_expr_genes,anno2,by="gene"))
most_expr_genes <- add_pseudo(most_expr_genes, pseudo.df)
write_tsv(most_expr_genes, file = "glossina_glands-design_most-expr_complete.tsv")
write_xlsx(most_expr_genes, path = "glossina_glands-design_most-expr-complete.xlsx")


select <- order(rowMeans(counts(dds_glands,normalized=TRUE)[,1:6]),
                decreasing=TRUE)[1:expr_cutoff]
df <- as.data.frame(colData(dds_glands[1:6])[,c("spiro")])
rownames(df) <- colnames(assay(ntd))
out <- pheatmap(assay(ntd)[select,], cluster_rows=TRUE, show_rownames=TRUE,
                cluster_cols=TRUE, annotation_col=df)
df <- as.data.frame(assay(ntd)[select,])
df$counts_mean <- rowMeans(df)
most_expr_genes <- as_tibble(cbind(rownames(df),df$counts_mean)) %>% dplyr::rename(gene=V1, counts_mean=V2)
most_expr_genes$gene <- gsub("^gene-", "", most_expr_genes$gene)
most_expr_genes <- unique(dplyr::left_join(most_expr_genes,anno2,by="gene"))
most_expr_genes <- add_pseudo(most_expr_genes, pseudo.df)
write_tsv(most_expr_genes, file = "glossina_glands-design_most-expr-neg_complete.tsv")
write_xlsx(most_expr_genes, path = "glossina_glands-design_most-expr-neg-complete.xlsx")

select <- order(rowMeans(counts(dds_glands,normalized=TRUE)[,6:12]),
                decreasing=TRUE)[1:expr_cutoff]
df <- as.data.frame(colData(dds_glands[1:6])[,c("spiro")])
rownames(df) <- colnames(assay(ntd))
out <- pheatmap(assay(ntd)[select,], cluster_rows=TRUE, show_rownames=TRUE,
                cluster_cols=TRUE, annotation_col=df)
df <- as.data.frame(assay(ntd)[select,])
df$counts_mean <- rowMeans(df)
most_expr_genes <- as_tibble(cbind(rownames(df),df$counts_mean)) %>% dplyr::rename(gene=V1, counts_mean=V2)
most_expr_genes$gene <- gsub("^gene-", "", most_expr_genes$gene)
most_expr_genes <- unique(dplyr::left_join(most_expr_genes,anno2,by="gene"))
most_expr_genes <- add_pseudo(most_expr_genes, pseudo.df)
write_tsv(most_expr_genes, file = "glossina_glands-design_most-expr-pos_complete.tsv")
write_xlsx(most_expr_genes, path = "glossina_glands-design_most-expr-pos-complete.xlsx")


ntd <- normTransform(dds_glands)
df <- as.data.frame(colData(dds_glands)[,"spiro"])
rownames(df) <- colnames(ntd)
select <- rownames(as.data.frame(resLFC_glands) %>%
                     filter(padj < 1 & abs(log2FoldChange) > 0.1))[1:40]
p <-pheatmap(assay(ntd)[select,], cluster_rows=TRUE, show_rownames=TRUE,
             cluster_cols=TRUE, annotation_col=df)
p
dev.off()
png('glossina_glands_heatmap-40genes.png', width=13,units="in", height=10, res=300)
pheatmap(assay(ntd)[select,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df)
dev.off()

#### PCA

vsd <- vst(dds_glands, blind=FALSE)
png('glossina_glands_PCA.png', width=13,units="in", height=10, res=300)
plotPCA(vsd, intgroup="spiro")
dev.off()



# Analisi con annotazioni de geni

anno_df <- read_tsv("./glossina_annotation_table.tsv")
sig_genes <- resLFC_glands_export%>%
  filter(padj < 0.05 & abs(log2FoldChange) > 0.5)
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
# usa la seguente linea per caricare direttamente il dataframe dei DEG dei testicoli posvs neg
# df <- read_tsv("./glossina-glands_annotated-significant-genes.tsv")
df$spiro <- ifelse(df$log2FoldChange < 0, "negative", "positive")

# check dell'espressione nei tessuti dei DEG significativi tra glands positivi e negativi 

neg_pos_DEGs <- read_tsv("./glossina_results_contrast_negatives-only.tsv")
# log2fc quando NEGATIVO è a favore dei glands, quando POSITIVO è a favore delle MAG
neg_pos_DEGs$gene <- gsub("^gene-", "", neg_pos_DEGs$gene)
neg_pos_DEGs$tissue <- ifelse(neg_pos_DEGs$log2FoldChange < 0, "glands", "testis")
neg_pos_DEGs$tissue_level <- ifelse(abs(neg_pos_DEGs$log2FoldChange) < 0.5, "same", "different")
neg_pos_DEGs <- neg_pos_DEGs %>%
  rename(log2FC_tissue = log2FoldChange) %>%
  rename(lfcSE_tissue = lfcSE) %>%
  rename(padj_tissue = padj) %>%
  dplyr::select(gene, tissue, tissue_level, log2FC_tissue, lfcSE_tissue, padj_tissue)

df <- left_join(df, neg_pos_DEGs, by = join_by(gene == gene))

write_tsv(df, file="./glossina-glands_annotated-significant-genes.tsv")
write_xlsx(df, path="./glossina-glands_annotated-significant-genes.xlsx")

# GO annotation of most significant expressed genes

gaf_file <- getGAF("./GCF_014805625.2_Yale_Gfus_2_genomic.gaf")
go_df <- as.data.frame(cbind(gaf_file@relations$elements,gaf_file@relations$sets))
colnames(go_df) <- c("gene","go_id")
go_df <- as_tibble(go_df)
DEGs <- sig_genes$gene
go_DEGs <- go_df[go_df$gene %in% DEGs,]
go_info <- as_tibble(get_names(go_DEGs$go_id)) # l'unico gene DEG non ha termini GO connessi provenienti dal .gaf file
go_DEGs_df <- unique(left_join(go_DEGs,go_info, by='go_id'))
# fare file con i termini di gene ontology

# enrichment analysis dei DEGs
### preparing the ranked list of genes using the shrinked lfc results for clusterProfiler
gsea_lfc <- sort(na.omit(resLFC_glands[resLFC_glands$padj<0.05,]$log2FoldChange), decreasing = TRUE)
names(gsea_lfc) <- gsub("^gene-", "", names(gsea_lfc))
upregulated <- names(gsea_lfc[gsea_lfc > 0.5])
downregulated <- names(gsea_lfc[gsea_lfc < -0.5])

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
BP <- final_res[final_res$domain == "BP",]
CC <- final_res[final_res$domain == "CC",]
MF <- final_res[final_res$domain == "MF",]

write_tsv(final_res, file="./glossina-glands_GSEA.tsv")
write_xlsx(final_res, path="./glossina-glands_GSEA.xlsx")

make_gsea_plot <- function(gsea_res, domain){
  # gsea_res has to have Term (GO term), enrichment (log(pval) of the topGO) and regulation (upregulated or downregulated)
  p <- ggplot(gsea_res, aes(x=reorder(Term, enrichment), y=enrichment, fill=regulation)) +
    stat_summary(geom = "bar", fun = mean, position = "dodge") +
    xlab(domain) +
    ylab("Enrichment") +
    ggtitle(substitute(paste(italic(Spiroplasma)," affected ", domain," in ",italic(Glossina)," glands"))) +
    scale_y_continuous(breaks = round(seq(min(gsea_res$enrichment), max(gsea_res$enrichment), by = 2), 1)) +
    theme_bw(base_size=26) +
    theme(
      legend.position='right',
      legend.background=element_rect(),
      plot.title=element_text(angle=0, size=28, vjust=1, hjust=0.5, face="bold"),
      axis.text.x=element_text(angle=0, size=20, hjust=1.10),
      axis.text.y=element_text(angle=0, size=20,  vjust=0.5),
      axis.title=element_text(size=24),
      legend.key=element_blank(),     #removes the border
      legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
      legend.text=element_text(size=20),  #Text size
      title=element_text(size=20)) +
    guides(fill=guide_legend(title="Regulation",override.aes=list(size=2.5))) +
    coord_flip()+
    scale_fill_manual(values=alpha(c("royalblue3","red3"),.3))
  print(p)
}

png('glossina_glands_gsea_MF.png', width=18,units="in", height=15, res=300)
make_gsea_plot(MF, "Molecular Functions")
dev.off()

png('glossina_glands_gsea_BP.png', width=18,units="in", height=15, res=300)
make_gsea_plot(BP, "Biological Processes")
dev.off()

png('glossina_glands_gsea_CC.png', width=18,units="in", height=15, res=300)
make_gsea_plot(CC, "Cellular Components")
dev.off()

png('glossina_glands_gsea.png', width=18,units="in", height=15, res=300)
ggplot(final_res, aes(x=reorder(Term, enrichment), y=enrichment, fill=domain)) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("GO terms") +
  ylab("Enrichment") +
  ggtitle(substitute(paste(italic(Spiroplasma)," affected domains in ",italic(Glossina)," glands"))) +
  scale_y_continuous(breaks = round(seq(min(final_res$enrichment), max(final_res$enrichment), by = 2), 1)) +
  theme_bw(base_size=26) +
  theme(
    legend.position='right',
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=28, vjust=1, hjust=0.5, face="bold"),
    axis.text.x=element_text(angle=0, size=20, hjust=1.10),
    axis.text.y=element_text(angle=0, size=20,  vjust=0.5),
    axis.title=element_text(size=24),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=20),  #Text size
    title=element_text(size=20)) +
  guides(fill=guide_legend(title="Regulation",override.aes=list(size=2.5))) +
  coord_flip()+
  scale_fill_manual(values=alpha(brewer.pal(n = 3, name = "Set2"),.7))
dev.off()

