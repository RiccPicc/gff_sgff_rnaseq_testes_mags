rm(list=ls(all=T))
library(DESeq2)
library(tidyverse)
library(writexl)
library(Biostrings)
library("pheatmap")
library(EnhancedVolcano)
library(BaseSet)
library(GOfuncR)
library(topGO)
library(enrichplot)
library("RColorBrewer")
library(ggplot2)

load("deseq2.dds.RData")

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
posdata <- sampledata %>% 
  filter(spiro=="positive")

idx <- c(1:6,13:18) # campioni da eliminare 
dds_pos <- dds[,-idx] # mantengo i campioni dei MAG
dds_pos <- DESeqDataSetFromMatrix(countData=assay(dds_pos), colData=posdata, design = ~ tissue)


######## filtering on number of samples

keep_pos <- rowSums(counts(dds_pos) >= 10) >= 6
dds_pos <- dds_pos[keep_pos,]

#### DE analysis for both

dds_pos <- DESeq(dds_pos)

################# results

res_pos = results(dds_pos)
res_pos = res_pos[order(res_pos$padj),]
res_pos_export = tibble(
  bind_cols(
    gene = row.names(res_pos),
    as_tibble(res_pos)
  )
)


write_tsv(res_pos_export, file = "glossina_pos-design_de-results.tsv")
library(writexl)
write_xlsx(res_pos_export, path = "glossina_pos-design_de-results.xlsx")

#### shrink for visualisation

resLFC_pos <- lfcShrink(dds_pos, coef="tissue_testis_vs_glands", type="apeglm")
resLFC_pos <- na.omit(resLFC_pos[order(abs(resLFC_pos$log2FoldChange), decreasing=T),]) 
resLFC_pos_export <- tibble(
  bind_cols(
    gene = row.names(resLFC_pos),
    as_tibble(resLFC_pos)
  )
)

###### checks

plotMA(resLFC_pos)

##### volcano plot
dev.off()
png('volcano-plot_pos-testis-glands.png', width=13,units="in", height=10, res=300)
EnhancedVolcano(resLFC_pos,
                lab = gsub("^gene-", "", rownames(resLFC_pos)),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlab = bquote(~Log[2]~ 'fold change'),
                title = bquote('Testis' ~ italic(versus) ~ 'MAGs'),
                subtitle = bquote('Negative to' ~ italic(Spiroplasma)),
                selectLab = c(),
                pCutoff = 0.05,
                FCcutoff = 2,
                ylim = c(0, -log10(10e-7)),
                xlim = c(-6, 6),
                pointSize = c(ifelse(abs(resLFC_pos$log2FoldChange)>2, 2.5, 1)),
                labSize = 3,
                caption = bquote(~Log[2]~ "fold change cutoff = 2; p-value cutoff = 0.05"),
                legendPosition = "None",
                colAlpha = 0.5,
                col = c("black", "grey30", "red4", "red3"),
                drawConnectors = TRUE,
                hline = c(10e-8),
                widthConnectors = 0.5)+
  ggplot2::coord_cartesian(xlim=c(-6, 6)) +
  ggplot2::scale_x_continuous(
    breaks=seq(-6,6, 1))
dev.off()
## the top gene is the same in the 2 designs
plotCounts(dds_pos, gene=which.min(res_pos$padj), intgroup="tissue")

ntd <- normTransform(dds_pos)
df <- as.data.frame(colData(dds_pos)[,"tissue"])
rownames(df) <- colnames(ntd)
select <- rownames(as.data.frame(resLFC_pos) %>%
                     filter(padj < 0.05 & abs(log2FoldChange) > 0.5))[1:40]
p <-pheatmap(assay(ntd)[select,], cluster_rows=TRUE, show_rownames=TRUE,
             cluster_cols=TRUE, annotation_col=df)
p
dev.off()
png('heatmap_pos_40genes.png', width=13,units="in", height=10, res=300)
pheatmap(assay(ntd)[select,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df)
dev.off()

#### PCA

vsd <- vst(dds_pos, blind=FALSE)
png('PCA_pos.png', width=13,units="in", height=10, res=300)
plotPCA(vsd, intgroup="spiro", ntop=80)
dev.off()


# Analisi con annotazioni de geni

anno_df <- read_tsv("./glossina_annotation_table.tsv")
sig_genes <- resLFC_pos_export%>%
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
  }))LOC119637269

df <- anno_sig_genes
write_tsv(df, file="./glossina-positives_annotated-significant-genes.tsv")
write_xlsx(df, path="./glossina-positives_annotated-significant-genes.xlsx")

up.df <- df[df$log2FoldChange > 2,]
down.df <- df[df$log2FoldChange < -2,]
write_tsv(up.df, file="./glossina-positives_annotated-upregulated-testes-genes.tsv")
write_xlsx(up.df, path="./glossina-positives_annotated-upregulated-testes-genes.xlsx")
write_tsv(down.df, file="./glossina-positives_annotated-downregulated-testes-genes.tsv")
write_xlsx(down.df, path="./glossina-positives_annotated-downregulated-testes-genes.xlsx")



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
gsea_lfc <- sort(na.omit(resLFC_pos[resLFC_pos$padj<0.05,]$log2FoldChange), decreasing = TRUE)
names(gsea_lfc) <- gsub("^gene-", "", names(gsea_lfc))
upregulated <- names(gsea_lfc[gsea_lfc > 2])
downregulated <- names(gsea_lfc[gsea_lfc < -2])

go_df_new <- go_df %>%
  group_by(gene) %>%
  dplyr::summarize(go_id = paste(go_id, collapse = ",")) %>%
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

make_gsea_plot <- function(gsea_res, domain){
  # gsea_res has to have Term (GO term), enrichment (log(pval) of the topGO) and regulation (upregulated or downregulated)
  gsea_res <- na.omit(gsea_res)
  p <- ggplot(gsea_res, aes(x=reorder(Term, enrichment), y=enrichment, fill=regulation)) +
    stat_summary(geom = "bar", fun = mean, position = "dodge") +
    xlab(domain) +
    ylab("Enrichment") +
    ggtitle(substitute(paste(italic(Spiroplasma), domain," influenced by DEGs in testis versus MAGs"))) +
    scale_y_continuous(breaks = round(seq(min(gsea_res$enrichment), max(gsea_res$enrichment), by = 2), 1)) +
    theme_bw(base_size=26) +
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
    guides(fill=guide_legend(title="Regulation",override.aes=list(size=2.5))) +
    coord_flip()+
    scale_fill_manual(values=alpha(c("royalblue3","red3"),.3))
  print(p)
}

png('glossina_posatives_gsea_MF.png', width=20,units="in", height=25, res=300)
make_gsea_plot(MF, "Molecular Functions")
dev.off()

png('glossina_posatives_gsea_BP.png', width=20,units="in", height=25, res=300)
make_gsea_plot(BP, "Biological Processes")
dev.off()

png('glossina_posatives_gsea_CC.png', width=20,units="in", height=25, res=300)
make_gsea_plot(CC, "Cellular Components")
dev.off()

png('glossina_posatives_gsea.png', width=18,units="in", height=15, res=300)
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
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=20),  #Text size
    title=element_text(size=20)) +
  guides(fill=guide_legend(title="Regulation",override.aes=list(size=2.5))) +
  coord_flip()+
  scale_fill_manual(values=alpha(brewer.pal(n = 3, name = "Set2"),.7))
dev.off()

write_tsv(final_res, file="./glossina-positives_GSEA.tsv")
write_xlsx(final_res, path="./glossina-positives_GSEA.xlsx")
