# script based on https://alexslemonade.github.io/refinebio-examples/03-rnaseq/pathway-analysis_rnaseq_02_gsea.html

# packages installation if required

bioc_install <- function(package){
  if (!(package %in% installed.packages())) {
    BiocManager::install(package, update = FALSE)
  }
}

bioc_install("clusterProfiler")
bioc_install("msigdbr")
bioc_install("org.Dm.eg.db")
bioc_install("pathview")
bioc_install("KEGGREST")

# libraries loading

library(clusterProfiler) # for the gene set analysis
library(msigdbr) # MSigDB gene sets in tidy format
library(org.Dm.eg.db) # Drosophila melanogaster annotations
library(magrittr) # for %>% piping
library(readxl) # to red excels
library(tidyverse)
library(pathview)
library(ggridges)
library(KEGGREST)
library('biomaRt')

#TISSUE DIFFERENCES #
# loading Drosophila orhologs symbols from Gff testes DEGs

gff2dmel <- read.csv("Gff2Dmel_proteins/Gff2Dmel.txt", sep="\t",header=F)
tissue.diff <- read_excel("glossina-negatives_annotated-significant-genes.xlsx")
gff.anno <- read.csv("glossina_annotation-full_table.tsv", sep="\t")
gff.anno <- gff.anno %>% dplyr::filter(grepl("LOC",Name)) %>% dplyr::select(gene_symbol, Name)
gff2dmel <- merge(gff2dmel, gff.anno, by.x = "V1", by.y = "gene_symbol") %>% 
  dplyr::rename(gff_symbol = V1, gff_id = Name, dmel_symbol = V2) %>% 
  dplyr::select(gff_symbol, gff_id, dmel_symbol)
neg.counts <- read_tsv("glossina_neg-counts.tsv")

neg.counts <- merge(gff2dmel, neg.counts, by.x = "gff_id", by.y = "gene")
tissue.diff <- merge(gff2dmel, tissue.diff, by.x = "gff_id", by.y = "gene")
ensembl <- useEnsembl(biomart = "ensembl", 
                      dataset = "dmelanogaster_gene_ensembl")
sig.tab.updated <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","entrezgene_id"),
                         filters = "ensembl_gene_id",
                         values = neg.counts$dmel_symbol,
                         mart = ensembl)
neg.counts <- merge(neg.counts, sig.tab.updated, by.x = "dmel_symbol", by.y = "ensembl_gene_id")
tissue.diff <- merge(tissue.diff, sig.tab.updated, by.x = "dmel_symbol", by.y = "ensembl_gene_id")
neg.counts <- neg.counts %>%
  dplyr::select(external_gene_name, everything()) %>%
  dplyr::select(-c(dmel_symbol, gff_id, gff_symbol, entrezgene_id)) %>%
  dplyr::distinct(external_gene_name, .keep_all = TRUE)
  
write_tsv(neg.counts, file = "glossina_neg-counts_dmel.tsv")

tissue.diff <- tissue.diff %>%
  dplyr::arrange(dplyr::desc(abs(log2FoldChange))) %>%
  dplyr::distinct(external_gene_name, .keep_all = TRUE)

write_tsv(tissue.diff, file = "Gff_tissue-diff-neg-sig.txt")

tissue.diff <- tissue.diff %>% 
  dplyr::filter(abs(log2FoldChange) > 2, padj < 0.05) %>%
  dplyr::mutate(regulation = ifelse(log2FoldChange>0, "upregulated", "downregulated")) %>%
  dplyr::select(gff_id, gff_symbol, dmel_symbol, log2FoldChange) %>% 
  dplyr::filter(!grepl("NA", gff_symbol)) %>% 
  na.omit()

dm_sets1 <- msigdbr(species = "Drosophila melanogaster", category = "C5")
dm_sets2 <- msigdbr(species = "Drosophila melanogaster", category = "C2")

# gene set enrichment analysis
any(duplicated(tissue.diff$dmel_symbol)) # check duplicates

tissue.diff <- tissue.diff %>%
  dplyr::arrange(dplyr::desc(abs(log2FoldChange))) %>%
  dplyr::distinct(dmel_symbol, .keep_all = TRUE)

tissue.diff.up <- tissue.diff %>% 
  dplyr::filter(log2FoldChange > 2)

tissue.diff.down <- tissue.diff %>% 
  dplyr::filter(log2FoldChange < -2) 

gsea_results <- function(gene_set, lfc_vector){
  return(
    GSEA(
      geneList = lfc_vector, # Ordered ranked gene list
      minGSSize = 2, # Minimum gene set size
      maxGSSize = 1000, # Maximum gene set set
      pvalueCutoff = 0.1, # p-value cutoff
      eps = 0, # Boundary for calculating the p value
      seed = TRUE, # Set seed to make results reproducible
      pAdjustMethod = "BH", # Benjamini-Hochberg correction
      TERM2GENE = dplyr::select(
        gene_set,
        gs_name,
        ensembl_gene
      )
    )
  )
}

do_gsea <- function(df){
  # Let's create a named vector ranked based on the log2 fold change values
  lfc_vector <- df$log2FoldChange
  names(lfc_vector) <- df$dmel_symbol
  
  # We need to sort the log2 fold change values in descending order here
  lfc_vector <- sort(lfc_vector, decreasing = TRUE)
  
  # Set the seed so our results are reproducible:
  set.seed(2020)
  
  dmel.BP <- msigdbr(species = "Drosophila melanogaster", category = "C5", subcategory="GO:BP")
  dmel.CC <- msigdbr(species = "Drosophila melanogaster", category = "C5", subcategory="GO:CC")
  dmel.MF <- msigdbr(species = "Drosophila melanogaster", category = "C5", subcategory="GO:MF")
  
  
  gsea_BP <- data.frame(gsea_results(dmel.BP, lfc_vector)@result)
  gsea_MF <- data.frame(gsea_results(dmel.MF, lfc_vector)@result)
  
  if (dim(gsea_MF)[1] == 0 & dim(gsea_BP)[1] == 0){
    res <- c()
  } 
  else if (dim(gsea_BP)[1] == 0){
    res <- gsea_MF
  }
  else if (dim(gsea_MF)[1] == 0){
    res <- gsea_BP
  } 
  else {
    res  <- rbind(gsea_BP, gsea_MF)
  }
  return(res)
}

res <- do_gsea(tissue.diff)
res_up <- do_gsea(as.data.frame(tissue.diff.up))
res_down <- do_gsea(tissue.diff.down)
# write_tsv(res, file = "GSEA_Gff_tissue-diff-neg.txt")


res <- read_tsv("GSEA_Gff_tissue-diff-neg.txt")
res$ID <- gsub("GOMF_", "", res$ID) 
res$ID <- gsub("GOBP_", "", res$ID) 
res$ID <- gsub("_", " ", res$ID)
png('GSEA_Gff_tissue-diff-neg.png', width=30,units="in", height=25, res=300)
ggplot(res, aes(x=reorder(ID, enrichmentScore), y=enrichmentScore, fill=regulation)) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("GO terms") +
  ylab("Enrichment") +
  ggtitle("Processes and functions different between Testes and MAGs") +
  scale_y_continuous(breaks = round(seq(min(res$enrichmentScore), max(res$enrichmentScore), by = 2), 1)) +
  theme_bw(base_size=26) +
  theme(
    legend.position='right',
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=28, vjust=1, hjust=0.5, face="bold"),
    axis.text.x=element_text(angle=0, size=0, hjust=1.10),
    axis.text.y=element_text(angle=0, size=30,  vjust=0.5),
    axis.title=element_text(size=30),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=30),  #Text size
    title=element_text(size=30)) +
  guides(fill=guide_legend(title="Regulation",override.aes=list(size=2.5))) +
  coord_flip()+
  scale_fill_manual(values=alpha(c("royalblue3","red3"),.3))

dev.off()

res[4,]$ID <- "RNA POLYMERASE II SPECIFIC DNA \nBINDING TRANSCRIPTION FACTOR BINDING"
png('Fig_2c.png', width=26,units="in", height=10, res=300)
ggplot(res, aes(x=reorder(ID, enrichmentScore), y=enrichmentScore, fill=regulation)) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("GO terms") +
  ylab("Enrichment") +
  scale_y_continuous(breaks = round(seq(-1,1, by = 0.2), 1)) +
  theme_bw(base_size=26) +
  theme(
    legend.position='right',
    legend.background=element_rect(),
    axis.text.x=element_text(angle=0, size=30, hjust=0.5),
    axis.text.y=element_text(angle=0, size=30,  vjust=0.5),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=30),  #Text size
    title=element_text(size=30)) +
  guides(fill=guide_legend(title="Regulation",override.aes=list(size=2.5))) +
  coord_flip()+
  scale_fill_manual(values=alpha(c("royalblue3","red3"),.3))

dev.off()

# cytoscape results



#SPIRO DIFFERENCES 
#### PART 1 #####
# loading Drosophila orhologs symbols from Gff testes DEGs

genes.df_1.5 <- read_xlsx("Tabelle_geni_diversi_fold-change.xlsx", sheet="l2FC1.5")
dmel.genes.df_1.5 <- genes.df_1.5 %>% dplyr::select(ID, Symbol, flybase, log2FoldChange) %>% filter(!grepl("NA",Symbol)) %>% na.omit()
genes.df_1 <- read_xlsx("Tabelle_geni_diversi_fold-change.xlsx", sheet="l2FC1")
dmel.genes.df_1 <- genes.df_1 %>% dplyr::select(ID, Symbol, flybase, log2FoldChange) %>% filter(!grepl("NA",Symbol)) %>% na.omit()
dmel.genes.df_1 <- rbind(dmel.genes.df_1.5,dmel.genes.df_1)
genes.df_0.5 <- read_xlsx("Tabelle_geni_diversi_fold-change.xlsx", sheet="l2FC0.5")
dmel.genes.df_0.5 <- genes.df_0.5 %>% dplyr::select(ID, Symbol, flybase, log2FoldChange) %>% filter(!grepl("NA",Symbol)) %>% na.omit()
dmel.genes.df_0.5 <- rbind(dmel.genes.df_1,dmel.genes.df_0.5)

dm_sets1 <- msigdbr(species = "Drosophila melanogaster", category = "C5")
dm_sets2 <- msigdbr(species = "Drosophila melanogaster", category = "C2")

# gene set enrichment analysis
any(duplicated(dmel.genes.df_0.5$flybase)) # check duplicates
# if false continue otherwise filter (see website)

# Let's create a named vector ranked based on the log2 fold change values
lfc_vector <- dmel.genes.df_0.5$log2FoldChange
names(lfc_vector) <- dmel.genes.df_0.5$flybase

# We need to sort the log2 fold change values in descending order here
lfc_vector <- sort(lfc_vector, decreasing = TRUE)

# Set the seed so our results are reproducible:
set.seed(2020)

gsea_results <- GSEA(
  geneList = lfc_vector, # Ordered ranked gene list
  minGSSize = 5, # Minimum gene set size
  maxGSSize = 1000, # Maximum gene set set
  pvalueCutoff = 1, # p-value cutoff
  eps = 0, # Boundary for calculating the p value
  seed = TRUE, # Set seed to make results reproducible
  pAdjustMethod = "BH", # Benjamini-Hochberg correction
  TERM2GENE = dplyr::select(
    dm_sets2,
    gs_name,
    gene_symbol
  )
)
gsea_result_df <- data.frame(gsea_results@result)
  

#### PART 2 ####

# loding file with orthologs ids form Gff to Dmel
gff2dmel <- read_tsv("./Gff2Dmel_proteins/Gff2Dmel.txt", col_names = c("Gff","Dmel","evalue"))
gff2dmel <-gff2dmel %>% 
  group_by(Gff) %>% 
  arrange(evalue) %>% 
  distinct(Gff, .keep_all = TRUE) %>% 
  ungroup()

# attach flybase id to all genes and relative log2fc 
l2fc.df <- read_xlsx("glossina_testes-design_de-results-complete.xlsx", sheet = "all_genes")
l2fc.df <- l2fc.df %>%
  left_join(gff2dmel %>% dplyr::select(Gff, Dmel, evalue),
            by = c("gene_id" = "Gff"))

dm_sets1 <- msigdbr(species = "Drosophila melanogaster", category = "C5") # GO
dm_sets2 <- msigdbr(species = "Drosophila melanogaster", category = "C2") # 

library(biomaRt)
ensembl <- useEnsembl(biomart = "ensembl", 
                      dataset = "dmelanogaster_gene_ensembl")
dmel_genes <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","entrezgene_id"),
                    filters = "ensembl_gene_id",
                    values = l2fc.df$Dmel,
                    mart = ensembl)
dmel_genes <- as_tibble(dmel_genes)

# gene set enrichment analysis
any(duplicated(l2fc.df$Dmel)) # check duplicates
# filter duplicates
l2fc.df <-l2fc.df %>% 
  group_by(Dmel) %>% 
  arrange(padj) %>% 
  distinct(Dmel, .keep_all = TRUE) %>% 
  ungroup()

# add Dmel gene symbols
l2fc.df <- l2fc.df %>%
  left_join(dmel_genes,
            by = c("Dmel" = "ensembl_gene_id"))

l2fc.df <-l2fc.df %>% 
  group_by(external_gene_name) %>% 
  arrange(padj) %>% 
  distinct(external_gene_name, .keep_all = TRUE) %>% 
  ungroup()

# Let's create a named vector ranked based on the log2 fold change values
lfc_vector <- l2fc.df$log2FoldChange
names(lfc_vector) <- l2fc.df$external_gene_name

# We need to sort the log2 fold change values in descending order here
lfc_vector <- sort(lfc_vector, decreasing = TRUE)

# Set the seed so our results are reproducible:
set.seed(2020)

gsea_results <- GSEA(
  geneList = lfc_vector, # Ordered ranked gene list
  TERM2GENE = dplyr::select(
    dm_sets1,
    gs_name,
    gene_symbol
  )
)
gsea_result_df <- data.frame(gsea_results@result)

# now let's see GSEA considering only significant genes
l2fc.sig <- l2fc.df %>% filter(padj < 0.05 & abs(log2FoldChange)>=0.5) 
# Let's create a named vector ranked based on the log2 fold change values
lfc_vector.sig <- l2fc.sig$log2FoldChange
names(lfc_vector.sig) <- l2fc.sig$external_gene_name

# We need to sort the log2 fold change values in descending order here
lfc_vector.sig <- sort(lfc_vector.sig, decreasing = TRUE)

# Set the seed so our results are reproducible:
set.seed(2020)

gsea_results.sig <- GSEA(
  geneList = lfc_vector.sig, # Ordered ranked gene list
  pvalueCutoff = 1,
  minGSSize = 10,
  pAdjustMethod = "fdr",
  TERM2GENE = dplyr::select(
    dm_sets1,
    gs_name,
    gene_symbol
  )
)
gsea_result_df.sig <- data.frame(gsea_results.sig@result)


#### PART 3 #### based on https://rpubs.com/jrgonzalezISGlobal/enrichment 
l2fc.sig <- l2fc.df %>% filter(padj < 0.05 & abs(log2FoldChange) > 0.5)
deGenes <- na.omit(l2fc.sig$entrezgene_id)
geneUniverse <- na.omit(l2fc.df$entrezgene_id)
ans.go <- enrichGO(gene = deGenes, 
                   keyType = "ENTREZID",
                   OrgDb ="org.Dm.eg.db",
                   ont = "all", 
                   universe = geneUniverse,
                   readable=TRUE,
                   pvalueCutoff = 0.05)
tab.go <- as.data.frame(ans.go)

dotplot(ans.go, showCategory = 15)


deGenes <- l2fc.sig %>% filter(entrezgene_id != "") %>% dplyr::select(entrezgene_id) %>% na.omit() %>% apply(.,1:2,as.character)
geneUniverse <- l2fc.df %>% filter(entrezgene_id != "") %>% dplyr::select(entrezgene_id) %>% na.omit() %>% apply(.,1:2,as.character)
ans.kegg <- enrichKEGG(gene = deGenes,
                       keyType = "ncbi-geneid",
                       organism = 'dme',
                       universe = geneUniverse,
                       pvalueCutoff = 0.05)
tab.kegg <- as.data.frame(ans.kegg)

dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)



#### PART 4 #### gseGO based on https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/
# first on the whole dataset
l2fc.df <- read_xlsx("glossina_testes-design_de-results-complete.xlsx", sheet = "all_genes")
l2fc.df <- l2fc.df %>%
  left_join(gff2dmel %>% dplyr::select(Gff, Dmel, evalue),
            by = c("gene_id" = "Gff"))
dmel_genes <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","entrezgene_id"),
                    filters = "ensembl_gene_id",
                    values = l2fc.df$Dmel,
                    mart = ensembl)
dmel_genes <- as_tibble(dmel_genes)

# gene set enrichment analysis
any(duplicated(l2fc.df$Dmel)) # check duplicates
# filter duplicates
l2fc.df <-l2fc.df %>% 
  group_by(Dmel) %>% 
  arrange(padj) %>% 
  distinct(Dmel, .keep_all = TRUE) %>% 
  ungroup()

# add Dmel gene symbols
l2fc.df <- l2fc.df %>%
  left_join(dmel_genes,
            by = c("Dmel" = "ensembl_gene_id"))

kegg_gene_list <- l2fc.df$log2FoldChange
names(kegg_gene_list) <- l2fc.df$entrezgene_id
kegg_gene_list<-na.omit(kegg_gene_list)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

kegg_organism = "dme"
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")

dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
emapplot(kk2)
cnetplot(kk2, categorySize="pvalue", foldChange=kegg_gene_list)
ridgeplot(kk2) + labs(x = "enrichment distribution")
dme <- pathview(gene.data=kegg_gene_list, pathway.id="dme00785", species = kegg_organism)
dme <- pathview(gene.data=kegg_gene_list, pathway.id="dme00785", species = kegg_organism, kegg.native = F)
knitr::include_graphics("dme00785.pathview.png")

# second on just significant genes
kegg_gene_list <- l2fc.sig$log2FoldChange
names(kegg_gene_list) <- l2fc.sig$entrezgene_id
kegg_gene_list<-na.omit(kegg_gene_list)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

kegg_organism = "dme"
kk3 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "fdr",
               keyType       = "ncbi-geneid")
# no terms found

# dotplot(kk3, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
# emapplot(kk3)
# cnetplot(kk3, categorySize="pvalue", foldChange=kegg_gene_list)
# ridgeplot(kk3) + labs(x = "enrichment distribution")
# dme <- pathview(gene.data=kegg_gene_list, pathway.id="dme00785", species = kegg_organism)
# dme <- pathview(gene.data=kegg_gene_list, pathway.id="dme00785", species = kegg_organism, kegg.native = F)
# knitr::include_graphics("dme00785.pathview.png")

# gsetGO
gene_list <- l2fc.sig$log2FoldChange
names(gene_list) <- l2fc.sig$entrezgene_id
gene_list<- gene_list %>% na.omit() %>% sort(decreasing = TRUE)

gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             OrgDb = "org.Dm.eg.db")
dotplot(gse, showCategory=20, split=".sign")

kegg_gene_list <- kegg_gene_list[abs(kegg_gene_list) > 0.5]
# let's see the pathways that might be influenced by our DEGs
did <- "dme01100"
dme <- pathview(gene.data=kegg_gene_list, pathway.id=did, species = kegg_organism)
knitr::include_graphics(paste0(did,".pathview.png"))

dme <- pathview(gene.data=kegg_gene_list, pathway.id="dme00785", species = kegg_organism)
knitr::include_graphics("dme00785.pathview.png")

dme <- pathview(gene.data=kegg_gene_list, pathway.id="dme03040", species = kegg_organism)
knitr::include_graphics("dme03040.pathview.png")

ctrl = F
if (ctr){
  dme_identifiers <- names(keggList("pathway", organism = "dme"))
  for (dme.path in dme_identifiers) {
    dme <- try(pathview(gene.data = kegg_gene_list, pathway.id = dme.path, species = "dme"), silent = TRUE)
    if (!inherits(dme, "try-error")) {
      outfile <- paste0(dme.path, ".pathview.png")
      png(outfile)
      print(dme)
      dev.off()
      knitr::include_graphics(outfile)
    }
  }
}

#### FINAL CONCLUSION ####

# ORA gives something, GSEA nothing
#### PART 3 #### based on https://rpubs.com/jrgonzalezISGlobal/enrichment 
pos <- l2fc.sig[l2fc.sig$log2FoldChange > 0.5,]
neg <- l2fc.sig[l2fc.sig$log2FoldChange < -0.5,]
upGenes <- pos %>% filter(entrezgene_id != "") %>% dplyr::select(entrezgene_id) %>% na.omit() %>% apply(.,1:2,as.character)
downGenes <- neg %>% filter(entrezgene_id != "") %>% dplyr::select(entrezgene_id) %>% na.omit() %>% apply(.,1:2,as.character)
geneUniverse <-  l2fc.df %>% filter(entrezgene_id != "") %>% dplyr::select(entrezgene_id) %>% na.omit() %>% apply(.,1:2,as.character)
ans.go.up <- enrichGO(gene = upGenes, 
                   keyType = "ENTREZID",
                   OrgDb ="org.Dm.eg.db",
                   ont = "BP", 
                   universe = geneUniverse,
                   readable=TRUE,
                   pvalueCutoff = 0.05,
                   pool=F)
tab.go.up <- as.data.frame(ans.go.up)
ans.go.down <- enrichGO(gene = downGenes, 
                      keyType = "ENTREZID",
                      OrgDb ="org.Dm.eg.db",
                      ont = "BP", 
                      universe = geneUniverse,
                      readable=TRUE,
                      pvalueCutoff = 0.05,
                      pool=F)
tab.go.down <- as.data.frame(ans.go.down)

dotplot(ans.go.down, showCategory = 15)
enrichMap(ans.go.down, vertex.label.cex=1.2, layout=igraph::layout.kamada.kawai)
plotGOgraph(ans.go.down)

gse_downGenes <- neg$log2FoldChange
names(gse_downGenes) <- neg$entrezgene_id
gse_downGenes <- gse_downGenes %>% na.omit() %>% sort(decreasing = TRUE)

gse <- gseGO(geneList=gse_downGenes, ont="all", OrgDb=org.Dm.eg.db, verbose=F, pvalueCutoff = 01)

