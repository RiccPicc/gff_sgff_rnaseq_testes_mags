# read GFF file with annotations in order to map our genes
library("tidyverse")
library("ape")

gff_file <- read.gff("GCF_014805625.2_Yale_Gfus_2_genomic.gff")

library(stringr)

attr_names <- gff_file$attributes[1]
attr_names <- sapply(strsplit(attr_names, ";"), "[")
attr_names <- sapply(strsplit(attr_names, "="), "[", 1)

attr_df <- gff_file$attributes %>%
  as_tibble()

attr_df <- attr_df %>%
  separate_rows(value, sep = ';') %>%
  separate(value, into=c("var","val"), sep = '=')

col_names <- unique(attr_df$var)

attr_df <- attr_df %>% mutate(ID = ifelse(var == "ID", val, NA)) %>% fill(ID)

# df <- attr_df %>%
#   group_by(ID) %>%
#   mutate(across(col_names, ~ ifelse(var == .x, val, NA))) %>%
#   ungroup() %>%
#   select(-var, -val)

df <- attr_df %>%
  group_by(ID) %>%
  mutate(
    Dbxref = ifelse(var == "Dbxref", val, NA),
    Name = ifelse(var == "Name", val, NA),
    chromosome = ifelse(var == "chromosome", val, NA),
    collection_date = ifelse(var == "collection-date", val, NA),
    gbkey = ifelse(var == "gbkey", val, NA),
    genome = ifelse(var == "genome", val, NA),
    isolate = ifelse(var == "isolate", val, NA),
    isolation_source = ifelse(var == "isolation-source", val, NA),
    mol_type = ifelse(var == "mol_type", val, NA),
    tissue_type = ifelse(var == "tissue-type", val, NA),
    gene = ifelse(var == "gene", val, NA),
    gene_biotype = ifelse(var == "gene_biotype", val, NA),
    Parent = ifelse(var == "Parent", val, NA),
    model_evidence = ifelse(var == "model_evidence", val, NA),
    product = ifelse(var == "product", val, NA),
    transcript_id = ifelse(var == "transcript_id", val, NA),
    protein_id = ifelse(var == "protein_id", val, NA),
    Note = ifelse(var == "Note", val, NA),
    exception = ifelse(var == "exception", val, NA),
    anticodon = ifelse(var == "anticodon", val, NA),
    inference = ifelse(var == "inference", val, NA),
    partial = ifelse(var == "partial", val, NA),
    start_range = ifelse(var == "start_range", val, NA),
    transl_except = ifelse(var == "transl_except", val, NA),     
    end_range = ifelse(var == "end_range", val, NA),
    pseudo = ifelse(var == "pseudo", val, NA),
    Target = ifelse(var == "Target", val, NA),
    for_remapping = ifelse(var == "for_remapping", val, NA),
    gap_count = ifelse(var == "gap_count", val, NA),
    num_ident = ifelse(var == "num_ident", val, NA),     
    num_mismatch = ifelse(var == "num_mismatch", val, NA),
    pct_coverage = ifelse(var == "pct_coverage", val, NA),
    pct_coverage_hiqual = ifelse(var == "pct_coverage_hiqual", val, NA),
    pct_identity_gap = ifelse(var == "pct_identity_gap", val, NA),
    pct_identity_ungap = ifelse(var == "pct_identity_ungap", val, NA), 
    rank = ifelse(var == "rank", val, NA),
    Gap = ifelse(var == "Gap", val, NA)
  ) %>%
  select(-var, -val)

df <- df %>% 
  group_by(ID) %>%
  summarize(Dbxref = first(na.omit(Dbxref)),
            Name = first(na.omit(Name)),
            chromosome = first(na.omit(chromosome)),
            collection_date = first(na.omit(collection_date)),
            gbkey = first(na.omit(gbkey)),
            genome = first(na.omit(genome)),
            isolate = first(na.omit(isolate)),
            isolation_source = first(na.omit(isolation_source)),
            mol_type = first(na.omit(mol_type)),
            tissue_type = first(na.omit(tissue_type)),
            gene = first(na.omit(gene)),
            gene_biotype = first(na.omit(gene_biotype)),
            Parent = first(na.omit(Parent)),
            model_evidence = first(na.omit(model_evidence)),
            product = first(na.omit(product)),
            transcript_id = first(na.omit(transcript_id)),
            protein_id = first(na.omit(protein_id)),
            Note = first(na.omit(Note)),
            exception = first(na.omit(exception)),
            anticodon = first(na.omit(anticodon)),
            inference = first(na.omit(inference)),
            partial = first(na.omit(partial)),
            start_range = first(na.omit(start_range)),
            transl_except = first(na.omit(transl_except)),     
            end_range = first(na.omit(end_range)),
            pseudo = first(na.omit(pseudo)),
            Target = first(na.omit(Target)),
            for_remapping = first(na.omit(for_remapping)),
            gap_count = first(na.omit(gap_count)),
            num_ident = first(na.omit(num_ident)),     
            num_mismatch = first(na.omit(num_mismatch)),
            pct_coverage = first(na.omit(pct_coverage)),
            pct_coverage_hiqual = first(na.omit(pct_coverage_hiqual)),
            pct_identity_gap = first(na.omit(pct_identity_gap)),
            pct_identity_ungap = first(na.omit(pct_identity_ungap)), 
            rank = first(na.omit(rank)),
            Gap = first(na.omit(Gap)))

write_tsv(df, file = "glossina_annotation_table.tsv")
library(writexl)
write_xlsx(df, path = "glossina_annotation_table.xlsx")

save.image("annotation_table.RData")

##########################################
gff_file <- read.gff("VectorBase-63_GfuscipesIAEA2018.gff")

attr_names <- gff_file$attributes[1]
attr_names <- sapply(strsplit(attr_names, ";"), "[")
attr_names <- sapply(strsplit(attr_names, "="), "[", 1)

attr_df <- gff_file$attributes %>%
  as_tibble()

attr_df <- attr_df %>%
  separate_rows(value, sep = ';') %>%
  separate(value, into=c("var","val"), sep = '=')

col_names <- unique(attr_df$var)

attr_df <- attr_df %>% mutate(ID = ifelse(var == "ID", val, NA)) %>% fill(ID)

df <- attr_df %>%
  group_by(ID) %>%
  mutate(
    description = ifelse(var == "description", val, NA),
    ebi_biotype = ifelse(var == "ebi_biotype", val, NA),
    Parent = ifelse(var == "Parent", val, NA),
    gene_ebi_biotype = ifelse(var == "gene_ebi_biotype", val, NA),
    gene_id = ifelse(var == "gene_id", val, NA),
    protein_source_id = ifelse(var == "protein_source_id", val, NA),
    Name = ifelse(var == "Name", val, NA)) %>%
  dplyr::select(-var, -val)

df <- df %>% 
  group_by(ID) %>%
  dplyr::summarize(description = dplyr::first(na.omit(description)),
            ebi_biotype = dplyr::first(na.omit(ebi_biotype)),
            Parent = dplyr::first(na.omit(Parent)),
            gene_ebi_biotype = dplyr::first(na.omit(gene_ebi_biotype)),
            gene_id = dplyr::first(na.omit(gene_id)),
            protein_source_id = dplyr::first(na.omit(protein_source_id)),
            Name = dplyr::first(na.omit(Name)))

write_tsv(df, file = "glossina_annotation2_table.tsv")
write_xlsx(df, path = "glossina_annotation2_table.xlsx")

save.image("annotation2_table.RData")

# per mergiare le due annotazioni bisogna fre un gsub dell'id degli "exon"

anno2 <- df
anno1 <- read_tsv("glossina_annotation_table.tsv")
anno2$ID <- gsub("^exon_", "", anno2$ID)
anno2$ID <- gsub("_t\\d+-E\\d+$", "", anno2$ID)
anno2 <- unique(anno2[grepl("^LOC", anno2$ID), ]  %>%
                  rename(gene_symbol = gene_id) %>% 
                  dplyr::select(-Parent) %>% 
                  dplyr::select_if(~!all(is.na(.))))

anno.full <- unique(merge(anno1,anno2,by.x = "gene", by.y = "ID"))
write_tsv(anno.full, file = "glossina_annotation-full_table.tsv")
write_xlsx(anno.full, path = "glossina_annotation-full_table.xlsx")

pseudo <- unique(na.omit(df[df$ebi_biotype == "pseudogene",]%>% 
                    dplyr::select_if(~!all(is.na(.)))))
write_tsv(pseudo, file = "glossina_pseudo-annotation_table.tsv")
write_xlsx(pseudo, path = "glossina_pseudo-annotation_table.xlsx")


