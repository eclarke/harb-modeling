
#' Common functions for MIRA 16S analysis [M01-01]
#' @author Erik Clarke

#' Loads the MIRA datasets into a phyloseq object
#' @param seqtab_fp the DADA2 sequence table
#' @param taxa_fp the DADA2 taxonomy table
#' @param specimen_fp specimen metadata table
#' @return a phyloseq object
load_mira_data <- function(
  seqtab_fp="seqtab.rds", taxa_fp="taxa.rds", specimen_fp="MIRA_Specimen_Table.csv")
{
  seqtab <- readRDS(seqtab_fp)
  taxa <- readRDS(taxa_fp)
  samples <- readr::read_csv(specimen_fp) %>%
    mutate_at(
      vars(specimen_id2, specimen_id, specimen_type, subject_id, flow_cell_id), 
      funs(factor)) %>%
    filter(!duplicated(specimen_id2)) %>%
    as.data.frame()
  rownames(samples) <- samples$specimen_id2
  seqtab2 <- seqtab[rownames(seqtab) %in% rownames(samples),]
  samples2 <- samples[rownames(samples) %in% rownames(seqtab2),]
  seqtab3 <- seqtab[match(rownames(seqtab2), rownames(samples2)),]
  taxa2 <- taxa[match(rownames(taxa),colnames(seqtab3)), ]
  seqs <- data.frame(seqs=rownames(taxa2), id=paste0("sv", seq_along(rownames(taxa2))))
  rownames(taxa2) <- seqs$id
  colnames(seqtab3) <- seqs$id
  list(
    ps = phyloseq(
      otu_table(seqtab3, taxa_are_rows = FALSE),
      sample_data(samples2),
      tax_table(taxa2)
    ),
    seqs = seqs
  )
}

#' Convert a phyloseq object to an agglomerated (melted) dataframe.
#' @param ps a phyloseq object
#' @param sample.col.name the sample_data() column corresponding rownames(otu_table()), as a string
#' @param otu_col the desired name of the otu id column
#' @param count_col the desired name of the count column
phyloseq_to_agglomerated <- function(ps, sample_col, otu_col, count_col) {
  otus <- otu_table(ps)
  otus <- otus[, colSums(otus) > 0]
  taxa <- tax_table(ps)
  taxa <- taxa[rownames(taxa) %in% colnames(otus), ]
  taxa.df <- as.data.frame(taxa@.Data)
  taxa.df$otu_id <- rownames(taxa)
  agg <- reshape2::melt(otus, varnames=c(sample_col, otu_col), value.name=count_col) %>%
    left_join(sample_data(ps)) %>%
    left_join(taxa.df)
}

#' Convert an agglomerated dataframe to a phyloseq object.
#' @param agg the agglomerated dataframe
#' @param otu_col the unquoted name of the otu id column in agg
#' @param sample_col the unquoted name of the sample id column in agg
#' @param count_col the unquoted name of the counts column in agg
#' @param ... remaining unquoted names of taxonomic ranks (i.e. Kingdom:Species)
agglomerated_to_phyloseq <- function(agg, otu_col, sample_col, count_col, ...) {
  otu_col <- enquo(otu_col)
  sample_col <- enquo(sample_col)
  count_col <- enquo(count_col)
  taxa_ranks <- quos(...)

  # Rebuild count matrix
  otu <- select(agg, !! otu_col, !! sample_col, !! count_col) %>%
    tidyr::spread(!! sample_col, !! count_col, fill=0) %>%
    filter(!is.na(!! otu_col)) %>%
    as.data.frame()
  rownames(otu) <- otu[[quo_name(otu_col)]]
  otu[[quo_name(otu_col)]] <- NULL
  otu.mat <- t(as.matrix(otu))
  
  # Rebuild taxonomy table
  taxa <- select(agg, !! otu_col, !!! taxa_ranks) %>%
    distinct() %>%
    as.data.frame()
  rownames(taxa) <- taxa[[quo_name(otu_col)]]
  taxa[[quo_name(otu_col)]] <- NULL
  taxa <- as.matrix(taxa)
  taxa <- taxa[match(rownames(taxa), colnames(otu.mat)),]
  
  # Rebuild sample table
  samples <- select(agg, -c(!! otu_col, !! count_col, !!! taxa_ranks)) %>%
    distinct() %>%
    as.data.frame()
  rownames(samples) <- samples[[quo_name(sample_col)]]
  samples <- samples[match(rownames(samples), rownames(otu.mat)), ]
  
  phyloseq(
    otu_table(otu.mat, taxa_are_rows = FALSE),
    sample_data(samples),
    tax_table(taxa)
  )
}
