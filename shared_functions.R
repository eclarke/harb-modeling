
#' Common functions for MIRA 16S analysis [M01-01]
#' @author Erik Clarke

#' Loads the MIRA datasets into a phyloseq object
#' @param seqtab_fp the DADA2 sequence table
#' @param taxa_fp the DADA2 taxonomy table
#' @param specimen_fp specimen metadata table
#' @return a phyloseq object
load_mira_data <- function(
  seqtab_fp="seqtab.rds", taxa_fp="taxa.rds", 
  meds_fp="MIRA_Medications_Table.csv", specimen_fp="MIRA_Specimen_Table.csv")
{
  # Load sequence table
  seqtab <- readRDS(seqtab_fp)
  # Load taxonomy table
  taxa <- readRDS(taxa_fp)
  # Load medication table
  meds <- readr::read_csv(meds_fp) %>%
    mutate()
  # Load sample metadata
  samples <- readr::read_csv(specimen_fp) %>%
    mutate_at(
      vars(
        specimen_id2,
        specimen_id,
        specimen_type,
        subject_id,
        flow_cell_id),
      funs(factor)
    ) %>%
    as.data.frame()
  # Ensure the sequence table and metadata have the same samples
  rownames(samples) <- samples$specimen_id2
  seqtab2 <- seqtab[rownames(seqtab) %in% rownames(samples),]
  samples2 <- samples[rownames(samples) %in% rownames(seqtab2),]
  seqtab3 <- seqtab[match(rownames(seqtab2), rownames(samples2)),]
  taxa2 <- taxa[match(rownames(taxa),colnames(seqtab3)), ]
  # Isolate and give non-integer names to the sequence variants ("sv00001", etc)
  seqs <- data.frame(
    seqs=rownames(taxa2), 
    id=paste0("sv", str_pad(seq_along(rownames(taxa2)), width = 5, pad=0)))
  # Set this name in the corresponding tables for the variants
  rownames(taxa2) <- seqs$id
  colnames(seqtab3) <- seqs$id
  list(
    ps = phyloseq(
      otu_table(seqtab3, taxa_are_rows = FALSE),
      sample_data(samples2),
      tax_table(taxa2)),
    samples = as_tibble(samples),
    seqs = as_tibble(seqs),
    meds = as_tibble(meds)
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

#' Calculates a "window" of effectiveness from a logical vector.
#' It proceeds stepwise through the vector; if the vector is true at that idx,
#' the new vector has every element from that index to that index + lag
#' marked as true.
#' @param lag how long the window lasts after x[i] == TRUE
effective_window <- function(x, lag=1) {
  l <- length(x)
  x2 <- rep(FALSE, l)
  for (i in 1:l) {
    if (!is.na(x[i]) & x[i]) {
      end = i + lag
      x2[i:end] <- TRUE
    }
  }
  x2[1:l]
}


pp_subject_1.5.0 <- function(post, subject, data, use_a_spec=F) {
#  browser()
  .d <- data
  .d$specimen_idx <- as.integer(as.factor(.d$specimen_id2))
  ps <- post %>% filter(subject_id == subject)
  stopifnot(nrow(ps) > 0)
  ps$specimen_idx <- as.integer(as.factor(ps$specimen_id2))
  cols <- n_distinct(ps$specimen_id2)
  rows <- nrow(ps)/cols
  pred_prob <- matrix(data=NA, nrow=rows, ncol=cols) 
  colnames(pred_prob) <- unique(ps$specimen_id2)
  
  p.link <- function(d, n, idx, prev_prop) {
    d2 <- filter(d, specimen_idx == idx)
    prob <- logistic(with(d2, mu + a_subj + a_spec*use_a_spec + b_lag * prev_prop + b_abx * on_abx))
    rbinom(n, unique(d2$total_reads), prob)/unique(d2$total_reads)
  }
  
  pred_prob[,1] <- p.link(ps, rows, 1, .d$lag_emp_prop[1])
  
  # Skip autoregression if only one timepoint
  if (cols > 1) {
    for (i in 2:cols) {
      pred_prob[,i] <- p.link(ps, rows, i, pred_prob[,i-1])
    }
  }
  
  pred.df <- gather(data.frame(pred_prob), key = specimen_id2)
  # levels(pred.df$key) <- levels(as.factor(ps$specimen_id2))

  set.seed(12)
  groups <- sample(1:rows, 100)
  pred.df.sub <- group_by(pred.df, specimen_id2) %>%
    mutate(group = seq_along(value)) %>%
    group_by(group) %>%
    filter(group %in% groups) %>%
    #filter(group %in% sample(1:n_distinct(group), 150)) %>%
    rename(pred_prob=value) %>%
    left_join(.d) %>%
    mutate(emp_prob = read_count/total_reads)
  
  dp = pred.df.sub %>% ungroup %>% distinct(specimen_id2, .keep_all = T)
  
  ggplot(pred.df.sub, aes(study_day, pred_prob, group=group)) +
    geom_line(alpha=0.1) +
    geom_line(data=dp, aes(y=emp_prob, group=NULL), color="red") +
    geom_vline(data=dp, aes(xintercept=study_day, alpha=on_abx), linetype=2) +
    scale_alpha_manual(values=c("TRUE"=0.6, "FALSE"=0), guide=FALSE) +
    scale_y_continuous(expand=c(0,0), labels=scales::percent) +
    ggtitle(subject)
}