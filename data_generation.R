#' ARB Data Generation Functions
#' 
#' Functions should be named with genXXXX, where XXXX is the version number.
#' 
#' Document the parameters and changes each makes to the structure, and the 
#' model numbers the generators work with.
#' 
#' Generic functions should not be versioned.


load_mira_phyloseq <- function(css_normalize=FALSE) {
  mira.all <- load_mira_data(
    seqtab_fp = here("../shared/data/seqtab.rds"),
    taxa_fp = here("../shared/data/taxa.rds"),
    meds_fp = here("../shared/data/MIRA_Medications_Table.csv"),
    specimen_fp = here("../shared/data/MIRA_Specimen_Table.csv")
  )
  
  seqs <- mira.all$seqs
  meds <- mira.all$meds
  
  mira <- mira.all$ps
  # Remove non-MIRA subjects from dataset (incl. blanks and controls)
  .nonmira <- is.na(sample_data(mira)$subject_id)
  print(sprintf("Removing %d non-MIRA samples...", sum(.nonmira)))
  mira <- prune_samples(!.nonmira, mira)
  # Remove culture samples (they break the subject/type/date unique constraint)
  .culture <- grepl("Culture", sample_data(mira)$specimen_type)
  print(sprintf("Removing %d culture samples...", sum(.culture)))
  mira <- prune_samples(!.culture, mira)
  # Remove empty samples
  .empty <- sample_sums(mira) == 0
  print(sprintf("Removing %d empty samples...", sum(.empty)))
  mira <- prune_samples(!.empty, mira)
  
  # Identify "duplicated" specimens (same subject, specimen type, and study day)
  sample_data(mira) <- sample_data(mira) %>% 
    group_by(subject_id, specimen_type, study_day) %>%
    # "specimen_id3" has the same value for duplicated specimens so phyloseq can 
    # use it as a grouping level
    mutate(specimen_id3 = case_when(
      n() > 1 ~ paste0(first(as.character(specimen_id2)), "_D"),
      TRUE ~ as.character(specimen_id2)
    )) %>%
    ungroup() %>% as.data.frame() %>%
    set_rownames(.$specimen_id2)
  
  # Sum abundances and merge sample table for duplicates
  mira <- phyloseq::merge_samples(mira, "specimen_id3", fun=sum)
  # Re-add the relevant metadata since merge_samples mangles factors and dates
  sample_data(mira) <- sample_data(mira) %>% 
    mutate(specimen_id2 = rownames(.)) %>%
    select(specimen_id2, specimen_id3) %>%
    left_join(mira.all$samples) %>%
    ungroup() %>% as.data.frame() %>%
    set_rownames(.$specimen_id2)
  
  # Restrict to only samples for which we have abx data
  .abx_specimens <- as.character(inner_join(sample_data(mira), meds)$specimen_id2)
  mira.abx <- prune_samples(
    sample_data(mira)$specimen_id2 %in% .abx_specimens, mira)
  
  # Normalize if requested
  if (css_normalize) {
    gt1_feature_samples <- rowSums(otu_table(mira) > 0) > 1
    message(sprintf("removing %d samples for having < 2 features", length(gt1_feature_samples)-sum(gt1_feature_samples)))
    mira.abx.gt1 <- prune_samples(names(which(gt1_feature_samples)), mira.abx)
    normalized <- t(metagenomeSeq::cumNormMat(t(mira.abx.gt1@otu_table@.Data)))
    otu_table(mira.abx.gt1) <- otu_table(normalized, taxa_are_rows = F)
    mira.abx <- mira.abx.gt1
  }
  
  return(list(mira.abx=mira.abx, meds=meds))
}

clean_data <- function(mira.abx, meds, .genus, .species) {
  # Converted to melted form
  agg <- phyloseq_to_agglomerated(mira.abx, "specimen_id2", "otu_id", "read_count") %>%
    group_by(specimen_id2) %>%
    mutate(total_reads = sum(read_count)) %>%
    ungroup()
  
  # Sanity checks to ensure specified lineage exists in data
  stopifnot(.genus %in% agg$Genus)
  if (!is.null(.species)) {
    stopifnot(.species %in% filter(agg, Genus==.genus)$Species)
  }
  
  # Clean up subject/timepoint data
  subjects <- sample_data(mira.abx) %>%
    group_by(subject_id) %>% 
    mutate(exit_date = max(collection_date)) %>%
    distinct(subject_id, enroll_date, exit_date) %>%
    right_join(meds) %>%
    group_by(subject_id) %>%
    mutate(study_day = as.integer(collection_date - enroll_date)) %>%
    mutate(exit_day = study_day[collection_date == exit_date]) %>%
    # Limit to only a week before enrollment date and nothing after
    filter(study_day > -3, collection_date <= exit_date)
  
  # Manually split MIRA_024 into two sub-subjects
  subjects <- subjects %>%
    mutate(subject_id2 = case_when(
      subject_id != "MIRA_024" ~ subject_id,
      study_day <= 33 ~ "MIRA_024a",
      study_day >= 73 ~ "MIRA_024b",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(subject_id2)) %>%
    mutate(study_day = case_when(
      subject_id2 == "MIRA_024b" ~ study_day - 73,
      TRUE ~ as.double(study_day)
    )) %>%
    mutate(exit_day = ifelse(subject_id2 == "MIRA_024b", exit_day-73, exit_day))
  
  # Clean up antibiotics data
  subjects <- subjects %>%
    ungroup() %>%
    separate_rows(abx_b) %>%
    # Handle coformulations
    mutate(abx_b2 = case_when(
      abx_b == "clavulanate" ~ "amoxicillin.clavulanate",
      abx_b == "sulbactam" ~ "ampicillin.sulbactam",
      abx_b == "tazobactam" ~ "piperacillin.tazobactam",
      abx_b == "trimethoprim" ~ "trimethoprim.sulfamethoxazole",
      TRUE ~ abx_b)) %>%
    group_by(subject_id, collection_date) %>%
    # Remove remaining part of coformulation if coformulation is present
    filter(
      (abx_b2 != "amoxicillin" & ("amoxicillin.clavulanate" %in% abx_b2)) | 
        !("amoxicillin.clavulanate" %in% abx_b2) |
        is.na(abx_b2)
    ) %>%
    filter(
      (abx_b2 != "ampicillin" & ("ampicillin.sulbactam" %in% abx_b2)) |
        !("ampicillin.sulbactam" %in% abx_b2) |
        is.na(abx_b2)
    ) %>%
    filter(
      (abx_b2 != "sulfamethoxazole" & ("trimethoprim.sulfamethoxazole" %in% abx_b2)) |
        !("trimethoprim.sulfamethoxazole" %in% abx_b2) |
        is.na(abx_b2)
    ) %>%
    ungroup() %>%
    filter((abx_b2 != "piperacillin") | is.na(abx_b2)) %>%
    mutate(abx_b = abx_b2) %>%
    select(-abx_b2) %>%
    mutate(abx_idx = as.numeric(as.factor(abx_b))) %>%
    mutate(abx_yn = ifelse(is.na(abx_b), 0, 1))
  
  if (is.null(.species)) {
    d <- agg %>%
      filter(!is.na(Genus)) %>%
      group_by(specimen_id2, total_reads, Kingdom, Phylum, Class, Order, Family, Genus) %>%
      summarize(read_count = sum(read_count)) %>%
      ungroup() %>%
      filter(Genus == .genus) %>%
      left_join(sample_data(mira.abx)) %>%
      select(-c(Kingdom:Family)) %>%
      droplevels()
  } else {
    d <- agg %>%
      filter(!is.na(Genus) & !is.na(Species)) %>%
      group_by(specimen_id2, total_reads, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
      summarize(read_count = sum(read_count)) %>%
      ungroup() %>%
      filter(Genus == .genus & Species == .species) %>%
      left_join(sample_data(mira.abx)) %>%
      select(-c(Kingdom:Family)) %>% 
      droplevels()
  }
  
  left_join(subjects, d)
}


gen01 <- function(d2, .specimen_type, min_times_abx=2, min_median_abx=2, min_reads_subj=1) {
  d3 <- d2 %>% filter(specimen_type==.specimen_type) %>%
    distinct(subject_id, specimen_id2, read_count, total_reads, study_day) %>%
    group_by(subject_id) %>%
    filter(sum(read_count) > min_reads_subj) %>%
    arrange(subject_id, study_day) %>%
    mutate(prev = lag(read_count)/lag(total_reads)) %>%
    right_join(d2) %>%
    filter(!is.na(prev)) %>%
    rename(reads=read_count, total=total_reads, specimens=specimen_id2, subjects=subject_id) %>%
    mutate(abx_b = ifelse(is.na(abx_b), "none", abx_b))
  
  abx_to_keep <- (d3 %>% group_by(abx_b, subjects) %>%
                    filter(abx_b != "none") %>%
                    summarize(n=n()) %>% 
                    mutate(n_subjects=n()) %>%
                    filter(n_subjects >= min_times_abx, median(n) >= min_median_abx) %>%
                    distinct(abx_b))$abx_b
  
  subjects_to_keep <- (filter(d3, abx_b %in% abx_to_keep) %>% distinct(subjects))$subjects 
  
  stopifnot(length(subjects_to_keep) > 0)
  
  d4 <- d3 %>% ungroup() %>%
    filter(subjects %in% subjects_to_keep)
  d_abx <- d4 %>%
    reshape2::dcast(subjects + specimens + reads + total + prev ~ abx_b, value.var="abx_yn", fill=0)
  # browser()
  result <- tidybayes::compose_data(d_abx[, c(1:5)])
  abx <- as.matrix(d_abx[, -c(1:5)])
  # Remove abxs seen less than the minimum, and the NA column
  abx <- abx[, abx_to_keep, drop=F]
  result$abx_b <- as.numeric(as.factor(colnames(abx)))
  result$n_abx <- ncol(abx)
  result$abx <- abx
  result$d4 <- d4
  result$d_abx <- d_abx[, c(1:5)]
  result
}

#' Frameshifts antibiotics by one day.
gen02 <- function(d, .specimen_type, min_times_abx=2, min_median_abx=2, min_reads_subj=1) {
  d <- d %>% filter(specimen_type==.specimen_type) %>%
    select(specimen_id2, subject_id, study_day, read_count, total_reads, abx_b)
  # Nest multiple antibiotics into a list-column
  d2 <- d %>% group_by(specimen_id2) %>%
    nest(abx_b, .key="abx_nested") %>%
    left_join(select(d, -abx_b)) %>% 
    distinct(specimen_id2, .keep_all = TRUE)
  # Remove subjects with fewer than the min. total reads for the target taxa
  d3 <- d2 %>% group_by(subject_id) %>%
    filter(sum(read_count) > min_reads_subj) %>%
    ungroup()
  message(sprintf(
    "Removed %d subjects for having <= %d reads on target", 
    n_distinct(d3$subject_id)-n_distinct(d2$subject_id), min_reads_subj))

  # Add lag columns for previous day's proportions and antibiotics
  d4 <- d3 %>%
    arrange(subject_id, study_day) %>%
    group_by(subject_id) %>%
    mutate(prev = lag(read_count)/lag(total_reads)) %>%
    mutate(abx_lag = lag(abx_nested, default = list(tibble(abx_b=NA)))) %>%
    unnest(abx_lag, .preserve=abx_nested) %>%
    rename(abx_lag = abx_b) %>%
    filter(!is.na(prev)) %>%
    rename(reads=read_count, total=total_reads, specimens=specimen_id2, subjects=subject_id) %>%
    mutate(abx_yn = ifelse(is.na(abx_lag), 0, 1)) %>% 
    mutate(abx_lag = ifelse(is.na(abx_lag), "none", abx_lag)) 

  
  abx_to_keep <- (d4 %>% group_by(abx_lag, subjects) %>%
                    filter(abx_lag != "none") %>%
                    summarize(n=n()) %>% 
                    mutate(n_subjects=n()) %>%
                    filter(n_subjects >= min_times_abx, median(n) >= min_median_abx) %>%
                    distinct(abx_lag))$abx_lag
  
  subjects_to_keep <- (filter(d4, abx_lag %in% abx_to_keep) %>% distinct(subjects))$subjects 
  
  stopifnot(length(subjects_to_keep) > 0)
  
  d5 <- d4 %>% ungroup() %>%
    filter(subjects %in% subjects_to_keep)
  d_abx <- d5 %>%
    reshape2::dcast(subjects + specimens + reads + total + prev ~ abx_lag, value.var="abx_yn", fill=0)
  result <- tidybayes::compose_data(d_abx[, c(1:5)])
  abx <- as.matrix(d_abx[, -c(1:5)])
  # Remove abxs seen less than the minimum, and the NA column
  abx <- abx[, abx_to_keep, drop=F]
  result$abx_b <- as.numeric(as.factor(colnames(abx)))
  result$n_abx <- ncol(abx)
  result$abx <- abx
  result$d5 <- d5
  result$d_abx <- d_abx[, c(1:5)]
  result
}

#' Frameshifts antibiotics by multiple days.
gen03 <- function(d, .specimen_type, min_times_abx=2, min_median_abx=2, min_reads_subj=1, .lag=1) {
  d <- d %>% filter(specimen_type==.specimen_type) %>%
    select(specimen_id2, subject_id, study_day, read_count, total_reads, abx_b)
  # Nest multiple antibiotics into a list-column
  d2 <- d %>% group_by(specimen_id2) %>%
    nest(abx_b, .key="abx_nested") %>%
    left_join(select(d, -abx_b)) %>% 
    distinct(specimen_id2, .keep_all = TRUE)
  # Remove subjects with fewer than the min. total reads for the target taxa
  d3 <- d2 %>% group_by(subject_id) %>%
    filter(sum(read_count) > min_reads_subj) %>%
    ungroup()
  message(sprintf(
    "Removed %d subjects for having <= %d reads on target", 
    n_distinct(d2$subject_id)-n_distinct(d3$subject_id), min_reads_subj))
  
  # Add lag columns for previous day's proportions and antibiotics
  d4 <- d3 %>%
    arrange(subject_id, study_day) %>%
    group_by(subject_id) %>%
    mutate(prev = lag(read_count)/lag(total_reads)) %>%
    mutate(abx_lag = lag(abx_nested, default = list(tibble(abx_b=NA)), lag=.lag)) %>%
    unnest(abx_lag, .preserve=abx_nested) %>%
    rename(abx_lag = abx_b) %>%
    filter(!is.na(prev)) %>%
    rename(reads=read_count, total=total_reads, specimens=specimen_id2, subjects=subject_id) %>%
    mutate(abx_yn = ifelse(is.na(abx_lag), 0, 1)) %>% 
    mutate(abx_lag = ifelse(is.na(abx_lag), "none", abx_lag)) 
  
  
  abx_to_keep <- (d4 %>% group_by(abx_lag, subjects) %>%
                    filter(abx_lag != "none") %>%
                    summarize(n=n()) %>% 
                    mutate(n_subjects=n()) %>%
                    filter(n_subjects >= min_times_abx, median(n) >= min_median_abx) %>%
                    distinct(abx_lag))$abx_lag
  
  subjects_to_keep <- (filter(d4, abx_lag %in% abx_to_keep) %>% distinct(subjects))$subjects 
  
  stopifnot(length(subjects_to_keep) > 0)
  
  d5 <- d4 %>% ungroup() %>%
    filter(subjects %in% subjects_to_keep)
  d_abx <- d5 %>%
    reshape2::dcast(subjects + specimens + reads + total + prev ~ abx_lag, value.var="abx_yn", fill=0)
  result <- tidybayes::compose_data(d_abx[, c(1:5)])
  abx <- as.matrix(d_abx[, -c(1:5)])
  # Remove abxs seen less than the minimum, and the NA column
  abx <- abx[, abx_to_keep, drop=F]
  result$abx_b <- as.numeric(as.factor(colnames(abx)))
  result$n_abx <- ncol(abx)
  result$abx <- abx
  result$d5 <- d5
  result$d_abx <- d_abx[, c(1:5)]
  result
}
