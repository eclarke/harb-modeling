#' ARB Data Generation Functions
#' 
#' Functions should be named with genXXXX, where XXXX is the version number.
#' 
#' Document the parameters and changes each makes to the structure, and the 
#' model numbers the generators work with.
#' 
#' Generic functions should not be versioned.

#' Loads and cleans MIRA data. 
#' This function uses relative paths and may not function outside the 
#' original directory structure.
#' 
#' @return a phyloseq experiment object with the relevant data
load_data <- function() {
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
  prune_samples(
    sample_data(mira)$specimen_id2 %in% .abx_specimens, mira)
}

#' Converts the phyloseq object from `load_data` into an agglomerated dataframe
#' and cleans up subject/antibiotic data.
#' 
#' @return a cleaned and agglomerated dataframe
clean_data <- function(mira.abx, .genus, .species=NULL) {
  if (missing(.genus)) {
    stop("Must specify genus")
  }
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
    mutate(abx_idx = as.numeric(as.factor(abx_b)))
  
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

gen_01 <- function(d2, .specimen_type, min_times_abx=2, min_median_abx=2, min_reads_subj=1) {
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
gen_02 <- function(d2, .specimen_type, min_times_abx=2, min_median_abx=2, min_reads_subj=1) {
  d3 <- d2 %>% filter(specimen_type==.specimen_type) %>%
    distinct(subject_id, specimen_id2, read_count, total_reads, study_day) %>%
    group_by(subject_id) %>%
    filter(sum(read_count) > min_reads_subj) %>%
    arrange(subject_id, study_day) %>%
    mutate(prev = lag(read_count)/lag(total_reads)) %>%
    mutate(abx_b = lag(abx_b)) %>%
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
