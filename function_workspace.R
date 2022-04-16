sample_info <- function(data) {
  # Create metadata tibble
  metadata <- tibble(colnames(get("se_soybean_cn_sub"))) %>%
    # Rename samples column
    rename(samples = 'colnames(get("se_soybean_cn_sub"))') %>%
    # Add column for replicates
    mutate(replicate = str_extract(samples, "(?<=\\.).*"))
    
}

