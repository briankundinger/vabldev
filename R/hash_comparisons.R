#' @export
#'
hash_comparisons <- function(cd,
         algorithm = c("vabl", "fabl", "BRL_hash"), R = 0,
         all_patterns = FALSE, store_pair_to_pattern = TRUE){


  if("BRL_hash" %in% algorithm){
    R <- 0
  }
  indicators <- cd[[1]]
  N <- dim(indicators)[1]
  n1 <- cd[[2]]
  n2 <- cd[[3]]

  levels <- cd[[4]]
  fields <- seq_along(cd[[4]])
  field_marker <- sapply(fields, function(x){
    rep(x, cd[[4]][x])
  }) %>%
    unlist(.) %>%
    as.vector(.)

  ids <- expand.grid(1:n1, 1:n2)
  rec1 <- ids[,1]
  rec2 <- ids[,2]

  Lf_vec<- (levels) %>%
    c(0, .) %>%
    cumsum()

  hash_vals <- purrr::imap(cd[[4]], ~hash_field(.x, .y, Lf_vec)) %>%
    unlist()

  hash <- sweep(indicators, 2, hash_vals, "*") %>%
    rowSums() + 1

  if(all_patterns == TRUE){

    unique_patterns <- possible_patterns_ohe(levels)
    unique_hashed <- sweep(unique_patterns, 2, hash_vals, "*") %>%
      rowSums() + 1
    P <- dim(unique_patterns)[1]
    hash_id <- match(hash, unique_hashed) %>%
      factor(levels = 1:P)

  } else {

    unique_hashed <- unique(hash)
    P <- length(unique_hashed)
    hash_id <- match(hash, unique_hashed) %>%
      factor(levels = 1:P)
    unique_patterns <- indicators[!duplicated(hash_id), ]
  }

  temp <- data.frame(rec1, rec2, hash_id)

  hash_count_list <- temp %>%
    group_by(rec2, hash_id, .drop = F) %>%
    count() %>%
    ungroup() %>%
    group_split(rec2) %>%
    purrr::map(~.x %>%
          select(n) %>%
          pull()
        )

  total_counts <- temp %>%
    group_by(hash_id, .drop = F) %>%
    count() %>%
    pull()

  pattern_lookup <- expand.grid(1:P, 1:n2) %>%
    data.frame() %>%
    setNames(., c("hash_id", "rec2"))

  pair_to_pattern <- NULL

  if("BRL_hash" %in% algorithm){
  pair_to_pattern <- temp %>%
    select(hash_id, rec2) %>%
    group_split(rec2, .keep = F) %>%
    lapply(., pull) %>%
    lapply(., as.double)
  }

  hash_to_file_1 <- temp %>%
    select(rec1, rec2, hash_id) %>%
    nest_by(rec2, hash_id, .keep = F) %>%
    mutate(hash_id = as.integer(hash_id)) %>%
    rowwise() %>%
    mutate(N = nrow(data))

  hash_to_file_1 <- left_join(x = pattern_lookup,
                              y = hash_to_file_1,
                              by = c("hash_id", "rec2"))

  hash_to_file_1$N[is.na(hash_to_file_1$N)] <- 0

  flags <- NULL

  if("vabl" %in% algorithm){

  flags <- hash_to_file_1 %>%
    filter(N ==1) %>%
    tidyr::unnest(data) %>%
    tidyr::complete(rec2 = unique(hash_to_file_1$rec2)) %>%
    select(-N) %>%
    setNames(c("rec2", "eligible_patterns", "eligible_records")) %>%
    group_split(rec2, .keep = F)
  }


  # TODO: Try reworking code, conduct SEI with original hash_to_file1 format
  if("fabl" %in% algorithm | "BRL_hash" %in% algorithm){
    hash_to_file_1 <- hash_to_file_1 %>%
      group_split(rec2) %>%
      purrr::map(., ~ .x %>%
                   group_split(hash_id)) %>%
      purrr::map(., ~purrr::map(.x, `[[`, "data")) %>%
      purrr::map(., ~purrr::map(., ~ unname(unlist(.x))))

    if(R > 0){
      hash_to_file_1 <- lapply(hash_to_file_1, function(z){
        purrr::map(z, ~sei(.x, R))
      })}


  }

  # if("fabl" %in% algorithm | "BRL_hash" %in% algorithm){
  #   thing <- hash_to_file_1 %>%
  #     group_split(rec2)
  #   if(R > 0){
  #     lapply(thing, function(rec){
  #
  #     })
  #     thing <- lapply(hash_to_file_1$data, function(z){
  #       purrr::map(z, ~sei(.x, R))
  #     })
  #   }
  #
  #
  # }

  if(!("fabl" %in% algorithm) & !("BRL_hash" %in% algorithm)){
    hash_to_file_1 <- NULL
  }




  patterns <- list(ohe = unique_patterns,
                   total_counts = total_counts,
                   #pattern_counts_by_record = pattern_counts_by_record,
                   #record_counts_by_pattern = record_counts_by_pattern,
                   hash_count_list = hash_count_list,
                   hash_to_file_1 = hash_to_file_1,
                   flags = flags,
                   field_marker = field_marker,
                   n1 = n1,
                   n2 = n2,
                   pair_to_pattern = pair_to_pattern)
  patterns

}

