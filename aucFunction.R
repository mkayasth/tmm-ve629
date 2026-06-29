run_harmonic_cv_selection_singscore <- function(expr_total, meta_total,
                                                expr_test_list = NULL,   # named list of expression matrices (optional)
                                                meta_test_list = NULL,   # named list of metadata data frames (optional)
                                                candidate_genes_up,
                                                candidate_genes_down,
                                                phenotype_col, batch_col,
                                                label_neg,
                                                min_per_subgroup,
                                                n_folds, n_repeats,
                                                n_cores = 1,
                                                max_genes = 20,
                                                pivot_genes = NULL,
                                                n_pivots    = 3L,
                                                lcb_conf   = 0.95,
                                                lcb_boot_R = 500,
                                                perm_R     = 500,
                                                lcb_tol    = 1e-3) {
  
  require(pROC)
  require(caret)
  require(singscore)
  require(GSEABase)
  
  # ── Input coercion: accept bare matrix/df, a list, or NULL (test sets optional)
  if (is.null(expr_test_list) || (is.list(expr_test_list) && !is.data.frame(expr_test_list) && length(expr_test_list) == 0)) {
    expr_test_list <- list()
  } else if (!is.list(expr_test_list) || is.data.frame(expr_test_list) || is.matrix(expr_test_list)) {
    expr_test_list <- list(TestSet1 = expr_test_list)
  }
  
  if (is.null(meta_test_list) || (is.list(meta_test_list) && !is.data.frame(meta_test_list) && length(meta_test_list) == 0)) {
    meta_test_list <- list()
  } else if (!is.list(meta_test_list) || is.data.frame(meta_test_list)) {
    meta_test_list <- list(TestSet1 = meta_test_list)
  }
  
  if (length(expr_test_list) != length(meta_test_list))
    stop("[VALIDATION] expr_test_list and meta_test_list must have the same length.")
  
  if (length(expr_test_list) > 0) {
    if (is.null(names(expr_test_list)))
      names(expr_test_list) <- paste0("TestSet", seq_along(expr_test_list))
    if (is.null(names(meta_test_list)))
      names(meta_test_list) <- names(expr_test_list)
  }
  
  has_test_sets  <- length(expr_test_list) > 0
  test_set_names <- if (has_test_sets) names(expr_test_list) else character(0)
  
  # ── Safe subset helper --------------------------------------------------------
  safe_subset <- function(df, idx) {
    out <- df[idx, , drop = FALSE]
    row.names(out) <- NULL
    out
  }
  
  set.seed(123)
  start_total <- Sys.time()
  message("\n[SYSTEM] --- Starting Robust Hybrid Selection ",
          "(Singscore + Per-Group Bootstrap LCB + Test Eval + Permutation P) ---")
  message(sprintf("[SYSTEM] Phenotype column : '%s'", phenotype_col))
  message(sprintf("[SYSTEM] Negative label   : '%s'", label_neg))
  if (has_test_sets) {
    message(sprintf("[SYSTEM] Test sets        : %s",
                    paste(test_set_names, collapse = ", ")))
  } else {
    message("[SYSTEM] Test sets        : (none — train-only mode)")
  }
  
  # ── Sanitise training metadata ------------------------------------------------
  meta_total <- as.data.frame(meta_total, stringsAsFactors = FALSE)
  row.names(meta_total) <- NULL
  meta_total <- safe_subset(meta_total, !duplicated(meta_total$SampleID))
  
  expr_total <- expr_total[, colnames(expr_total) %in% meta_total$SampleID, drop = FALSE]
  expr_total <- expr_total[, match(meta_total$SampleID, colnames(expr_total)), drop = FALSE]
  
  # ── Sanitise test set list (only when test sets were provided) ---------------
  if (has_test_sets) {
    meta_test_list <- lapply(meta_test_list, function(m) {
      m <- as.data.frame(m, stringsAsFactors = FALSE)
      row.names(m) <- NULL
      safe_subset(m, !duplicated(m$SampleID))
    })
    
    expr_test_list <- mapply(function(expr, meta) {
      expr <- expr[, colnames(expr) %in% meta$SampleID, drop = FALSE]
      expr[, match(meta$SampleID, colnames(expr)), drop = FALSE]
    }, expr_test_list, meta_test_list, SIMPLIFY = FALSE)
  }
  
  # ── Derive comparisons dynamically from phenotype column ---------------------
  all_pheno_levels <- unique(meta_total[[phenotype_col]])
  pos_levels       <- setdiff(all_pheno_levels, label_neg)
  
  if (length(pos_levels) == 0)
    stop(sprintf(
      "[VALIDATION] No positive phenotype levels found after excluding label_neg='%s'.\n",
      "  Levels present: %s", label_neg, paste(all_pheno_levels, collapse = ", ")
    ))
  
  # Named list: short_key -> phenotype_value
  # Short keys are used in column names / messages (spaces replaced with "_")
  make_key <- function(x) gsub("[^A-Za-z0-9_]", "_", x)
  comparisons     <- setNames(as.list(pos_levels), make_key(pos_levels))
  comparison_keys <- names(comparisons)
  
  message(sprintf("[SYSTEM] Positive groups (%d): %s",
                  length(pos_levels), paste(pos_levels, collapse = ", ")))
  
  # ── Helpers ------------------------------------------------------------------
  singscore_score <- function(up_genes, down_genes, expr) {
    ranked     <- rankGenes(expr)
    up_genes   <- up_genes[up_genes %in% rownames(expr)]
    down_genes <- down_genes[down_genes %in% rownames(expr)]
    
    if (length(up_genes) > 0 && length(down_genes) > 0) {
      scores <- simpleScore(rankData    = ranked,
                            upSet       = GSEABase::GeneSet(up_genes),
                            downSet     = GSEABase::GeneSet(down_genes),
                            centerScore = TRUE)
    } else if (length(up_genes) > 0) {
      scores <- simpleScore(rankData    = ranked,
                            upSet       = GSEABase::GeneSet(up_genes),
                            centerScore = TRUE)
    } else if (length(down_genes) > 0) {
      scores <- simpleScore(rankData    = ranked,
                            upSet       = GSEABase::GeneSet(down_genes),
                            centerScore = TRUE)
      scores$TotalScore <- -scores$TotalScore
    } else {
      stop("[singscore_score] Both up_genes and down_genes are empty.")
    }
    setNames(scores$TotalScore, rownames(scores))
  }
  
  harmonic_mean <- function(x) {
    x <- x[!is.na(x) & is.finite(x)]
    if (length(x) == 0 || any(x <= 0)) return(0)
    length(x) / sum(1 / x)
  }
  
  bootstrap_lcb <- function(fold_aucs, conf = 0.95, R = 500, seed = 42) {
    set.seed(seed)
    n          <- length(fold_aucs)
    boot_means <- replicate(R, mean(sample(fold_aucs, n, replace = TRUE)))
    quantile(boot_means, probs = 1 - conf, names = FALSE)
  }
  
  # ── Core CV worker: returns named vector of per-group LCBs ------------------
  run_cv <- function(trial_up, trial_down) {
    vapply(comparison_keys, function(grp_key) {
      target_grp  <- comparisons[[grp_key]]
      samples_idx <- which(meta_total[[phenotype_col]] %in% c(label_neg, target_grp))
      
      m_sub        <- safe_subset(meta_total, samples_idx)
      m_sub$target <- ifelse(m_sub[[phenotype_col]] == label_neg, 0L, 1L)
      
      strata     <- interaction(m_sub$target, m_sub[[batch_col]], drop = TRUE)
      strata_tbl <- split(seq_along(strata), strata)
      min_req    <- min_per_subgroup * n_folds
      
      os_local <- unlist(lapply(strata_tbl, function(s) {
        if (length(s) < min_req) {
          set.seed(100 + length(s))
          c(s, sample(s, min_req - length(s), replace = TRUE))
        } else s
      }), use.names = FALSE)
      
      m_sub_os      <- safe_subset(m_sub, os_local)
      global_col_os <- samples_idx[os_local]
      strata_os     <- interaction(m_sub_os$target, m_sub_os[[batch_col]], drop = TRUE)
      
      fold_aucs <- unlist(lapply(seq_len(n_repeats), function(r) {
        folds <- caret::createFolds(strata_os, k = n_folds, returnTrain = FALSE)
        vapply(folds, function(test_local) {
          if (length(unique(m_sub_os$target[test_local])) < 2) return(NA_real_)
          unique_col_idx    <- unique(global_col_os[c(setdiff(seq_along(strata_os),
                                                              test_local), test_local)])
          expr_fold         <- expr_total[, unique_col_idx, drop = FALSE]
          all_scores        <- singscore_score(trial_up, trial_down, expr_fold)
          test_sample_names <- colnames(expr_total)[global_col_os[test_local]]
          test_scores_fold  <- all_scores[test_sample_names]
          test_labels       <- m_sub_os$target[test_local]
          if (length(unique(test_labels)) < 2) return(NA_real_)
          roc_obj <- pROC::roc(test_labels, test_scores_fold, quiet = TRUE)
          auc_raw <- as.numeric(pROC::auc(roc_obj))
          max(auc_raw, 1 - auc_raw)
        }, FUN.VALUE = numeric(1))
      }), use.names = FALSE)
      
      fold_aucs <- fold_aucs[!is.na(fold_aucs)]
      if (length(fold_aucs) == 0) return(0)
      grp_seed <- 42 + which(comparison_keys == grp_key)
      bootstrap_lcb(fold_aucs, conf = lcb_conf, R = lcb_boot_R, seed = grp_seed)
    }, FUN.VALUE = numeric(1))
  }
  
  # ── Validate gene pools against training AND all test sets -------------------
  validate_pool <- function(genes, label) {
    miss_tr <- setdiff(genes, rownames(expr_total))
    if (length(miss_tr) > 0)
      warning(sprintf("[VALIDATION] %d %s gene(s) missing from expr_total: %s",
                      length(miss_tr), label, paste(head(miss_tr, 5), collapse = ", ")))
    
    keep <- intersect(genes, rownames(expr_total))
    if (has_test_sets) {
      for (ts in test_set_names) {
        miss_te <- setdiff(keep, rownames(expr_test_list[[ts]]))
        if (length(miss_te) > 0)
          warning(sprintf("[VALIDATION] %d %s gene(s) missing from expr_test '%s': %s",
                          length(miss_te), label, ts,
                          paste(head(miss_te, 5), collapse = ", ")))
        keep <- intersect(keep, rownames(expr_test_list[[ts]]))
      }
    }
    keep
  }
  
  pool_up   <- validate_pool(candidate_genes_up,   "UP")
  pool_down <- validate_pool(candidate_genes_down, "DOWN")
  
  # ── Parse pivot_genes (named char vector: c(GENE1="UP", GENE2="DOWN")) -------
  n_pivots <- as.integer(n_pivots)
  stopifnot(n_pivots >= 1L, n_pivots <= 10L)
  
  parse_pivot_genes <- function(pivot_genes, pool_up, pool_down, n_pivots) {
    if (is.null(pivot_genes) || length(pivot_genes) == 0)
      return(list(up = character(0), down = character(0),
                  pool_up = pool_up, pool_down = pool_down, n_explicit = 0L))
    
    if (!is.character(pivot_genes) || is.null(names(pivot_genes)))
      stop(paste0("[PIVOT] pivot_genes must be a named character vector.\n",
                  "  Example: c(GENE1 = 'UP', GENE2 = 'DOWN')"))
    
    directions <- toupper(unname(pivot_genes))
    genes      <- names(pivot_genes)
    
    bad_dir <- !directions %in% c("UP", "DOWN")
    if (any(bad_dir))
      stop(sprintf("[PIVOT] Invalid direction(s) for: %s. Use 'UP' or 'DOWN'.",
                   paste(genes[bad_dir], collapse = ", ")))
    
    if (length(genes) > n_pivots) {
      warning(sprintf("[PIVOT] %d pivot gene(s) supplied but n_pivots=%d; using first %d.",
                      length(genes), n_pivots, n_pivots))
      genes      <- genes[seq_len(n_pivots)]
      directions <- directions[seq_len(n_pivots)]
    }
    
    all_expr_genes <- if (has_test_sets) {
      Reduce(intersect, c(list(rownames(expr_total)), lapply(expr_test_list, rownames)))
    } else {
      rownames(expr_total)
    }
    missing <- !genes %in% all_expr_genes
    if (any(missing))
      warning(sprintf("[PIVOT] Pivot gene(s) absent from one or more expression matrices, skipped: %s",
                      paste(genes[missing], collapse = ", ")))
    genes      <- genes[!missing]
    directions <- directions[!missing]
    
    if (length(genes) == 0) {
      warning("[PIVOT] No valid pivot genes remain. All pivots will be auto-selected.")
      return(list(up = character(0), down = character(0),
                  pool_up = pool_up, pool_down = pool_down, n_explicit = 0L))
    }
    
    up_genes   <- genes[directions == "UP"]
    down_genes <- genes[directions == "DOWN"]
    
    not_in_up   <- setdiff(up_genes,   pool_up)
    not_in_down <- setdiff(down_genes, pool_down)
    if (length(not_in_up)   > 0)
      warning(sprintf("[PIVOT] UP pivot gene(s) not in candidate_genes_up pool: %s",
                      paste(not_in_up, collapse = ", ")))
    if (length(not_in_down) > 0)
      warning(sprintf("[PIVOT] DOWN pivot gene(s) not in candidate_genes_down pool: %s",
                      paste(not_in_down, collapse = ", ")))
    
    message(sprintf("[PIVOT] Locking %d user-supplied pivot gene(s):", length(genes)))
    for (k in seq_along(genes))
      message(sprintf("         %-20s  ->  %s", genes[k], directions[k]))
    
    list(up         = up_genes,
         down       = down_genes,
         pool_up    = setdiff(pool_up,   up_genes),
         pool_down  = setdiff(pool_down, down_genes),
         n_explicit = length(genes))
  }
  
  piv <- parse_pivot_genes(pivot_genes, pool_up, pool_down, n_pivots)
  
  selected_up       <- piv$up
  selected_down     <- piv$down
  pool_up           <- piv$pool_up
  pool_down         <- piv$pool_down
  n_explicit_pivots <- piv$n_explicit
  
  # ── Auto-select remaining pivots via greedy boot-LCB ------------------------
  n_auto_pivots <- n_pivots - n_explicit_pivots
  
  if (n_auto_pivots > 0) {
    message(sprintf("\n[PIVOT PHASE] Auto-selecting %d pivot gene(s) via greedy boot-LCB ...",
                    n_auto_pivots))
    
    for (pv in seq_len(n_auto_pivots)) {
      cand_list_pv <- list()
      if (length(pool_up)   > 0)
        cand_list_pv[["up"]]   <- data.frame(gene = pool_up,   direction = "UP",
                                             stringsAsFactors = FALSE)
      if (length(pool_down) > 0)
        cand_list_pv[["down"]] <- data.frame(gene = pool_down, direction = "DOWN",
                                             stringsAsFactors = FALSE)
      
      if (length(cand_list_pv) == 0) { message("[PIVOT] No candidates left."); break }
      pv_candidates <- do.call(rbind, cand_list_pv)
      
      pv_results <- lapply(seq_len(nrow(pv_candidates)), function(ci) {
        cand_gene  <- pv_candidates$gene[ci]
        cand_dir   <- pv_candidates$direction[ci]
        trial_up   <- if (cand_dir == "UP")   c(selected_up,   cand_gene) else selected_up
        trial_down <- if (cand_dir == "DOWN") c(selected_down, cand_gene) else selected_down
        if (length(trial_up) + length(trial_down) == 0) return(NULL)
        lcbs <- run_cv(trial_up, trial_down)
        list(gene = cand_gene, direction = cand_dir, hm_lcb = harmonic_mean(lcbs))
      })
      pv_results <- Filter(Negate(is.null), pv_results)
      if (length(pv_results) == 0) break
      
      hm_vals <- vapply(pv_results, `[[`, "hm_lcb", FUN.VALUE = numeric(1))
      best_pv <- pv_results[[which.max(hm_vals)]]
      
      if (best_pv$direction == "UP") {
        selected_up <- c(selected_up, best_pv$gene)
        pool_up     <- setdiff(pool_up, best_pv$gene)
      } else {
        selected_down <- c(selected_down, best_pv$gene)
        pool_down     <- setdiff(pool_down, best_pv$gene)
      }
      
      message(sprintf("[PIVOT %d/%d] Locked: %-20s  ->  %s  |  HM-LCB: %.4f",
                      n_explicit_pivots + pv, n_pivots,
                      best_pv$gene, best_pv$direction, best_pv$hm_lcb))
    }
  }
  
  message(sprintf(
    "\n[PIVOT] Final pivot set (%d gene(s)):\n  UP  : %s\n  DOWN: %s",
    length(selected_up) + length(selected_down),
    if (length(selected_up)   > 0) paste(selected_up,   collapse = ", ") else "(none)",
    if (length(selected_down) > 0) paste(selected_down, collapse = ", ") else "(none)"
  ))
  
  # ── Build history data frame with dynamic column names ----------------------
  # Per-comparison-group train columns: lcb_<key>, sd_<key>
  # Per-comparison-group test columns (per test set): <ts>_auc_<key>, <ts>_sd_<key>, <ts>_pval_<key>
  train_lcb_cols <- paste0("lcb_", comparison_keys)
  train_sd_cols  <- paste0("sd_",  comparison_keys)
  
  if (has_test_sets) {
    test_auc_cols  <- as.vector(outer(test_set_names, comparison_keys,
                                      function(ts, ck) paste0(ts, "_auc_",  ck)))
    test_sd_cols_h <- as.vector(outer(test_set_names, comparison_keys,
                                      function(ts, ck) paste0(ts, "_sd_",   ck)))
    test_pval_cols <- as.vector(outer(test_set_names, comparison_keys,
                                      function(ts, ck) paste0(ts, "_pval_", ck)))
    test_hm_cols   <- paste0(test_set_names, "_hm_auc")
    test_hmp_cols  <- paste0(test_set_names, "_pval_hm_auc")
  } else {
    test_auc_cols <- test_sd_cols_h <- test_pval_cols <-
      test_hm_cols <- test_hmp_cols <- character(0)
  }
  
  fixed_cols <- c("size", "gene_added", "direction", "hm_lcb", "sd_pooled", "method")
  all_cols   <- c(fixed_cols,
                  train_lcb_cols, train_sd_cols,
                  test_hm_cols,   test_auc_cols,
                  test_sd_cols_h, test_pval_cols, test_hmp_cols)
  
  history           <- as.data.frame(
    matrix(NA_real_, nrow = max_genes, ncol = length(all_cols)),
    stringsAsFactors = FALSE
  )
  colnames(history) <- all_cols
  history$gene_added <- NA_character_
  history$direction  <- NA_character_
  history$method     <- NA_character_
  history$size       <- NA_integer_
  
  # Seed history rows for pivot genes
  all_pivots   <- c(selected_up, selected_down)
  all_piv_dirs <- c(rep("UP", length(selected_up)), rep("DOWN", length(selected_down)))
  for (k in seq_along(all_pivots)) {
    history[k, "size"]       <- k
    history[k, "gene_added"] <- all_pivots[k]
    history[k, "direction"]  <- all_piv_dirs[k]
    history[k, "method"]     <- if (k <= n_explicit_pivots) "USER-PIVOT" else "AUTO-PIVOT"
  }
  
  # ── Test-set evaluation helper (operates over all test sets) ----------------
  evaluate_test_set <- function(up_genes, down_genes, boot_R = 200,
                                perm_R = 500, seed = 99) {
    results <- list()
    
    for (ts in test_set_names) {
      expr_te <- expr_test_list[[ts]]
      meta_te <- meta_test_list[[ts]]
      
      te_scores    <- singscore_score(up_genes, down_genes, expr_te)
      grp_auc      <- setNames(rep(NA_real_, length(comparison_keys)), comparison_keys)
      grp_sd       <- grp_auc
      grp_pval     <- grp_auc
      
      for (gi in seq_along(comparison_keys)) {
        grp_key    <- comparison_keys[gi]
        target_grp <- comparisons[[grp_key]]
        grp_rows   <- which(meta_te[[phenotype_col]] %in% c(label_neg, target_grp))
        
        if (length(grp_rows) < 2 ||
            length(unique(meta_te[[phenotype_col]][grp_rows])) < 2) next
        
        grp_ids    <- meta_te$SampleID[grp_rows]
        grp_scores <- te_scores[grp_ids]
        grp_labels <- ifelse(meta_te[[phenotype_col]][grp_rows] == label_neg, 0L, 1L)
        
        roc_obs <- pROC::roc(grp_labels, grp_scores, quiet = TRUE)
        auc_obs <- as.numeric(pROC::auc(roc_obs))
        auc_obs <- max(auc_obs, 1 - auc_obs)
        
        grp_seed <- seed + gi
        n_te     <- length(grp_rows)
        
        set.seed(grp_seed)
        boot_aucs <- replicate(boot_R, {
          b_idx <- sample(n_te, n_te, replace = TRUE)
          if (length(unique(grp_labels[b_idx])) < 2) return(NA_real_)
          r <- pROC::roc(grp_labels[b_idx], grp_scores[b_idx], quiet = TRUE)
          b <- as.numeric(pROC::auc(r)); max(b, 1 - b)
        })
        boot_aucs <- boot_aucs[!is.na(boot_aucs)]
        
        set.seed(grp_seed + 1000)
        perm_aucs <- replicate(perm_R, {
          perm_labels <- sample(grp_labels)
          if (length(unique(perm_labels)) < 2) return(NA_real_)
          r <- pROC::roc(perm_labels, grp_scores, quiet = TRUE)
          b <- as.numeric(pROC::auc(r)); max(b, 1 - b)
        })
        perm_aucs <- perm_aucs[!is.na(perm_aucs)]
        
        grp_auc[gi]  <- auc_obs
        grp_sd[gi]   <- if (length(boot_aucs) > 1) sd(boot_aucs) else NA_real_
        grp_pval[gi] <- if (length(perm_aucs) > 0) mean(perm_aucs >= auc_obs) else NA_real_
      }
      
      results[[ts]] <- list(
        hm_auc      = harmonic_mean(grp_auc),
        auc         = grp_auc,
        sd          = grp_sd,
        pval        = grp_pval,
        pval_hm_auc = if (all(is.na(grp_pval))) NA_real_
        else min(grp_pval, na.rm = TRUE)
      )
    }
    results
  }
  
  # ── Greedy forward selection loop -------------------------------------------
  start_i <- length(selected_up) + length(selected_down) + 1L
  
  if (start_i > max_genes)
    message("[WARN] n_pivots >= max_genes. No greedy expansion will occur.")
  
  for (i in start_i:max_genes) {
    message(sprintf("\n[GREEDY] Step %d / %d", i, max_genes))
    
    cand_list <- list()
    if (length(pool_up)   > 0)
      cand_list[["up"]]   <- data.frame(gene = pool_up,   direction = "UP",
                                        stringsAsFactors = FALSE)
    if (length(pool_down) > 0)
      cand_list[["down"]] <- data.frame(gene = pool_down, direction = "DOWN",
                                        stringsAsFactors = FALSE)
    
    if (length(cand_list) == 0) { message("[WARN] No candidates left."); break }
    candidates <- do.call(rbind, cand_list)
    
    iter_results <- lapply(seq_len(nrow(candidates)), function(ci) {
      cand_gene  <- candidates$gene[ci]
      cand_dir   <- candidates$direction[ci]
      trial_up   <- if (cand_dir == "UP")   c(selected_up,   cand_gene) else selected_up
      trial_down <- if (cand_dir == "DOWN") c(selected_down, cand_gene) else selected_down
      if (length(trial_up) + length(trial_down) == 0) return(NULL)
      
      lcbs <- run_cv(trial_up, trial_down)
      
      # Also collect per-group SD over folds for history
      per_grp_sd    <- setNames(rep(NA_real_, length(comparison_keys)), comparison_keys)
      all_fold_aucs <- numeric(0)
      
      for (gi in seq_along(comparison_keys)) {
        grp_key    <- comparison_keys[gi]
        target_grp <- comparisons[[grp_key]]
        samples_idx <- which(meta_total[[phenotype_col]] %in% c(label_neg, target_grp))
        
        m_sub        <- safe_subset(meta_total, samples_idx)
        m_sub$target <- ifelse(m_sub[[phenotype_col]] == label_neg, 0L, 1L)
        strata       <- interaction(m_sub$target, m_sub[[batch_col]], drop = TRUE)
        strata_tbl   <- split(seq_along(strata), strata)
        min_req      <- min_per_subgroup * n_folds
        
        os_local <- unlist(lapply(strata_tbl, function(s) {
          if (length(s) < min_req) {
            set.seed(100 + length(s))
            c(s, sample(s, min_req - length(s), replace = TRUE))
          } else s
        }), use.names = FALSE)
        
        m_sub_os      <- safe_subset(m_sub, os_local)
        global_col_os <- samples_idx[os_local]
        strata_os     <- interaction(m_sub_os$target, m_sub_os[[batch_col]], drop = TRUE)
        
        fold_aucs <- unlist(lapply(seq_len(n_repeats), function(r) {
          folds <- caret::createFolds(strata_os, k = n_folds, returnTrain = FALSE)
          vapply(folds, function(test_local) {
            if (length(unique(m_sub_os$target[test_local])) < 2) return(NA_real_)
            unique_col_idx    <- unique(global_col_os[c(setdiff(seq_along(strata_os),
                                                                test_local), test_local)])
            expr_fold         <- expr_total[, unique_col_idx, drop = FALSE]
            all_scores        <- singscore_score(trial_up, trial_down, expr_fold)
            test_sample_names <- colnames(expr_total)[global_col_os[test_local]]
            test_scores_fold  <- all_scores[test_sample_names]
            test_labels       <- m_sub_os$target[test_local]
            if (length(unique(test_labels)) < 2) return(NA_real_)
            roc_obj <- pROC::roc(test_labels, test_scores_fold, quiet = TRUE)
            auc_raw <- as.numeric(pROC::auc(roc_obj))
            max(auc_raw, 1 - auc_raw)
          }, FUN.VALUE = numeric(1))
        }), use.names = FALSE)
        
        fold_aucs           <- fold_aucs[!is.na(fold_aucs)]
        per_grp_sd[gi]      <- if (length(fold_aucs) > 1) sd(fold_aucs) else NA_real_
        all_fold_aucs       <- c(all_fold_aucs, fold_aucs)
      }
      
      list(gene      = cand_gene,
           direction = cand_dir,
           hm_lcb    = harmonic_mean(lcbs),
           sd_pooled = if (length(all_fold_aucs) > 1) sd(all_fold_aucs) else NA_real_,
           lcbs      = lcbs,           # named vector, length == n comparisons
           grp_sds   = per_grp_sd)     # named vector, length == n comparisons
    })
    
    iter_results <- Filter(Negate(is.null), iter_results)
    if (length(iter_results) == 0) { message("[WARN] All NULL."); break }
    
    res_df <- do.call(rbind, lapply(iter_results, function(x) {
      base <- data.frame(gene = x$gene, direction = x$direction,
                         hm_lcb = x$hm_lcb, sd_pooled = x$sd_pooled,
                         stringsAsFactors = FALSE, row.names = NULL)
      for (ck in comparison_keys) {
        base[[paste0("lcb_", ck)]] <- x$lcbs[ck]
        base[[paste0("sd_",  ck)]] <- x$grp_sds[ck]
      }
      base
    }))
    
    # ── Dual-criterion selection: LCB-first, SD-minimisation on saturation ------
    max_lcb     <- max(res_df$hm_lcb, na.rm = TRUE)
    saturated   <- (max_lcb - res_df$hm_lcb) <= lcb_tol   # candidates within tol of best
    
    if (sum(saturated, na.rm = TRUE) > 1L) {
      # All saturated candidates have indistinguishably high LCB:
      # break the tie by minimising pooled SD across all folds and groups.
      sat_df    <- res_df[saturated, , drop = FALSE]
      best_idx  <- which(saturated)[which.min(sat_df$sd_pooled)]
      selection_criterion <- sprintf("SD-MIN (%.0f candidates saturated within lcb_tol=%.4f)",
                                     sum(saturated), lcb_tol)
    } else {
      best_idx  <- which.max(res_df$hm_lcb)
      selection_criterion <- "HM-LCB-MAX"
    }
    
    best_gene <- res_df$gene[best_idx]
    best_dir  <- res_df$direction[best_idx]
    
    if (best_dir == "UP") {
      selected_up <- c(selected_up, best_gene)
      pool_up     <- setdiff(pool_up, best_gene)
    } else {
      selected_down <- c(selected_down, best_gene)
      pool_down     <- setdiff(pool_down, best_gene)
    }
    
    # Write fixed train columns
    history[i, "size"]       <- i
    history[i, "gene_added"] <- best_gene
    history[i, "direction"]  <- best_dir
    history[i, "hm_lcb"]     <- res_df$hm_lcb[best_idx]
    history[i, "sd_pooled"]  <- res_df$sd_pooled[best_idx]
    history[i, "method"]     <- paste0("SINGSCORE-BOOTLCB|", selection_criterion)
    
    for (ck in comparison_keys) {
      history[i, paste0("lcb_", ck)] <- res_df[[paste0("lcb_", ck)]][best_idx]
      history[i, paste0("sd_",  ck)] <- res_df[[paste0("sd_",  ck)]][best_idx]
    }
    
    # ── Test evaluation (only when test sets were provided) --------------------
    if (has_test_sets) {
      message(sprintf(" >>> [SIZE %d] Best: %s (%s) | Evaluating on test set(s)...",
                      i, best_gene, best_dir))
      
      test_eval <- evaluate_test_set(selected_up, selected_down,
                                     boot_R = 200, perm_R = perm_R, seed = 99 + i)
      
      for (ts in test_set_names) {
        te <- test_eval[[ts]]
        history[i, paste0(ts, "_hm_auc")]      <- te$hm_auc
        history[i, paste0(ts, "_pval_hm_auc")] <- te$pval_hm_auc
        for (ck in comparison_keys) {
          history[i, paste0(ts, "_auc_",  ck)] <- te$auc[ck]
          history[i, paste0(ts, "_sd_",   ck)] <- te$sd[ck]
          history[i, paste0(ts, "_pval_", ck)] <- te$pval[ck]
        }
      }
    }
    
    # ── Pretty print progress --------------------------------------------------
    train_lines <- paste(vapply(comparison_keys, function(ck) {
      sprintf("             %-20s -> LCB: %0.4f | SD: %0.4f",
              ck,
              history[i, paste0("lcb_", ck)],
              history[i, paste0("sd_",  ck)])
    }, character(1)), collapse = "\n")
    
    if (has_test_sets) {
      test_lines <- paste(vapply(test_set_names, function(ts) {
        te_hdr <- sprintf("     TEST [%s]  | HM-AUC: %0.4f  (p_min: %0.4f)",
                          ts, test_eval[[ts]]$hm_auc, test_eval[[ts]]$pval_hm_auc)
        te_grp <- paste(vapply(comparison_keys, function(ck) {
          sprintf("             %-20s -> AUC: %0.4f | SD: %0.4f | p: %0.4f",
                  ck,
                  test_eval[[ts]]$auc[ck],
                  test_eval[[ts]]$sd[ck],
                  test_eval[[ts]]$pval[ck])
        }, character(1)), collapse = "\n")
        paste(te_hdr, te_grp, sep = "\n")
      }, character(1)), collapse = "\n")
    } else {
      test_lines <- "     TEST   | (no test sets provided)"
    }
    
    message(sprintf(
      paste0(" >>> [SIZE %d] Added: %s (%s)  [criterion: %s]\n",
             "     TRAIN  | HM-LCB: %0.4f  | SD(pooled): %0.4f\n",
             "%s\n",
             "%s"),
      i, best_gene, best_dir, selection_criterion,
      history[i, "hm_lcb"], history[i, "sd_pooled"],
      train_lines, test_lines
    ))
  }
  
  
  end_total <- Sys.time()
  message(sprintf("\n[SYSTEM] Done. Elapsed: %0.2f mins.",
                  as.numeric(difftime(end_total, start_total, units = "mins"))))
  
  # ── Return ------------------------------------------------------------------
  list(
    pivot_genes      = c(
      setNames(rep("UP",   length(selected_up)),   selected_up),
      setNames(rep("DOWN", length(selected_down)), selected_down)
    )[seq_len(min(length(selected_up) + length(selected_down), n_pivots))],
    final_genes_up   = selected_up,
    final_genes_down = selected_down,
    comparisons      = comparisons,   # the dynamic key -> phenotype_value map
    history          = history[!is.na(history$size), , drop = FALSE]
  )
}