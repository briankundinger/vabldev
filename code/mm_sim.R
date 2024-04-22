library(vabldev)

#k <- SLURM
k <- 2
k_vec <- 1:k
n1 <- 2^6 # target
n2 <- 2^6 # base
overlap <- n2 / (2^k_vec)
S <- 200
burn <- S * .1


m <- rep(c(.95, .05), 5)
u <- c(.01, .99, .01, .99,
       1/30, 1- 1/30, 1/12, 1 - 1/12, 1/15, 1 - 1/15)

levels <- rep(2, 5)
cd <- simulate_comparisons_mm(m, u, levels, n1, n2, overlap)
hash <- hash_comparisons(cd, all_patterns = F)

out_fl <- variational_fastlink(hash)
estimate_fl <- estimate_links_fl(out_fl, hash)
Z_hat <- data.frame(id_1 = estimate_fl$fs_linkages$a,
                    id_2 = estimate_fl$fs_linkages$b)
eval_fl <- evaluate_links(Z_hat, cd$Z_true_long, n1, "pairs")

out_fabl <- fabl(hash, S = S, burn = burn)
estimate_fabl <- estimate_links(out_fabl, hash)
Z_hat <- make_Zhat_pairs(estimate_fabl$Z_hat)
eval_fabl <- evaluate_links(Z_hat, cd$Z_true_long, n1, "pairs")

out_mm <- fabl_mm(hash, S = S, burn = burn)
estimate_fabl <- estimate_links_mm(out_mm, hash)
Z_hat <- data.frame(id_1 = estimate_fabl$Z_hat$target_id,
                    id_2 = estimate_fabl$Z_hat$base_id)
eval_mm <- evaluate_links(Z_hat, cd$Z_true_long, n1, "pairs")
