n1 <- 20
n2 <- 10
overlap <- n2/2
S = 100
burn = S * .1

Z_true <- rep(0, n2)
Z_true[1:overlap] <- 1:overlap

show_progress <- T
fast = F
R <- NULL
all_patterns <- T

m <- rep(c(.95, .05), 5)
u <- c(.01, .99, .01, .99,
       1/30, 1- 1/30, 1/12, 1 - 1/12, 1/15, 1 - 1/15)

levels <- rep(2, 5)
cd <- simulate_comparisons(m, u, levels, n1, n2, overlap)
cd_mm <- simulate_comparisons_mm(m, u, levels, n1, n2, c(overlap, 1))
cd_new <- swap_base_file(cd_mm)
hash <- hash_comparisons(cd_new)
out <- fabl(hash, S = S, burn = burn)
result <- estimate_links(out, hash)
result$Z_hat


hash_mm <- hash_comparisons(cd_mm)
out <- fabl_mm(hash_mm, S = S, burn = burn)
result <- estimate_links_mm(out, hash_mm)
result$Z_hat
