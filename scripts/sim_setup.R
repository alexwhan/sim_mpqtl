library(dplyr)
library(tidyr)
library(mpMap)
lengths <- sample(250:350, 21)
nmarkers <- round(lengths * rnorm(21, 1.5))
map <- sim.map(lengths, nmarkers, include.x = FALSE)
expand.grid(1:7, c("A", "B", "D"), stringsAsFactors = FALSE) -> nm
names(map) <- paste0(nm$Var1, nm$Var2)
ped <- sim.mpped(4, 3, 15, 6, 20)

pick_loc <- function(map, names) {
  locs <- lapply(names, function(x) {
    chrom <- map[[x]]
    locs <- seq(from = 0, to = max(chrom), by = 0.1)
    loc_out <- sample(locs, 1)
    return(loc_out)
  })
  return(unlist(locs))
}
set.seed(123)
qtl_chroms <- sample(length(map), 30, TRUE)
qtl_locs <- pick_loc(map, qtl_chroms)
qtl_sizes <- rexp(30, rate = 0.1) + 3
qtl_sizes <- qtl_sizes / (sum(qtl_sizes))

make_effects <- function(nfounders) {
  f3 <- round(rnorm(nfounders - 1), 2)
  f4 <- sum(f3) - 2*sum(f3)
  f_out <- c(f3, f4)
  return(f_out)
}

set.seed(123)
lapply(rep(4, times = 30), make_effects) -> effects_list
effects_mat <- effects_list %>% unlist %>% matrix(30, 4, TRUE)

qtl_model <- cbind(qtl_chroms, qtl_locs, 
                   effects_mat)

cr <- sim.mpcross(map, ped, qtl_model, seed = 123)
effects_df <- effects_mat %>% 
  as_data_frame() %>% 
  setNames(paste0("f", 1:4)) %>% 
  mutate(qtlid = paste0("QTL", 1:30)) %>% 
  gather(founder, feffect, -qtlid) %>% 
  mutate(founder = as.numeric(sub("f", "", founder)))

set.seed(123)
gen_var_prop <- rexp(30, rate = 0.1) ^ 1.3
gen_var_prop[gen_var_prop < 3] <- gen_var_prop[gen_var_prop < 3] * 20
gen_var_prop[gen_var_prop < 8] <- gen_var_prop[gen_var_prop < 8] * 4
gen_var_prop <- round(gen_var_prop / sum(gen_var_prop), 3)
gen_var_df <- data_frame(
  qtlid = paste0("QTL", 1:30),
  gen_var_prop = gen_var_prop
)

effects_df <- effects_df %>% 
  left_join(gen_var_df)

qtl_geno <- cr$qtlgeno$finals[,1:30] %>% 
  as_data_frame() %>% 
  mutate(fid = 1:900) %>% 
  gather(qtlid, founder, -fid) %>% 
  mutate(qtlid = sub("ibd1_", "", qtlid))

gen_effects <- effects_df %>% 
  full_join(qtl_geno) %>% 
  mutate(allele_cont = feffect * gen_var_prop)

gen_effects %>% 
  group_by(fid) %>% 
  summarise(all_sum = sum(allele_cont)) -> temp
cr.dat <- mpprob(cr, 3)
crq <- mpIM(object = cr.dat, ncov = 0, responsename = "pheno")
