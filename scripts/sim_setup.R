library(dplyr)
library(tidyr)
library(mpMap)
library(qtl)
library(mpwgaim)
library(DiGGer)
library(ggplot2)
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
cr_int <- mpcross2int(cr, gen.type = "interval")
saveRDS(cr_int, "data/mp4_int.rds")
effects_df <- effects_mat %>% 
  as.data.frame() %>% 
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
  as.data.frame() %>% 
  mutate(fid = 1:900) %>% 
  gather(qtlid, founder, -fid) %>% 
  mutate(qtlid = sub("ibd1_", "", qtlid))

gen_effects <- effects_df %>% 
  full_join(qtl_geno) %>% 
  mutate(allele_cont = feffect * gen_var_prop)

gen_effects %>% 
  group_by(fid) %>% 
  summarise(all_sum = sum(allele_cont)) %>% 
  mutate(id = as.factor(paste0("L", fid))) -> phen_gen

saveRDS(phen_gen, "data/phen_gen.rds")
# cr.dat <- mpprob(cr, 3)
# crq <- mpIM(object = cr.dat, ncov = 0, responsename = "pheno")

#Make phase1 design:
#Field. 20 rows, 70 ranges
#900 RILs, 500 reps

ph1.des <- des.prep00(900, 20, 70, rep(c(2, 1), c(500, 400)), phen_gen$id,
                      ribs = 20, cibs = 35, tgrp = rep(c(2, 1), c(500, 400)))
ph1_out <- run(ph1.des)
ph1_des_out <- ph1_out$dlist %>% 
  as.tbl() %>% 
  rename(row = ROW, range = RANGE, id = ID)
ph1_row_effects <- data_frame(
  row = 1:20,
  ph1_row_effect = seq(from = -0.25, to = 0.25, length = 20)
)
ph1_range_effects <- data_frame(
  range = 1:70,
  ph1_range_effect = seq(from = -0.4, to = 0.4, length = 70)
)
ph1_pheno <- ph1_des_out %>% 
  left_join(phen_gen) %>% 
  left_join(ph1_row_effects) %>% 
  left_join(ph1_range_effects) %>% 
  mutate(residual = rnorm(1400, sd = 0.1),
         ph1_pheno = all_sum + ph1_row_effect + ph1_range_effect + residual)

saveRDS(ph1_pheno, file = "data/ph1_pheno.rds")
ggplot(ph1_des, aes(row, range)) + geom_tile(aes(fill = ph1_pheno))
