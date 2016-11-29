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
ped <- sim.mpped(4, 3, 160, 6, 2)


pick_loc <- function(map, n, min_dist) {
  locs <- data_frame(chrom = rep(NA_character_, times = n), 
                     pos = rep(NA_real_, times = n))
  for(i in seq_len(n)) {
    print(i)
    k <- 0
    while(k == 0) {
      new_chrom <- sample(names(map), 1)
      new_pos <- sample(seq(from = 0, to = max(map[[new_chrom]]), by = 0.1), 1)
      if(any(locs$chrom == new_chrom, na.rm = TRUE)) {
        old_pos <- as.numeric(locs$pos[locs$chrom == new_chrom])
        print(old_pos)
        print(new_pos)
        if(all(abs(old_pos - new_pos) > min_dist, na.rm = TRUE)) {
          locs[i,1] <- new_chrom
          locs[i, 2] <- new_pos
          k <- 1
        } 
      } else {
        locs[i,1] <- new_chrom
        locs[i, 2] <- new_pos
        k <- 1
      }
    }
  }
  return(locs)
}
set.seed(123)
qtl_locs <- pick_loc(map, 20, 50)
ggplot(qtl_locs, aes(chrom, pos)) + geom_point()
qtl_sizes <- rexp(20, rate = 0.1) + 3
qtl_sizes <- qtl_sizes / (sum(qtl_sizes))
plot(qtl_sizes)
make_effects <- function(nfounders) {
  f3 <- round(rnorm(nfounders - 1), 2)
  f4 <- sum(f3) - 2*sum(f3)
  f_out <- c(f3, f4)
  return(f_out)
}

set.seed(123)
lapply(rep(4, times = 20), make_effects) -> effects_list
effects_mat <- effects_list %>% unlist %>% matrix(20, 4, TRUE)

qtl_model <- cbind(qtl_locs, 
                   effects_mat)

cr <- sim.mpcross(map, ped, qtl_model, seed = 123)
qtl_model$qtlid <- paste0("QTL", 1:20)
cr_int <- mpcross2int(cr, gen.type = "mpInterval")
saveRDS(cr_int, "data/mp4_int.rds")
effects_df <- effects_mat %>% 
  as.data.frame() %>% 
  setNames(paste0("f", 1:4)) %>% 
  mutate(qtlid = paste0("QTL", 1:20)) %>% 
  gather(founder, feffect, -qtlid) %>% 
  mutate(founder = as.numeric(sub("f", "", founder)))

set.seed(123)
gen_var_prop <- rexp(20, rate = 0.1) ^ 1.3
gen_var_prop[gen_var_prop < 3] <- gen_var_prop[gen_var_prop < 3] * 20
gen_var_prop[gen_var_prop < 8] <- gen_var_prop[gen_var_prop < 8] * 4
gen_var_prop <- round(gen_var_prop / sum(gen_var_prop), 3)
gen_var_df <- data_frame(
  qtlid = paste0("QTL", 1:20),
  gen_var_prop = gen_var_prop
)

effects_df <- effects_df %>% 
  left_join(gen_var_df)

qtl_geno <- cr$qtlgeno$finals[,1:20] %>% 
  as.data.frame() %>% 
  mutate(fid = 1:960) %>% 
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

#Make phase1 design:
#Field. 20 rows, 70 ranges
#1200 RILs, 500 reps

ph1.des <- des.prep00(960, 30, 50, rep(c(2, 1), c(540, 420)), phen_gen$id,
                      ribs = 30, cibs = 25, tgrp = rep(c(2, 1), c(540, 420)))
ph1_outprep <- run(ph1.des)
ph1_des_outprep <- ph1_outprep$dlist %>%
  as.tbl() %>%
  rename(row = ROW, range = RANGE, id = ID)
ph1_row_effects <- data_frame(
  row = 1:30,
  ph1_row_effect = seq(from = -0.25, to = 0.25, length = 30)
)
  ph1_range_effects <- data_frame(
  range = 1:50,
  ph1_range_effect = seq(from = -0.4, to = 0.4, length = 50)
)
  ph1_row_effects2 <- data_frame(
    row = 1:30,
    ph1_row_effect = seq(from = -2, to = 2, length = 30)
  )
  ph1_range_effects2 <- data_frame(
    range = 1:50,
    ph1_range_effect = seq(from = -4, to = 4, length = 50)
  )
  
  ph1_phenoprep <- ph1_des_outprep %>%
  left_join(phen_gen) %>%
  left_join(ph1_row_effects) %>%
  left_join(ph1_range_effects) %>%
  mutate(residual = rnorm(1500, sd = 0.1),
         ph1_pheno = all_sum + ph1_row_effect + ph1_range_effect + residual)

saveRDS(ph1_pheno, file = "data/ph1_phenoprep.rds")

ph1_phenoprep2 <- ph1_des_outprep %>%
  left_join(phen_gen) %>%
  left_join(ph1_row_effects) %>%
  left_join(ph1_range_effects) %>%
  mutate(residual = rnorm(1500, sd = 0.5),
         ph1_pheno = all_sum + ph1_row_effect + ph1_range_effect + residual)

saveRDS(ph1_phenoprep2, file = "data/ph1_phenoprep2.rds")

# ggplot(ph1_pheno, aes(row, range)) + geom_tile(aes(fill = ph1_pheno))


#Make double rep design

ph1.desdub <- DiGGer(750, 30, 50, 30, 25,
                     RowsInBlockSequence = c(15),
                     ColumnsInBlockSequence = c(50),
                     TreatmentName = sample(phen_gen$id, 750, replace = FALSE),
                     TreatmentRepeatsPerReplicate = rep(1, 750))

treatment_fix <- data_frame(wrong_id = factor(1:750),
                            id = sample(phen_gen$id, 750, replace = FALSE))
ph1_outdub <- ph1.desdub$idsgn %>%
  as_data_frame() %>% 
  mutate(row = 1:30) %>% 
  gather(range, id, -row) %>% 
  mutate(range = as.numeric(sub("V", "", range)),
         id = factor(id),
         rep = ceiling(range / 25),
         block = ceiling(row / 15)) %>% 
  rename(wrong_id = id) %>% 
  left_join(treatment_fix)


ph1_outdub %>% 
  group_by(id) %>% 
  mutate(nrow = length(unique(row)),
            nrange = length(unique(range))) %>% 
  filter(nrow == 1 | nrange == 1) %>% 
  arrange(id)

ggplot(ph1_outdub, aes(row, range)) + geom_tile(colour = "black", fill = "transparent") +
  geom_line(aes(group = id))

ph1_phenodub <- ph1_outdub %>%
  left_join(phen_gen) %>%
  left_join(ph1_row_effects) %>%
  left_join(ph1_range_effects) %>%
  mutate(residual = rnorm(1500, sd = 0.1),
         ph1_pheno = all_sum + ph1_row_effect + ph1_range_effect + residual,
         id = factor(id))

saveRDS(ph1_phenodub, file = "data/ph1_phenodub.rds")

ph1_phenodub2 <- ph1_outdub %>%
  left_join(phen_gen) %>%
  left_join(ph1_row_effects) %>%
  left_join(ph1_range_effects) %>%
  mutate(residual = rnorm(1500, sd = 0.5),
         ph1_pheno = all_sum + ph1_row_effect + ph1_range_effect + residual,
         id = factor(id))

saveRDS(ph1_phenodub2, file = "data/ph1_phenodub2.rds")

#Make triple rep design
ph1.destrip <- DiGGer(500, 30, 50, 10, 50, 
                      RowsInBlockSequence = 30,
                      ColumnsInBlockSequence = 25,
                      TreatmentName = sample(phen_gen$id, 500, replace = FALSE))
treatment_fix2 <- data_frame(wrong_id = factor(1:500),
                            id = sample(phen_gen$id, 500, replace = FALSE))

ph1_outtrip <- ph1.destrip$idsgn %>% 
  as_data_frame() %>% 
  mutate(row = 1:30) %>% 
  gather(range, id, -row) %>% 
  mutate(range = as.numeric(sub("V", "", range)),
         id = factor(id),
         rep = ceiling(row / 10),
         block = ceiling(range / 25))%>% 
  rename(wrong_id = id) %>% 
  left_join(treatment_fix2)

ph1_outtrip %>% 
  group_by(id) %>% 
  mutate(nrow = length(unique(row)),
         nrange = length(unique(range))) %>% 
  filter(nrow == 1 | nrange == 1) %>% 
  arrange(id)

ggplot(ph1_outtrip, aes(row, range)) + geom_tile(colour = "black", fill = "transparent") +
  geom_line(aes(group = id))

ph1_phenotrip <- ph1_outtrip %>%
  left_join(phen_gen) %>%
  left_join(ph1_row_effects) %>%
  left_join(ph1_range_effects) %>%
  mutate(residual = rnorm(1500, sd = 0.1),
         ph1_pheno = all_sum + ph1_row_effect + ph1_range_effect + residual,
         id = factor(id))

saveRDS(ph1_phenotrip, file = "data/ph1_phenotrip.rds")

ph1_phenotrip2 <- ph1_outtrip %>%
  left_join(phen_gen) %>%
  left_join(ph1_row_effects) %>%
  left_join(ph1_range_effects) %>%
  mutate(residual = rnorm(1500, sd = 0.5),
         ph1_pheno = all_sum + ph1_row_effect + ph1_range_effect + residual,
         id = factor(id))

saveRDS(ph1_phenotrip2, file = "data/ph1_phenotrip2.rds")

##################
# Bigger spatial #
##################

ph1_row_effects2 <- data_frame(
  row = 1:30,
  ph1_row_effect = seq(from = -2, to = 2, length = 30)
)
ph1_range_effects2 <- data_frame(
  range = 1:50,
  ph1_range_effect = seq(from = -4, to = 4, length = 50)
)

ph1_phenoprep_bigspatial <- ph1_des_outprep %>%
  left_join(phen_gen) %>%
  left_join(ph1_row_effects2) %>%
  left_join(ph1_range_effects2) %>%
  mutate(residual = rnorm(1500, sd = 0.1),
         ph1_pheno = all_sum + ph1_row_effect + ph1_range_effect + residual)

saveRDS(ph1_pheno_bigspatial, file = "data/ph1_phenoprep_bigspatial.rds")

ph1_phenoprep_bigspatial2 <- ph1_des_outprep %>%
  left_join(phen_gen) %>%
  left_join(ph1_row_effects2) %>%
  left_join(ph1_range_effects2) %>%
  mutate(residual = rnorm(1500, sd = 0.5),
         ph1_pheno = all_sum + ph1_row_effect + ph1_range_effect + residual)

saveRDS(ph1_phenoprep_bigspatial2, file = "data/ph1_phenoprep_bigspatial2.rds")

# ggplot(ph1_pheno, aes(row, range)) + geom_tile(aes(fill = ph1_pheno))


#Make double rep design

treatment_fix <- data_frame(wrong_id = factor(1:750),
                            id = sample(phen_gen$id, 750, replace = FALSE))

ph1_phenodub_bigspatial <- ph1_outdub %>%
  left_join(phen_gen) %>%
  left_join(ph1_row_effects2) %>%
  left_join(ph1_range_effects2) %>%
  mutate(residual = rnorm(1500, sd = 0.1),
         ph1_pheno = all_sum + ph1_row_effect + ph1_range_effect + residual,
         id = factor(id))

saveRDS(ph1_phenodub_bigspatial, file = "data/ph1_phenodub_bigspatial.rds")

ph1_phenodub_bigspatial2 <- ph1_outdub %>%
  left_join(phen_gen) %>%
  left_join(ph1_row_effects2) %>%
  left_join(ph1_range_effects2) %>%
  mutate(residual = rnorm(1500, sd = 0.5),
         ph1_pheno = all_sum + ph1_row_effect + ph1_range_effect + residual,
         id = factor(id))

saveRDS(ph1_phenodub_bigspatial2, file = "data/ph1_phenodub_bigspatial2.rds")

#Make triple rep design
ph1_phenotrip_bigspatial <- ph1_outtrip %>%
  left_join(phen_gen) %>%
  left_join(ph1_row_effects2) %>%
  left_join(ph1_range_effects2) %>%
  mutate(residual = rnorm(1500, sd = 0.1),
         ph1_pheno = all_sum + ph1_row_effect + ph1_range_effect + residual,
         id = factor(id))

saveRDS(ph1_phenotrip_bigspatial, file = "data/ph1_phenotrip_bigspatial.rds")

ph1_phenotrip_bigspatial2 <- ph1_outtrip %>%
  left_join(phen_gen) %>%
  left_join(ph1_row_effects2) %>%
  left_join(ph1_range_effects2) %>%
  mutate(residual = rnorm(1500, sd = 0.5),
         ph1_pheno = all_sum + ph1_row_effect + ph1_range_effect + residual,
         id = factor(id))

saveRDS(ph1_phenotrip_bigspatial2, file = "data/ph1_phenotrip_bigspatial2.rds")





###########
# Phase 2 #
###########

#prep
phen_gen <- readRDS("data/phen_gen.rds")
ph1.des <- des.prep00(960, 30, 50, rep(c(2, 1), c(540, 420)), phen_gen$id,
                      ribs = 30, cibs = 25, tgrp = rep(c(2, 1), c(540, 420)))

ph2.desprep <- des.prep00(1499, 17, 100, rep(c(2, 1), c(201, 1298)),
                          tnam = as.character(1:1499),# phen_gen$id,
                          ribs = c(17), cibs = c(50), 
                          tgrp = rep(c(2, 1), c(201, 1298)))

ph2_outprep <- run(ph2.desprep)


