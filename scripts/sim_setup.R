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





############
# 2 Phases #
############
#prep
phen_gen <- readRDS("data/phen_gen.rds")

#ph1 750 trts
# 20 x 50 des

p2p1_row_effects <- data_frame(
  row = 1:20,
  p1_row_effect = seq(from = -0.5, to = 0.5, length = 20)
)
p2p1_range_effects <- data_frame(
  range = 1:50,
  p1_range_effects = seq(from = -0.75, to = 0.75, length = 50)
)

p2ph1.des <- des.prep00(750, 20, 50, rep(c(2, 1), c(250, 500)), sample(phen_gen$id, 750),
                      ribs = c(20), cibs = c(25), tgrp = rep(c(2, 1), c(250, 500)))
saveRDS(p2ph1.des, "data/p2ph1.des")
p2ph1_out <- run(p2ph1.des)
saveRDS(p2ph1_out, "data/p2ph1_out")

phph1_df <- p2ph1_out$dlist %>% 
  tbl_df() %>% 
  rename(row = ROW,
         range = RANGE) %>% 
  mutate(id = as.numeric(ID)) %>% 
  mutate(block = ceiling(range / 25)) %>% 
  left_join(p2p1_row_effects) %>% 
  left_join(p2p1_range_effects) %>% 
  group_by(id) %>% 
  mutate(p1rep = n()) %>% 
  ungroup() %>% 
  mutate(p1_plotid = 1:nrow(.)) %>% 
  mutate()

#Phase2

#ph2 10 x 120 des

p2p2_row_effects <- data_frame(
  p2_row = 1:10,
  p2_row_effect = rnorm(10, sd = 0.15) + rnorm(10, sd = 0.05)
)
p2p2_range_effects <- data_frame(
  p2_range = 1:120,
  p2_range_effects = seq(from = -0.4, to = 0.4, length = 120) + rnorm(120, sd = 0.075)
)
p1_doubles <- phph1_df$id[phph1_df$p1rep == 2] %>% unique
p1_singles <- phph1_df$id[phph1_df$p1rep == 1]
p2_doubles <- sample(p1_singles, 200)
p2_singles <- p1_singles[!p1_singles %in% p2_doubles]

p2ph2.des <- des.prep00(750, 10, 120, rep(c(2, 1), c(450, 300)),
                        c(p1_doubles, p2_doubles, p2_singles), 
                        ribs = c(10), cibs = 60, tgrp = rep(c(2, 1), c(450, 300)))
saveRDS(p2ph2.des, "data/p2ph2.des")
p2ph2_out <- run(p2ph2.des)
saveRDS(p2ph2_out, "data/p2ph2_out")

p2ph2_df <- p2ph2_out$dlist %>% 
  tbl_df() %>% 
  rename(p2_row = ROW,
         p2_range = RANGE) %>% 
  mutate(id = as.numeric(as.character(ID)),
         p2_block = ceiling(p2_range / 60)) %>% 
  left_join(p2p2_row_effects) %>% 
  left_join(p2p2_range_effects) %>% 
  mutate(p2_plotid = 1:nrow(.)) %>% 
  group_by(ID) %>% 
  mutate(p2rep = n())

phph1_reps <- phph1_df %>%
  rename(p1_block = block) %>% 
  select(p1_row = row, p1_range = range, id, contains("p1"), ID) %>% 
  filter(p1rep == 2)

p2ph2_reps <- p2ph2_df %>% 
  ungroup() %>% 
  select(contains("p2"), id) %>% 
  filter(id %in% phph1_reps$id)

p1p2_joinreps <- phph1_reps %>% 
  arrange(id) %>% 
  bind_cols(p2ph2_reps %>% 
              arrange(id) %>% 
              select(-id))

phph1_sing <-  phph1_df %>%
  rename(p1_block = block) %>% 
  select(p1_row = row, p1_range = range, id, contains("p1"), ID) %>% 
  filter(p1rep == 1) %>% 
  full_join(p2ph2_df %>% 
              ungroup() %>% 
              select(contains("p2"), id) %>% 
              filter(!id %in% phph1_reps$id))

p1p2_df <- p1p2_joinreps %>% 
  arrange(id) %>% 
  bind_rows(phph1_sing %>% 
              arrange(id)) 

#Phase 3
#ph3 15 x 100 des

p2p3_row_effects <- data_frame(
  p3_row = 1:15,
  p3_row_effect = rnorm(15, sd = 0.1)
)
p2p3_range_effects <- data_frame(
  p3_range = 1:100,
  p3_range_effects = seq(from = -0.3, to = 0.3, length = 100)
)

p2ph3.des <- DiGGer(750, 15, 100, 15, 50, TreatmentName = unique(p2ph2_df$id))

p2ph3_df <- p2ph3.des$idsgn %>%
  as.data.frame %>% 
  tbl_df() %>% 
  mutate(p3_row = 1:15) %>% 
  gather(p3_range, id, -p3_row) %>% 
  mutate(p3_range = as.numeric(sub("V", "", p3_range)),
         id = factor(id),
         p3_block = ceiling(p3_range / 50)) %>% 
  left_join(p2p3_row_effects) %>% 
  left_join(p2p3_range_effects) %>% 
  rename(id_merge = id)

p1p2_df <- p1p2_df %>% 
  mutate(id_merge = as.factor(as.numeric(as.factor(id)))) 

p1p2_dubs <- p1p2_df %>% 
  filter(p2rep == 2) %>% 
  arrange(id_merge) 

p1p3_join1 <- p2ph3_df %>% 
  filter(id_merge %in% p1p2_dubs$id_merge) %>% 
  arrange(id_merge) %>% 
  bind_cols(p1p2_dubs)

p1p3_join2 <- p2ph3_df %>% 
  filter(!id_merge %in% p1p2_dubs$id_merge) %>% 
  left_join(p1p2_df)

p1p3_df <- p1p3_join1 %>% 
  bind_rows(p1p3_join2) %>% 
  mutate(p3_spatial = p3_row_effect + p3_range_effects +
           p2_row_effect + p2_range_effects +
           p1_row_effect + p1_range_effects) %>% 
  left_join(phen_gen %>% rename(ID = id)) %>% 
  mutate(p3_pheno = p3_spatial + all_sum +
           rnorm(1500, sd = 0.2))

saveRDS(p1p3_df, "data/p1p3_df.rds")
