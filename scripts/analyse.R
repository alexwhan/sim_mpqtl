library(mpwgaim)
library(asreml)

ph1_pheno <- readRDS("data/ph1_pheno.rds")
mp4_int <- readRDS("data/mp4_int.rds")
ph1.asr <- asreml(ph1_pheno ~ 1, ~ id + row + range, data = ph1_pheno)
summary(ph1.asr)$varcomp

gen_mpw <- mpwgaim(ph1.asr, ph1_pheno, mp4_int, "id")
ph1_pheno <- as.data.frame(ph1_pheno)
gen_pred <- predict(gen_mod, classify = "id")
gen_pred$predictions$pvals
