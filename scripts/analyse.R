library(mpwgaim)
library(asreml)
mp4_int <- readRDS("data/mp4_int.rds")

ph1_phenoprep <- readRDS("data/ph1_phenoprep.rds")
ph1.asr <- asreml(ph1_pheno ~ 1, ~ id + row + range + block, data = ph1_phenoprep)
ph1.asrprep <- mpwgaim(ph1.asr, ph1_phenoprep, mp4_int, "id")

ph1_phenodub <- readRDS("data/ph1_phenodub.rds")
ph1.asrdub <- asreml(ph1_pheno ~ 1, ~ id + row + range + rep + block, data = ph1_phenodub)
ph1.wgdub <- mpwgaim(ph1.asrdub, ph1_phenodub, mp4_int, "id")

ph1_phenotrip <- readRDS("data/ph1_phenotrip.rds")
ph1.asrtrip <- asreml(ph1_pheno ~ 1, ~ id + row + range + rep + block, data = ph1_phenotrip)
ph1.wgtrip <- mpwgaim(ph1.asrtrip, ph1_phenotrip, mp4_int, "id")
