ped8 <- ped

qtl_chroms <- sample(length(map), 30, TRUE)
qtl_locs <- pick_loc(map, qtl_chroms)
set.seed(123)
qtl_sizes <- rexp(30, rate = 0.1) + 3
qtl_sizes <- qtl_sizes / (sum(qtl_sizes))

qtl_model8 <- qtl_model
cr8 <- cr
cr.dat8 <- cr.dat