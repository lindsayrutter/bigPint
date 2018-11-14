set.seed(1)
data(soybean_cn)
keep <- sample(1:nrow(soybean_cn), nrow(soybean_cn)/10, replace=FALSE)
soybean_cn_sub <- soybean_cn[keep,]
save(soybean_cn_sub, file = "soybean_cn_sub.rda")
