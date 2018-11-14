set.seed(1)
data(soybean_ir)
keep <- sample(1:nrow(soybean_ir), nrow(soybean_ir)/10, replace=FALSE)
soybean_ir_sub <- soybean_ir[keep,]
save(soybean_ir_sub, file = "soybean_ir_sub.rda")
