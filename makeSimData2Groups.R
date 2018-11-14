library(PROPER)
set.seed(1)

fun.lfc=function(x) rnorm(x, mean=0, sd=3)
simOptions = RNAseq.SimOptions.2grp(ngenes=10000, p.DE = 0.5, lfc=fun.lfc, lBaselineExpr="bottomly")
data = simRNAseq(simOptions, n1=3, n2=3)
simData = as.data.frame(data[["counts"]])
simData$ID = paste0("ID", 1:nrow(simData))
simData = simData[,c(7,1:6)]
colnames(simData) = c("ID", paste0(rep(c("A.", "B."), each = 3), 1:3))

simData[,-1] <- simData[,-1] + rep(runif(6,1,3), nrow(simData))
sdat2G <- simData
save(sdat2G, file = "sdat2G.rda")
