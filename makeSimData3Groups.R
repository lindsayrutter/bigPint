library(PROPER)
set.seed(1)

fun.lfc=function(x) rnorm(x, mean=0, sd=3)
simOptions = RNAseq.SimOptions.2grp(ngenes=8000, p.DE = 0.3, lfc=fun.lfc, lBaselineExpr="bottomly")
data = simRNAseq(simOptions, n1=3, n2=3)
fun.lfc=function(x) rnorm(x, mean=0, sd=2)
simOptions = RNAseq.SimOptions.2grp(ngenes=8000, p.DE = 0.3, lfc=fun.lfc, lBaselineExpr="bottomly")
data2 = simRNAseq(simOptions, n1=3, n2=3)
simData = as.data.frame(data[["counts"]])
simData2 = as.data.frame(data2[["counts"]])
simData = cbind(simData, simData2[,1:3])
simData$ID = paste0("ID", 1:nrow(simData))
simData = simData[ , c(ncol(simData),1:(ncol(simData)-1))]
colnames(simData) = c("ID", paste0(rep(c("C.", "D.", "E."), each = 3), 1:3))

#simData[,-1] <- log(simData[,-1]+1)
sdat3G <- simData
save(sdat3G, file = "sdat3G.rda")
