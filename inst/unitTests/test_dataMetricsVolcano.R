# This tests sixteen main types of error thrown in
# helperTestDataMetricsVolcano.R when a user inputs a dataMetrics object
# that does not fit the expected form

ID <- function(){
    paste0("ID", 1:100)
}
IDRep <- function(){
    paste0("ID", c(1:99, 1))
}
Nu <- function(){
    rnorm(1:100)
}
Le <- function(){
    sample(letters, 100, replace = TRUE)
}
Ch <- function(d){
    d$ID = as.character(d$ID)
    return(as.data.frame(d))
}
dME <- function(){
    dME <- data.frame(ID = ID(), FDR = Nu(), logFC = Nu(), PValue = Nu(),
    Qvar = Nu())
    dME$ID = as.character(dME$ID)
    return(dME)
}

# Test for case of dataMetrics with one list element
d = data.frame(ID = ID(), A.1 = Nu(), A.2 = Nu(), B.1 = Nu(), B.2 = Nu())
d = Ch(d)
dM = list()
dM[["A_B"]] = dME()

dMFail1 = d
dMFail2 = dM
dMFail2[["A_B"]] = ID()
dMFail3 = dM
dMFail3[["C_D"]] = d
dMFail4 = dM
colnames(dMFail4[["A_B"]])[1] = "notID"
dMFail5 = dM
dMFail5[["A_B"]][1] = as.factor(dMFail5[["A_B"]][1])
dMFail6 = dM
dMFail6[["A_B"]][1] = Le()
dMFail7 = dM
names(dMFail7) = c("A.B")
dMFail8 = dM
names(dMFail8) = c("A_A")
dMFail9 = dM
names(dMFail9) = c("C_D")
dMFail10 = dM
dMFail10[["A_B"]]["FDR"] = Le()
dMFail11 = dM
dMFail11[["A_B"]]["PValue"] = Le()
dMFail12 = dM
dMFail12[["A_B"]]["logFC"] = Le()
dMFail13 = dM
dMFail13[["A_B"]][["ID"]][[100]] = "ID101"

test_data <- function() {
    checkException(helperTestDataMetricsVolcano(d, dMFail1, "FDR",
    "PValue", "logFC"))
    checkException(helperTestDataMetricsVolcano(d, dMFail2, "FDR",
    "PValue", "logFC"))
    checkException(helperTestDataMetricsVolcano(d, dMFail3, "FDR",
    "PValue", "logFC"))                                                
    checkException(helperTestDataMetricsVolcano(d, dMFail4, "FDR",
    "PValue", "logFC"))                      
    checkException(helperTestDataMetricsVolcano(d, dMFail5, "FDR",
    "PValue", "logFC"))
    checkException(helperTestDataMetricsVolcano(d, dMFail6, "FDR",
    "PValue", "logFC"))                                                
    checkException(helperTestDataMetricsVolcano(d, dMFail7, "FDR",
    "PValue", "logFC"))                                                  
    checkException(helperTestDataMetricsVolcano(d, dMFail8, "FDR",
    "PValue", "logFC"))                                                   
    checkException(helperTestDataMetricsVolcano(d, dMFail9, "FDR",
    "PValue", "logFC"))   
    checkException(helperTestDataMetricsVolcano(d, dM, "noVar", "PValue",
    "logFC"))   
    checkException(helperTestDataMetricsVolcano(d, dM, "FDR", "PValue",
    "noVar"))
    checkException(helperTestDataMetricsVolcano(d, dM, "FDR", "noVar",
    "logFC"))   
    checkException(helperTestDataMetricsVolcano(d, dMFail10, "FDR",
    "PValue", "logFC"))
    checkException(helperTestDataMetricsVolcano(d, dMFail11, "FDR",
    "PValue", "logFC"))
    checkException(helperTestDataMetricsVolcano(d, dMFail12, "FDR",
    "PValue", "logFC"))
    checkException(helperTestDataMetricsVolcano(d, dMFail13, "FDR",
    "PValue", "logFC"))
}

# Test for case of dataMetrics with three list elements
d = data.frame(ID = ID(), A.1 = Nu(), A.2 = Nu(), B.1 = Nu(), B.2 = Nu(),
C.1 = Nu(), C.2 = Nu(), C.3 = Nu())
d = Ch(d)
dM = list()
dM[["A_B"]] = dME()
dM[["A_C"]] = dME()
dM[["B_C"]] = dME()

dMFail2 = dM
dMFail2[["A_B"]] = ID()
dMFail3 = dM
dMFail3[["E_F"]] = dME()
dMFail4 = dM
colnames(dMFail4[["A_B"]])[1] = "notID"
dMFail5 = dM
dMFail5[["A_B"]][1] = as.factor(dMFail5[["A_B"]][1])
dMFail6 = dM
dMFail6[["A_B"]][1] = Le()
dMFail7 = dM
names(dMFail7) = c("A.B")
dMFail9 = dM
names(dMFail9) = c("A_B", "A_E", "B_E")
dMFail10 = dM
dMFail10[["A_B"]]["FDR"] = Le()
dMFail11 = dM
dMFail11[["A_B"]]["PValue"] = Le()
dMFail12 = dM
dMFail12[["A_B"]]["logFC"] = Le()
dMFail13 = dM
dMFail13[["A_B"]][["ID"]][[100]] = "ID101"

test_data <- function() {
    checkException(helperTestDataMetricsVolcano(d, dMFail2, "FDR",
    "PValue", "logFC"))
    checkException(helperTestDataMetricsVolcano(d, dMFail3, "FDR",
    "PValue", "logFC"))                                                
    checkException(helperTestDataMetricsVolcano(d, dMFail4, "FDR",
    "PValue", "logFC"))                                                 
    checkException(helperTestDataMetricsVolcano(d, dMFail5, "FDR",
    "PValue", "logFC"))                                                   
    checkException(helperTestDataMetricsVolcano(d, dMFail6, "FDR",
    "PValue", "logFC"))                                                 
    checkException(helperTestDataMetricsVolcano(d, dMFail7, "FDR",
    "PValue", "logFC"))                                                     
    checkException(helperTestDataMetricsVolcano(d, dMFail9, "FDR",
    "PValue", "logFC")) 
    checkException(helperTestDataMetricsVolcano(d, dM, "noVar", "PValue",
    "logFC"))   
    checkException(helperTestDataMetricsVolcano(d, dM, "FDR", "PValue",
    "noVar"))
    checkException(helperTestDataMetricsVolcano(d, dM, "FDR", "noVar",
    "logFC"))   
    checkException(helperTestDataMetricsVolcano(d, dMFail10, "FDR",
    "PValue", "logFC"))
    checkException(helperTestDataMetricsVolcano(d, dMFail11, "FDR",
    "PValue", "logFC"))
    checkException(helperTestDataMetricsVolcano(d, dMFail12, "FDR",
    "PValue", "logFC"))
    checkException(helperTestDataMetricsVolcano(d, dMFail13, "FDR",
    "PValue", "logFC"))
}
