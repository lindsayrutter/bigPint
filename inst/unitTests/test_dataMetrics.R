# This tests the twelve main types of error thrown in helperTestDataMetrics.R
# when a user inputs a dataMetrics object that does not fit the expected form

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
    sample(letters, 100, replace = T)
}
Ch <- function(d){
    d$ID = as.character(d$ID)
    return(as.data.frame(d))
}
dME <- function(){
    dME <- data.frame(ID = ID(), FDR = Nu(), logFC = Nu(), PVal = Nu(),
    Qvar = Nu())
    dME$ID = as.character(dME$ID)
    return(dME)
}

# Test for case of dataMetrics with one list element
d = data.frame(ID = ID(), A.1 = Nu(), A.2 = Nu(), B.1 = Nu(), B.2 = Nu())
d = Ch(d)
dM = list()
dM[["A_B"]] = dME

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
dMFail11 = dM
dMFail11[["A_B"]]["FDR"] = Le()
dMFail12 = dM
dMFail12[["A_B"]][["ID"]][[100]] = "ID101"

test_data <- function() {
    checkException(helperTestDataMetrics(d, dMFail1, "FDR"))
    checkException(helperTestDataMetrics(d, dMFail2, "FDR"))
    checkException(helperTestDataMetrics(d, dMFail3, "FDR"))
    checkException(helperTestDataMetrics(d, dMFail4, "FDR"))
    checkException(helperTestDataMetrics(d, dMFail5, "FDR"))
    checkException(helperTestDataMetrics(d, dMFail6, "FDR"))
    checkException(helperTestDataMetrics(d, dMFail7, "FDR"))
    checkException(helperTestDataMetrics(d, dMFail8, "FDR"))
    checkException(helperTestDataMetrics(d, dMFail9, "FDR"))
    checkException(helperTestDataMetrics(d, dM, "noVar"))
    checkException(helperTestDataMetrics(d, dMFail11, "FDR"))
    checkException(helperTestDataMetrics(d, dMFail12, "FDR"))
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
dMFail11 = dM
dMFail11[["A_B"]]["FDR"] = Le()
dMFail12 = dM
dMFail12[["A_B"]][["ID"]][[100]] = "ID101"
dMFail13 = dM
names(dMFail13) = c("A_A", "B_B", "C_C")

test_data <- function() {
    checkException(helperTestDataMetrics(d, dMFail2, "FDR"))
    checkException(helperTestDataMetrics(d, dMFail3, "FDR"))
    checkException(helperTestDataMetrics(d, dMFail4, "FDR"))
    checkException(helperTestDataMetrics(d, dMFail5, "FDR"))
    checkException(helperTestDataMetrics(d, dMFail6, "FDR"))
    checkException(helperTestDataMetrics(d, dMFail7, "FDR"))
    checkException(helperTestDataMetrics(d, dMFail9, "FDR"))
    checkException(helperTestDataMetrics(d, dM, "noVar"))
    checkException(helperTestDataMetrics(d, dMFail11, "FDR"))
    checkException(helperTestDataMetrics(d, dMFail12, "FDR"))
    checkException(helperTestDataMetrics(d, dMFail13, "FDR"))
}