# This tests the twelve main types of error thrown in
# helperTestDataMetricsLitreApp.R when a user inputs a dataMetrics object
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
dME = data.frame(ID = ID(), FDR = Nu(), logFC = Nu(), PVal = Nu(),
Qvar = Nu())
dME = Ch(dME)
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
dMFail10 = dM
dMFail10[["A_B"]][["ID"]][[100]] = "ID101"

test_data <- function() {
    checkException(helperTestDataMetricsLitreApp(d, dMFail1))
    checkException(helperTestDataMetricsLitreApp(d, dMFail2))
    checkException(helperTestDataMetricsLitreApp(d, dMFail3))
    checkException(helperTestDataMetricsLitreApp(d, dMFail4))
    checkException(helperTestDataMetricsLitreApp(d, dMFail5))
    checkException(helperTestDataMetricsLitreApp(d, dMFail6))
    checkException(helperTestDataMetricsLitreApp(d, dMFail7))
    checkException(helperTestDataMetricsLitreApp(d, dMFail8))
    checkException(helperTestDataMetricsLitreApp(d, dMFail9))
    checkException(helperTestDataMetricsLitreApp(d, dMFail10))
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
dMFail8 = dM
names(dMFail8) = c("A_B", "A_E", "B_E")
dMFail9 = dM
dMFail9[["A_B"]][["ID"]][[100]] = "ID101"
dMFail10 = dM
names(dMFail10) = c("A_A", "B_B", "C_C")

test_data <- function() {
    checkException(helperTestDataMetricsLitreApp(d, dMFail2))
    checkException(helperTestDataMetricsLitreApp(d, dMFail3))
    checkException(helperTestDataMetricsLitreApp(d, dMFail4))
    checkException(helperTestDataMetricsLitreApp(d, dMFail5))
    checkException(helperTestDataMetricsLitreApp(d, dMFail6))
    checkException(helperTestDataMetricsLitreApp(d, dMFail7))
    checkException(helperTestDataMetricsLitreApp(d, dMFail8))
    checkException(helperTestDataMetricsLitreApp(d, dMFail9))
    checkException(helperTestDataMetricsLitreApp(d, dMFail10))
}
