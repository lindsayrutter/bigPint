# This tests the nine main types of error thrown in helperTestData.R
# when a user inputs a data object that does not fit the expected form

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

dFail1 = list(ID = ID())
dFail = data.frame(ID = ID(), A_1 = Nu(), A_2 = Nu(), B.1 = Nu(),
B.2 = Nu())
dFail2 = Ch(dFail)
dFail3 = data.frame(notID = ID())
dFail4 = data.frame(ID = ID())
dFail5 = data.frame(ID = IDRep())
dFail5 = Ch(dFail5)
dFail6 = data.frame(ID = ID(), A.1 = Nu(), A.2 = Nu(), B.1 = Nu())
dFail6 = Ch(dFail6)
dFail = data.frame(ID = ID(), A.1 = Le(), A.2 = Le(), B.1 = Le(),
B.2 = Le())
dFail7 = Ch(dFail)
dFail = data.frame(ID = ID(), A.1 = Nu(), A.2 = Nu(), A.3 = Nu(),
A.4 = Nu())
dFail8 = Ch(dFail)
dFail = data.frame(ID = ID(), A.1 = Nu(), A.2 = Nu(), A.3 = Nu(),
B.1 = Nu())
dFail9 = Ch(dFail)

test_data <- function() {
    checkException(helperTestData(dFail1))
    checkException(helperTestData(dFail2))
    checkException(helperTestData(dFail3))
    checkException(helperTestData(dFail4))  
    checkException(helperTestData(dFail5))
    checkException(helperTestData(dFail6))
    checkException(helperTestData(dFail7))
    checkException(helperTestData(dFail8))
    checkException(helperTestData(dFail9))
}
