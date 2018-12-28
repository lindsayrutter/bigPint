helperSMPrint <- function(ret, outDir, fName){

seqVec = seq_along(ret)

invisible(lapply(seqVec, function(x){
    fileName = paste0(outDir, "/", names(ret[x]), fName)
    jpeg(filename=fileName, height=900, width=900)
    print(ret[[x]])
    dev.off()
}))

}
