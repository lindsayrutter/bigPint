helperSMPrint <- function(ret, outDir, fName){

seqVec = seq(1, length(ret))
for (i in seq_along(seqVec)){
    fileName = paste0(outDir, "/", names(ret[i]), fName)
    jpeg(filename=fileName, height=900, width=900)
    print(ret[[i]])
    dev.off()
}
}
