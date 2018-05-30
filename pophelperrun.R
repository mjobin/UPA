library(pophelper)



args <- commandArgs(TRUE)

if (length(args) < 1) {
    stop('Missing argument(s).')
}

file <- as.character(args[1])
basefile <- tools::file_path_sans_ext(file)


famfile <- paste(basefile, ".fam", sep = "")
famtbl=read.table(famfile)
famlist <- famtbl$V1
famdf <- data.frame(sapply(famlist,as.character), stringsAsFactors=FALSE)
names(famdf)[1] <- "grp"

clist <- list(
"standard_12"=c("#2121D9","#9999FF","#DF0101","#04B404","#FFFB23","#FF9326","#A945FF","#0089B2","#B26314","#610B5E","#FE2E9A","#BFF217"),
"rich.colors"=pophelper:::getColours(13))


setwd(paste(basefile, "/BEST", sep = ""))
sfiles <- list.files(pattern="\\.Q$",full.names=TRUE)
slist <- readQ(sfiles)

# plotQ(slist,imgoutput="join", sppos="left", showindlab=F, grplab=famdf,  outputfilename="plotq", imgtype="png")

plotQ(slist,imgoutput="join", outputfilename="plotq", imgtype="png")