args <- commandArgs(TRUE)

if (length(args) < 1) {
    stop('Missing argument(s).') 
}
infile <- args[1]
basefile <- tools::file_path_sans_ext(infile)
pdf(paste(basefile, ".pdf", sep = ""))
cvs <- read.csv(infile, header=TRUE)
cvsmeans <- colMeans(cvs)
plot(cvsmeans,xlab="K", ylab="Mean CV")
lines(cvsmeans)
dev.off()
