
library(ggplot2)
library(scatterplot3d)
library(SNPRelate)


args <- commandArgs(TRUE)

if (length(args) < 1) {
    stop('Missing argument(s).') 
}
print(args)
filebase <- args[1]
threads <- as.integer(args[2])

fileped <- paste(filebase, ".ped", sep = "")

print(fileped)

snpgdsPED2GDS(paste(filebase, ".ped", sep = ""), paste(filebase, ".map", sep = ""), paste(filebase, ".gds", sep = ""))

var <- snpgdsOpen(paste(filebase, ".gds", sep = ""))


varpca <- snpgdsPCA(var, autosome.only=FALSE, num.thread=threads)
sample.id <- read.gdsn(index.gdsn(var, "sample.id"))
pop_code <- read.gdsn(index.gdsn(var, "sample.annot/family"))
head(cbind(sample.id, pop_code))

pc.p <- varpca$varprop*100
pc.percent <- head(round(pc.p, 2))



vartab <- data.frame(sample.id = varpca$sample.id,
    pop = factor(pop_code)[match(varpca$sample.id, sample.id)],
    EV1 = varpca$eigenvect[,1],    
    EV2 = varpca$eigenvect[,2],

	EV3 = varpca$eigenvect[,3],
    # EV4 = varpca$eigenvect[,4],
    stringsAsFactors = FALSE)
    
    
pdf(paste(filebase, "-pca-1-2.pdf", sep = ""))
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(vartab$EV1, vartab$EV2, pch=as.integer(vartab$pop), col=as.integer(vartab$pop), xlab=paste("PC1: ", pc.percent[1], "%", sep = ""), ylab=paste("PC2: ", pc.percent[2], "%", sep = ""), main=filebase)
legend("topright", inset=c(-0.2,0), legend=levels(vartab$pop), pch=1:nlevels(vartab$pop), col=1:nlevels(vartab$pop))
dev.off()


# pdf(paste(filebase, "-pca-3d.pdf", sep = ""))
# scatterplot3d(vartab$EV1, vartab$EV2, vartab$EV3, pch=as.integer(vartab$pop), color=as.integer(vartab$pop), xlab=paste("PC1: ", pc.percent[1], "%", sep = ""), ylab=paste("PC2: ", pc.percent[2], "%", sep = ""), zlab=paste("PC3: ", pc.percent[3], "%", sep = ""))
# legend("topright", inset=c(-0.2,0), legend=levels(vartab$pop), pch=1:nlevels(vartab$pop), col=1:nlevels(vartab$pop))
# dev.off()

# pdf(paste(filebase, "-pca-3-4.pdf", sep = ""))
# plot(vartab$EV3, vartab$EV4, pch=as.integer(vartab$pop), col=as.integer(vartab$pop), xlab=paste("PC3: ", pc.percent[3], "%", sep = ""), ylab=paste("PC4: ", pc.percent[4], "%", sep = ""), main=filebase)
# legend("topright", inset=c(-0.2,0), legend=levels(vartab$pop), pch=1:nlevels(vartab$pop), col=1:nlevels(vartab$pop))
# dev.off()










