
    data <- read.csv("exposure.pvalue.csv")
    data <- data[,c("SNP", "chr.exposure", "pos.exposure", "pval.exposure")]
    colnames(data) <- c("SNP", "CHR", "BP", "pvalue")

    # Manhattan plot
    library(CMplot)
    CMplot(data, plot.type="m",
           LOG10=TRUE, threshold=5e-08, threshold.lwd=3, threshold.lty=1, signal.cex=0.2,
           chr.den.col=NULL, cex=0.2, bin.size=1e5, ylim=c(0,50),
           file="pdf", file.output=TRUE, width=15, height=9, verbose=TRUE)
    CMplot(data, plot.type="c",
           LOG10=TRUE, threshold=5e-08, threshold.lwd=3, threshold.lty=1, signal.cex=0.2,
           chr.den.col=NULL, cex=0.2, bin.size=1e5, ylim=c(0,100),
           file="pdf", file.output=TRUE, width=7, height=7, verbose=TRUE)
    