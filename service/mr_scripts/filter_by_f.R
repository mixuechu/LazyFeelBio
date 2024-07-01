
    dat <- read.csv("exposure.LD.csv")
    Ffilter <- 10
    N <- dat[1, "samplesize.exposure"]
    dat <- transform(dat, R2 = 2*((beta.exposure)^2)*eaf.exposure*(1-eaf.exposure))
    dat <- transform(dat, F = (N-2)*R2/(1-R2))
    outTab <- dat[dat$F > Ffilter,]
    write.csv(dat, "exposure.F.csv", row.names=FALSE)
    