######################################################################
##
## a list of functions for testing DML from Bisulfite seq data in methrix format
##
######################################################################

############################################
## wrapper function for DML test for methrix
############################################

DMLtest.methrix <- function(meth, group1, group2, equal.disp=FALSE, smoothing=FALSE,
                    smoothing.span=500, ncores) {

    ## determine number of cores used in parallel computing setting
    if(missing(ncores)) {
        if(.Platform$OS.type == "windows" | Sys.info()['sysname'] == "Windows")
            ## Windows, use single core
            ncores = 1
        else # for Mac and linux
            ncores = max(detectCores() - 3, 1)
    }

    if(ncores > detectCores())
        stop("ncores exceeds the number of CPU cores on the system.")

    ## grab two group data
    tmp <- getBSseqIndex(colnames(meth), group1, group2)
    ## sampleNames don't work on methrix. It should be either colnames, or we should implement this function.
    meth1 <- meth[,tmp$group1]
    meth2 <- meth[,tmp$group2]

    ## remove loci with all 0 coverages in a condition
    ## It's not required that all replicates have coverage.
    ## But there must be some coverage from at least one replicate.
    n1 <- as.array(get_matrix(meth1, "C"))
    n1[is.na(n1)] <- 0
    n2 <- as.array(get_matrix(meth2, "C"))
    n2[is.na(n2)] <- 0
    ## remove if a rather long strech (like 200 bps) of regions have no coverage
    allpos <- meth@elementMetadata$start
    ix1 <- hasCoverage(n1, allpos)
    ix2 <- hasCoverage(n2, allpos)
    ix <- ix1 & ix2
    meth1 <- meth1[ix]
    meth2 <- meth2[ix]

    ## Check the consistence of inputs
    nreps1 <- dim(meth1)[2]
    nreps2<- dim(meth2)[2]
    if( (nreps1==1 | nreps2==1) & !equal.disp ) { ## singel replicate case, and unequal dispersion in two groups
        if( !smoothing )
            stop("There is no biological replicates in at least one condition. Please set smoothing=TRUE or equal.disp=TRUE and retry.")
    }

    if(!smoothing) { ## no smoothing.
        dmls <- DMLtest.noSmoothMethrix(meth1, meth2, equal.disp, ncores)
    } else { ## smoothing version
        dmls <- DMLtest.SmoothMethrix(meth1, meth2, equal.disp, smoothing.span, ncores)
    }

    class(dmls) = c("DMLtest", class(dmls))
    invisible(dmls)
}

######################################
## test DML without smoothing
######################################
DMLtest.noSmoothMethrix <- function(meth1, meth2, equal.disp, ncores) {
    ## grab counts
#coverage
    n1 <- as.array(get_matrix(meth1, type="C"))
# number of methylated reads
    x1 <- as.array(get_matrix(meth1, type="C")*get_matrix(meth1, type="M"))
    n1[is.na(n1)] <- 0
    x1[is.na(x1)] <- 0

    #coverage
    n2 <- as.array(get_matrix(meth2, type="C"))
    # number of methylated reads
    x2 <- as.array(get_matrix(meth2, type="C")*get_matrix(meth2, type="M"))
    n2[is.na(n2)] <- 0
    x2[is.na(x2)] <- 0

    nreps1 <- ncol(x1)
    nreps2 <- ncol(x2)

    ## estimate means
    estprob1 <- compute.mean.noSmooth(x1, n1)
    estprob2 <- compute.mean.noSmooth(x2, n2)

    ## estimate dispersion
    ## - this part is slow. Could be computed parallely. Will implement later.
    cat("Estimating dispersion for each CpG site, this will take a while ...\n")
    if(equal.disp | nreps1==1 | nreps2==1) {
        phi1 <- phi2 <- est.dispersion.BSseq(cbind(x1,x2), cbind(n1,n2),
                                             cbind(estprob1, estprob2), ncores)
    } else {
        phi1 <- est.dispersion.BSseq(x1, n1, estprob1, ncores)
        phi2 <- est.dispersion.BSseq(x2, n2, estprob2, ncores)
    }


    ## weight the counts
    wt1 <- 1 / (1+(n1-1)*phi1);    wt1 <- wt1 / mean(wt1)
    wt2 <- 1 / (1+(n2-1)*phi2);    wt2 <- wt2 / mean(wt2)
    x1.wt <- x1*wt1
    n1.wt <- n1*wt1
    x2.wt <- x2*wt2
    n2.wt <- n2*wt2

    ## re-estimate means
    estprob1 <- compute.mean.noSmooth(x1.wt, n1.wt)
    estprob2 <- compute.mean.noSmooth(x2.wt, n2.wt)

    ## perform Wald test
    allchr <- as.character(meth1@elementMetadata$chr)
    allpos <- meth1@elementMetadata$start
    wald <- waldTest.DML(x1.wt, n1.wt, estprob1, phi1, x2.wt, n2.wt, estprob2, phi2,
                         smoothing=FALSE, allchr=allchr, allpos=allpos)

    return(wald)

}

######################################
## test DML with smoothing
######################################
DMLtest.SmoothMethrix <- function(meth1, meth2, equal.disp, smoothing.span, ncores) {
    ## grab counts
  ## grab counts
  #coverage
  n1 <- as.array(get_matrix(meth1, type="C"))
  # number of methylated reads
  x1 <- as.array(get_matrix(meth1, type="C")*get_matrix(meth1, type="M"))
  n1[is.na(n1)] <- 0
  x1[is.na(x1)] <- 0

  #coverage
  n2 <- as.array(get_matrix(meth2, type="C"))
  # number of methylated reads
  x2 <- as.array(get_matrix(meth2, type="C")*get_matrix(meth2, type="M"))
  n2[is.na(n2)] <- 0
  x2[is.na(x2)] <- 0

    nreps1 <- ncol(x1)
    nreps2 <- ncol(x2)
    allchr <- as.character(meth1@elementMetadata$chr)
    allpos <- meth1@elementMetadata$start

    ## Smoothing
    cat("Smoothing ...\n")
    estprob1 <- compute.mean.Smooth(x1, n1, allchr, allpos, smoothing.span)
    estprob2 <- compute.mean.Smooth(x2, n2, allchr, allpos, smoothing.span)

    ## estimate priors from counts
    cat("Estimating dispersion for each CpG site, this will take a while ...\n")
    if(equal.disp) {
        phi1 <- phi2 <- est.dispersion.BSseq(cbind(x1,x2), cbind(n1,n2), cbind(estprob1, estprob2), ncores)
    } else {
        phi1 <- est.dispersion.BSseq(x1, n1, estprob1, ncores)
        phi2 <- est.dispersion.BSseq(x2, n2, estprob2, ncores)
    }

    ## update counts - weight by dispersion
    wt1 <- 1 / (1+(n1-1)*phi1); wt1 <- wt1 / mean(wt1, na.rm=TRUE)
    wt2 <- 1 / (1+(n2-1)*phi2); wt2 <- wt2 / mean(wt2, na.rm=TRUE)
    x1.wt <- x1*wt1
    n1.wt <- n1*wt1
    x2.wt <- x2*wt2
    n2.wt <- n2*wt2

    ## re-estimate means
    estprob1 <- compute.mean.Smooth(x1.wt, n1.wt, allchr, allpos, smoothing.span)
    estprob2 <- compute.mean.Smooth(x2.wt, n2.wt, allchr, allpos, smoothing.span)

    cat("Computing test statistics ...\n")
    wald <- waldTest.DML(x1.wt, n1.wt, estprob1, phi1, x2.wt, n2.wt, estprob2, phi2,
                         smoothing=TRUE, smoothing.span, allchr=allchr, allpos=allpos)
##     wald <- waldTest.DML(x1, n1, estprob1, phi1, x2, n2, estprob2, phi2,
##                          smoothing=TRUE, smoothing.span, allchr=allchr, allpos=allpos)

    return(wald)
}

