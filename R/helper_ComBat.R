####### ComBat and ComBat-Seq Helper Functions ######

# Check for confounding
confounding_check <- function(design, n.batch, test_name){
    if(qr(design)$rank < ncol(design)) {
        ## if(ncol(design)<=(n.batch)){
        #stop("Batch variables are redundant! Remove one or more of the
        #batch variables so they are no longer confounded")}
        if(ncol(design)==(n.batch+1)) {
            stop("The covariate is confounded with batch! Remove ",
                "the covariate and rerun ", test_name)
        }
        if(ncol(design)>(n.batch+1)) {
            if((qr(design[,-c(seq_len(n.batch))])$rank<ncol(
                design[,-c(seq_len(n.batch))]))){
                stop("The covariates are confounded! Please remove one or ",
                    "more of the covariates so the design is not confounded")
            } else {
                stop("At least one covariate is confounded with batch! Please ",
                    "remove confounded covariates and rerun ", test_name)
            }
        }
    }
}

# Check for NAs
NA_check <- function(dat){
    NAs <- any(is.na(dat))
    if(NAs){
        message(c('Found', sum(is.na(dat)), 'Missing Data Values'), sep=' ')}
    return(NAs)
}

####### ComBat Helper functions ######

# Check for only one batch
one_batch_check <- function(batch){
    if(length(dim(batch))>1){
        stop("This version of ComBat only allows one batch variable")
    } 
}
## find genes with zero variance in any of the batches
zero_variance <- function(dat, batch){
    batch <- as.factor(batch)
    zero.rows.lst <- lapply(levels(batch), function(batch_level){
        if(sum(batch==batch_level)>1){
            return(which(apply(dat[, batch==batch_level], 1, 
                function(x){var(x)==0})))
        }else{
            return(which(rep(1,3)==2))
        }
    })
    zero.rows <- Reduce(union, zero.rows.lst)
    keep.rows <- setdiff(seq_len(nrow(dat)), zero.rows)
    
    dat.orig <- NULL
    
    if (length(zero.rows) > 0) {
        message(sprintf("Found %d genes with uniform expression within a ", 
            "single batch (all zeros); these will not be adjusted for batch.",
            length(zero.rows)))
        # keep a copy of the original data matrix and remove zero var rows
        dat.orig <- dat
        dat <- dat[keep.rows, ]
    }
    
    return(list(batch = batch, dat = dat, zero.rows = zero.rows,
        keep.rows = keep.rows, dat.orig = dat.orig))
}

## Make a set of batch indicators/characteristics and create design matrix
batch_design <- function(batch, mean.only, ref.batch, mod){
    if(any(table(batch)==1)){mean.only <- TRUE}
    if(mean.only==TRUE){message("Using the 'mean only' version of ComBat")}
    
    batchmod <- model.matrix(~-1+batch)  
    if (!is.null(ref.batch)){
        if (!(ref.batch%in%levels(batch))) {
            stop("reference level ref.batch is not one of the levels of the ",
                "batch variable")}
        message("Using batch =", ref.batch, 
            "as a reference batch (this batch won't change)")
        ref <- which(levels(as.factor(batch))==ref.batch) # find the reference
        batchmod[,ref] <- 1
    } else {
        ref <- NULL
    }
    message("Found", nlevels(batch), "batches")
    
    ## A few other characteristics on the batches
    n.batch <- nlevels(batch)
    batches <- list()
    for (i in seq_len(n.batch)) {
        batches[[i]] <- which(batch == levels(batch)[i])
    } # list of samples in each batch  
    n.batches <- vapply(batches, length, numeric(1))
    if(any(n.batches==1)){
        mean.only <- TRUE
        message("Note: one batch has only one sample, setting mean.only=TRUE")
    }
    n.array <- sum(n.batches)
    
    ## combine batch variable and covariates
    design <- cbind(batchmod, mod)
    
    ## check for intercept in covariates, and drop if present
    check <- apply(design, 2, function(x) all(x == 1))
    if(!is.null(ref)){
        check[ref] <- FALSE
    } ## except don't throw away the reference batch indicator
    design <- as.matrix(design[,!check])
    
    ## Number of covariates or covariate levels
    message("Adjusting for", 
        ncol(design)-ncol(batchmod), 
        'covariate(s) or covariate level(s)')
    
    return(list(ref = ref, n.batch = n.batch, batches = batches,
        n.batches = n.batches, n.array = n.array, design = design, 
        mean.only = mean.only))
}

## Standardize Data across genes
standardize_data <- function(NAs, design, dat){
    message('Standardizing Data across genes')
    if (!NAs){
        B.hat <- solve(crossprod(design), tcrossprod(t(design),
            as.matrix(dat)))
    } else {
        B.hat <- apply(dat, 1, Beta.NA, design) # FIXME
    }
    
    return(B.hat)
}

ref_batch_variables_update <- function(ref.batch, ref, B.hat, n.batches,
    n.array, n.batch, NAs, batches, dat, design) {
    if(!is.null(ref.batch)){
        grand.mean <- t(B.hat[ref, ])
    } else {
        grand.mean <- crossprod(n.batches/n.array, B.hat[seq_len(n.batch),])
    }
    
    if (!NAs){
        if(!is.null(ref.batch)) {
            ref.dat <- dat[, batches[[ref]]]
            var.pooled <- ((ref.dat - t(design[batches[[ref]], ] %*% 
                    B.hat))^2) %*% rep(1/n.batches[ref],n.batches[ref]) # FIXME
        } else {
            var.pooled <- ((dat-t(design %*% B.hat))^2) %*% 
                rep(1/n.array,n.array) # FIXME
        }
    } else {
        if(!is.null(ref.batch)) {
            ref.dat <- dat[, batches[[ref]]]
            var.pooled <- rowVars(ref.dat-t(design[batches[[ref]], ]%*%B.hat),
                na.rm=TRUE)
        } else {
            var.pooled <- rowVars(dat-t(design %*% B.hat), na.rm=TRUE)
        }
    }
    
    stand.mean <- t(grand.mean) %*% t(rep(1,n.array)) # FIXME
    if(!is.null(design)){
        tmp <- design
        tmp[,c(seq_len(n.batch))] <- 0
        stand.mean <- stand.mean+t(tmp %*% B.hat) #FIXME
    }  
    s.data <- (dat-stand.mean)/(sqrt(var.pooled) %*% t(rep(1,n.array))) # FIXME
    
    return(list(var.pooled = var.pooled, stand.mean = stand.mean,
        s.data = s.data))
}

regression_batch_parameters <- function(design, n.batch, NAs, s.data, batches,
    mean.only){
    message("Fitting L/S model and finding priors")
    batch.design <- design[, seq_len(n.batch)]
    if (!NAs){
        gamma.hat <- solve(crossprod(batch.design), tcrossprod(t(batch.design),
            as.matrix(s.data)))
    } else{
        gamma.hat <- apply(s.data, 1, Beta.NA, batch.design) # FIXME
    }
    
    delta.hat <- NULL
    for (i in batches){
        if(mean.only==TRUE) {
            delta.hat <- rbind(delta.hat,rep(1,nrow(s.data))) 
        } else {
            delta.hat <- rbind(delta.hat, rowVars(s.data[,i], na.rm=TRUE))
        }
    }
    
    return(list(batch.design = batch.design, gamma.hat = gamma.hat, 
        delta.hat = delta.hat))
}

## Find priors; Plot emperical and parametric priors
priors_calculations <- function(gamma.hat, delta.hat, prior.plots, par.prior){
    ##Find Priors
    gamma.bar <- rowMeans(gamma.hat)
    t2 <- rowVars(gamma.hat)
    a.prior <- apply(delta.hat, 1, aprior) # FIXME 
    b.prior <- apply(delta.hat, 1, bprior) # FIXME
    
    ## Plot empirical and parametric priors
    
    if (prior.plots && par.prior) {
        old_pars <- par(no.readonly = TRUE)
        on.exit(par(old_pars))
        par(mfrow=c(2,2))
        
        ## Top left
        tmp <- density(gamma.hat[1,])
        plot(tmp,  type='l',
            main=expression(paste("Density Plot of First Batch ", hat(gamma))))
        xx <- seq(min(tmp$x), max(tmp$x), length=100)
        lines(xx, dnorm(xx, gamma.bar[1], sqrt(t2[1])), col=2)
        
        ## Top Right
        qqnorm(gamma.hat[1,], 
            main = expression(paste("Normal Q-Q Plot of First Batch ",
                hat(gamma))))
        qqline(gamma.hat[1,], col=2)
        
        ## Bottom Left
        tmp <- density(delta.hat[1,])
        xx <- seq(min(tmp$x), max(tmp$x), length=100)
        tmp1 <- list(x=xx, y=dinvgamma(xx, a.prior[1], b.prior[1]))
        plot(tmp, typ="l", ylim=c(0, max(tmp$y, tmp1$y)),
            main=expression(paste("Density Plot of First Batch ", hat(delta))))
        lines(tmp1, col=2)
        
        ## Bottom Right
        invgam <- 1/qgamma(1-ppoints(ncol(delta.hat)), a.prior[1], b.prior[1])
        qqplot(invgam, delta.hat[1,],
            main=expression(paste("Inverse Gamma Q-Q Plot of First Batch ",
                hat(delta))),
            ylab="Sample Quantiles", xlab="Theoretical Quantiles")
        lines(c(0, max(invgam)), c(0, max(invgam)), col=2)
    }
    
    return(list(a.prior = a.prior, b.prior = b.prior, gamma.bar = gamma.bar,
        t2 = t2))
}

## Find EB batch adjustments
EB_batch_adjustment <- function(n.batch, s.data, par.prior, mean.only,
    gamma.hat, gamma.bar, t2, batches, delta.hat, a.prior, b.prior, ref,
    ref.batch, BPPARAM){
    
    gamma.star <- delta.star <- matrix(NA, nrow=n.batch, ncol=nrow(s.data))
    if (par.prior) {
        message("Finding parametric adjustments")
        results <- bplapply(seq_len(n.batch), function(i) {
            if (mean.only) {
                gamma.star <- postmean(gamma.hat[i,], gamma.bar[i], 1, 1, t2[i])
                delta.star <- rep(1, nrow(s.data))
            }
            else {
                temp <- it.sol(s.data[, batches[[i]]], gamma.hat[i, ],
                    delta.hat[i, ], gamma.bar[i], t2[i], a.prior[i],
                    b.prior[i])
                gamma.star <- temp[1, ]
                delta.star <- temp[2, ]
            }
            list(gamma.star=gamma.star, delta.star=delta.star)
        }, BPPARAM = BPPARAM)
        for (i in seq_len(n.batch)) {
            gamma.star[i,] <- results[[i]]$gamma.star
            delta.star[i,] <- results[[i]]$delta.star
        }
    }
    else {
        message("Finding nonparametric adjustments")
        results <- bplapply(seq_len(n.batch), function(i) {
            if (mean.only) {
                delta.hat[i, ] <- 1
            }
            temp <- int.eprior(as.matrix(s.data[, batches[[i]]]),
                gamma.hat[i, ], delta.hat[i, ])
            list(gamma.star=temp[1,], delta.star=temp[2,])
        }, BPPARAM = BPPARAM)
        for (i in seq_len(n.batch)) {
            gamma.star[i,] <- results[[i]]$gamma.star
            delta.star[i,] <- results[[i]]$delta.star
        }
    }
    
    if(!is.null(ref.batch)){
        gamma.star[ref,] <- 0  ## set reference batch mean equal to 0
        delta.star[ref,] <- 1  ## set reference batch variance equal to 1
    }
    
    return(list(gamma.star = gamma.star, delta.star = delta.star))
}

## Normalize the Data ###
normalize_data <- function(s.data, batches, batch.design, gamma.star,
    delta.star, n.batches, var.pooled, n.array, stand.mean, ref.batch, ref, dat,
    zero.rows, keep.rows, dat.orig){
    message("Adjusting the Data\n")
    
    bayesdata <- s.data
    j <- 1
    for (i in batches){
        bayesdata[,i] <- (bayesdata[,i] - t(batch.design[i,] %*% gamma.star))/
            (sqrt(delta.star[j,])%*%t(rep(1,n.batches[j]))) # FIXME
        j <- j+1
    }
    
    bayesdata <- (bayesdata*(sqrt(var.pooled) %*% t(rep(1,n.array)))) + 
        stand.mean # FIXME
    
    ## Do not change ref batch at all in reference version
    if(!is.null(ref.batch)){
        bayesdata[, batches[[ref]]] <- dat[, batches[[ref]]]
    }
    
    ## put genes with 0 variance in any batch back in data
    if (length(zero.rows) > 0) {
        dat.orig[keep.rows, ] <- bayesdata
        bayesdata <- dat.orig
    }
    
    return(bayesdata)
}