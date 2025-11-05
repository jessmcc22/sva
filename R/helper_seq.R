    ########  Preparation  ######## 
preparation <- function(batch, counts){
    ## Does not support 1 sample per batch yet
    batch <- as.factor(batch)
    if(any(table(batch)<=1)){
        stop("ComBat-seq doesn't support 1 sample per batch yet")
    }
    
    ## Remove genes with only 0 counts in any batch
    keep_lst <- lapply(levels(batch), function(b){
        which(apply(counts[, batch==b], 1, function(x){!all(x==0)}))
    })
    keep <- Reduce(intersect, keep_lst)
    rm <- setdiff(seq_len(nrow(counts)), keep)
    countsOri <- counts
    counts <- counts[keep, ]
    
    # require bioconductor 3.7, edgeR 3.22.1
    dge_obj <- DGEList(counts=counts)
    
    ## Prepare characteristics on batches
    n_batch <- nlevels(batch)  # number of batches
    batches_ind <- lapply(seq_len(n_batch),
        function(i){
            which(batch == levels(batch)[i])}) # list of samples in each batch  
    n_batches <- vapply(batches_ind, length, numeric(1))
    n_sample <- sum(n_batches)
    message("Found ", n_batch, ' batches')
    
    return(list(batch = batch, keep = keep, rm = rm, countsOri = countsOri, 
        counts = counts, dge_obj = dge_obj, n_batch = n_batch,
        batches_ind = batches_ind, n_batches = n_batches, n_sample = n_sample))
}


#### Make Design Matrix
create_design_matrix <- function(batch, group, counts, covar_mod, full_mod){
    # batch
    batchmod <- model.matrix(~-1+batch)  # colnames: levels(batch)
    # covariate
    group <- as.factor(group)
    if(full_mod & nlevels(group)>1){
        message("Using full model in ComBat-seq.")
        mod <- model.matrix(~group)
    }else{
        message("Using null model in ComBat-seq.")
        mod <- model.matrix(~1, data=as.data.frame(t(counts)))
    }
    # drop intercept in covariate model
    if(!is.null(covar_mod)){
        if(is.data.frame(covar_mod)){
            covar_mod <- do.call(cbind,
                lapply(seq_len(ncol(covar_mod)), 
                    function(i){
                        model.matrix(~covar_mod[,i])
                    }
                )
            )
        }
        covar_mod <- covar_mod[, !apply(covar_mod, 2, function(x){all(x==1)})]
    }
    # bind with biological condition of interest
    mod <- cbind(mod, covar_mod)
    # combine
    design <- cbind(batchmod, mod)
    
    ## Check for intercept in covariates, and drop if present
    check <- apply(design, 2, function(x) all(x == 1))
    
    design <- as.matrix(design[,!check])
    message("Adjusting for ", 
        ncol(design)-ncol(batchmod),
        ' covariate(s) or covariate level(s)')
    return(list(design = design, batchmod = batchmod, mod = mod))
}

####  Expand a vector into matrix (columns as the original vector)
vec2mat <- function(vec, n_times){
    return(matrix(rep(vec, n_times), ncol=n_times, byrow=FALSE))
}

########  Estimate gene-wise dispersions within each batch  ########
estimate_gene_wise_dispersions <- function(n_batch, n_batches, design, batchmod,
    mod, batches_ind, counts, batch){
    message("Estimating dispersions")
    ## Estimate common dispersion within each batch as an initial value
    disp_common <- vapply(seq_len(n_batch), 
        function(i){
            if((n_batches[i] <= ncol(design)-ncol(batchmod)+1) | 
                    qr(mod[batches_ind[[i]], ])$rank < ncol(mod)){ 
                # not enough residual degree of freedom
                return(estimateGLMCommonDisp(counts[, batches_ind[[i]]], 
                    design=NULL, 
                    subset=nrow(counts)))
            }else{
                return(estimateGLMCommonDisp(counts[, batches_ind[[i]]], 
                    design=mod[batches_ind[[i]], ], 
                    subset=nrow(counts)))
            }
        },
        numeric(1))
    
    ## Estimate gene-wise dispersion within each batch 
    genewise_disp_lst <- lapply(seq_len(n_batch), function(j){
        if((n_batches[j] <= ncol(design)-ncol(batchmod)+1) | 
                qr(mod[batches_ind[[j]], ])$rank < ncol(mod)){
            # not enough residual degrees of freedom - use the common dispersion
            return(rep(disp_common[j], nrow(counts)))
        }else{
            return(estimateGLMTagwiseDisp(counts[, batches_ind[[j]]], 
                design = mod[batches_ind[[j]], ], 
                dispersion = disp_common[j], prior.df=0))
        }
    })
    names(genewise_disp_lst) <- paste0('batch', levels(batch))
    
    ## construct dispersion matrix
    phi_matrix <- matrix(NA, nrow=nrow(counts), ncol=ncol(counts))
    for(k in seq_len(n_batch)){
        phi_matrix[, batches_ind[[k]]] <- vec2mat(genewise_disp_lst[[k]],
            n_batches[k]) 
    }
    
    return(list(genewise_disp_lst = genewise_disp_lst, phi_matrix = phi_matrix))
}

########  Estimate parameters from NB GLM  ########
est_parameters_NB_GLM <- function(dge_obj, design, phi_matrix, n_batch, 
    n_batches, n_sample, counts, genewise_disp_lst){
    message("Fitting the GLM model")
    glm_f <- glmFit(dge_obj, 
        design = design, 
        dispersion = phi_matrix, 
        prior.count = 1e-4) # no intercept - nonEstimable;
    # compute offset (library sizes) within function
    
    alpha_g <- glm_f$coefficients[, seq_len(n_batch)] %*% 
        as.matrix(n_batches/n_sample) # compute intercept as batch-size-weighted
    # average from batches
    new_offset <- t(vec2mat(getOffset(dge_obj), 
        nrow(counts))) +   # original offset - sample (library) size
        vec2mat(alpha_g, ncol(counts)) # new offset - gene background expression
    
    # getOffset(dge_obj) is the same as 
    # log(dge_obj$samples$lib.size)
    
    glm_f2 <- glmFit.default(dge_obj$counts,
        design = design, 
        dispersion = phi_matrix, 
        offset = new_offset, 
        prior.count = 1e-4) 
    
    gamma_hat <- glm_f2$coefficients[, seq_len(n_batch)]
    mu_hat <- glm_f2$fitted.values
    phi_hat <- do.call(cbind, genewise_disp_lst)
    
    return(list(gamma_hat = gamma_hat, mu_hat = mu_hat, phi_hat = phi_hat))
}

########  In each batch, compute posterior estimation through   ######## 
########  Monte-Carlo integration  ########  
monte_carlo_est <- function(shrink, counts, batches_ind, mu_hat, gamma_hat, 
    shrink.disp, phi_hat, gene.subset.n, batch, n_batch, n_batches) {
    if(shrink){
        message("Apply shrinkage - computing posterior estimates for ", 
            "parameters")
        monte_carlo_res <- lapply(seq_len(n_batch), function(ii){
            if(ii==1){
                mcres <- monte_carlo_int_NB(dat = counts[, batches_ind[[ii]]], 
                    mu = mu_hat[, batches_ind[[ii]]], 
                    gamma = gamma_hat[, ii],
                    phi = phi_hat[, ii],
                    gene.subset.n = gene.subset.n)
            }else{
                invisible(capture.output(mcres <- monte_carlo_int_NB(
                    dat=counts[, batches_ind[[ii]]],
                    mu=mu_hat[, batches_ind[[ii]]], 
                    gamma=gamma_hat[, ii], 
                    phi=phi_hat[, ii], 
                    gene.subset.n=gene.subset.n)))
            }
            return(mcres)
        })
        names(monte_carlo_res) <- paste0('batch', levels(batch))
        
        gamma_star_mat <- lapply(monte_carlo_res, function(res){res$gamma_star})
        gamma_star_mat <- do.call(cbind, gamma_star_mat)
        phi_star_mat <- lapply(monte_carlo_res, function(res){res$phi_star})
        phi_star_mat <- do.call(cbind, phi_star_mat)
        
        if(!shrink.disp){
            message("Apply shrinkage to mean only")
            phi_star_mat <- phi_hat
        }
    }else{
        message("Shrinkage off - using GLM estimates for parameters")
        gamma_star_mat <- gamma_hat
        phi_star_mat <- phi_hat
    }
    
    ########  Obtain adjusted batch-free distribution  ########
    mu_star <- matrix(NA, nrow=nrow(counts), ncol=ncol(counts))
    for(jj in seq_len(n_batch)){
        mu_star[, batches_ind[[jj]]] <- exp(log(mu_hat[, batches_ind[[jj]]]) - 
                vec2mat(gamma_star_mat[, jj], n_batches[jj]))
    }
    phi_star <- rowMeans(phi_star_mat)
    
    return(list(mu_star = mu_star, phi_star = phi_star))
}

####  Monte Carlo integration functions
monte_carlo_int_NB <- function(dat, mu, gamma, phi, gene.subset.n){
    weights <- pos_res <- list()
    for(i in seq_len(nrow(dat))){
        m <- mu[-i,!is.na(dat[i,])]
        x <- dat[i,!is.na(dat[i,])]
        gamma_sub <- gamma[-i]
        phi_sub <- phi[-i]
        
        # take a subset of genes to do integration - save time
        if(!is.null(gene.subset.n) & is.numeric(gene.subset.n) & 
                length(gene.subset.n)==1){
            if(i==1){
                message(sprintf("Using %s random genes for Monte Carlo ",
                "integration", gene.subset.n))
            }
            mcint_ind <- sample(seq_len((nrow(dat)-1)), gene.subset.n,
                replace=FALSE)
            m <- m[mcint_ind, ];
            gamma_sub <- gamma_sub[mcint_ind]; 
            phi_sub <- phi_sub[mcint_ind]
            G_sub <- gene.subset.n
        }else{
            if(i==1){
                message("Using all genes for Monte Carlo integration; the ",
                "function runs very slow for large number of genes")}
            G_sub <- nrow(dat)-1
        }
        
        LH <- vapply(seq_len(G_sub), 
            function(j){prod(dnbinom(x, mu=m[j,], size=1/phi_sub[j]))},
            numeric(1))
        LH[is.nan(LH)] <- 0; 
        if(sum(LH)==0 | is.na(sum(LH))){
            pos_res[[i]] <- c(gamma.star=as.numeric(gamma[i]),
                phi.star=as.numeric(phi[i]))
        }else{
            pos_res[[i]] <- c(gamma.star=sum(gamma_sub*LH)/sum(LH),
                phi.star=sum(phi_sub*LH)/sum(LH))
        }
        
        weights[[i]] <- as.matrix(LH/sum(LH))
    }
    pos_res <- do.call(rbind, pos_res)
    weights <- do.call(cbind, weights)
    res <- list(gamma_star=pos_res[, "gamma.star"],
        phi_star=pos_res[, "phi.star"], weights=weights)
    return(res)
} 

####  Match quantiles
match_quantiles <- function(counts_sub, old_mu, old_phi, new_mu, new_phi){
    new_counts_sub <- matrix(NA, nrow=nrow(counts_sub), ncol=ncol(counts_sub))
    for(a in seq_len(nrow(counts_sub))){
        for(b in seq_len(ncol(counts_sub))){
            if(counts_sub[a, b] <= 1){
                new_counts_sub[a,b] <- counts_sub[a, b]
            }else{
                tmp_p <- pnbinom(counts_sub[a, b]-1, mu=old_mu[a, b], 
                    size=1/old_phi[a])
                if(abs(tmp_p-1)<1e-4){
                    new_counts_sub[a,b] <- counts_sub[a, b]  
                    # for outlier count, if p==1, will return Inf values -> use 
                    # original count instead
                }else{
                    new_counts_sub[a,b] <- 1+qnbinom(tmp_p, mu=new_mu[a, b],
                        size=1/new_phi[a])
                }
            }
        }
    }
    return(new_counts_sub)
}

mapDisp <- function(old_mu, new_mu, old_phi, divider){
    new_phi <- matrix(NA, nrow=nrow(old_mu), ncol=ncol(old_mu))
    for(a in seq_len(nrow(old_mu))){
        for(b in seq_len(ncol(old_mu))){
            old_var <- old_mu[a, b] + old_mu[a, b]^2 * old_phi[a, b]
            new_var <- old_var / (divider[a, b]^2)
            new_phi[a, b] <- (new_var - new_mu[a, b]) / (new_mu[a, b]^2)
        }
    }
    return(new_phi)
}

adjust_data <- function(counts, n_batch, batches_ind, mu_hat, phi_hat, mu_star,
    phi_star){
    message("Adjusting the data")
    adjust_counts <- matrix(NA, nrow=nrow(counts), ncol=ncol(counts))
    for(kk in seq_len(n_batch)){
        counts_sub <- counts[, batches_ind[[kk]]]
        old_mu <- mu_hat[, batches_ind[[kk]]]
        old_phi <- phi_hat[, kk]
        new_mu <- mu_star[, batches_ind[[kk]]]
        new_phi <- phi_star
        adjust_counts[, batches_ind[[kk]]] <- match_quantiles(
            counts_sub=counts_sub, 
            old_mu=old_mu, 
            old_phi=old_phi, 
            new_mu=new_mu, 
            new_phi=new_phi)
    }
    return(adjust_counts)
}

add_0_counts <- function(countsOri, adjust_counts, keep, rm){
    adjust_counts_whole <- matrix(NA,
        nrow=nrow(countsOri),
        ncol=ncol(countsOri))
    dimnames(adjust_counts_whole) <- dimnames(countsOri)
    adjust_counts_whole[keep, ] <- adjust_counts
    adjust_counts_whole[rm, ] <- countsOri[rm, ]
    return(adjust_counts_whole)
}