#' Adjust for batch effects using an empirical Bayes framework
#'
#' ComBat allows users to adjust for batch effects in datasets where the batch
#' covariate is known, using methodology described in Johnson et al. 2007. It
#' uses either parametric or non-parametric empirical Bayes frameworks for
#' adjusting data for batch effects.  Users are returned an expression matrix
#' that has been corrected for batch effects. The input data are assumed to be
#' cleaned and normalized before batch effect removal.
#'
#' @param dat Genomic measure matrix (dimensions probe x sample) - for example,
#'   expression matrix
#' @param batch {Batch covariate (only one batch allowed)}
#' @param mod Model matrix for outcome of interest and other covariates besides
#'   batch
#' @param par.prior (Optional) TRUE indicates parametric adjustments will be
#'   used, FALSE indicates non-parametric adjustments will be used
#' @param prior.plots (Optional) TRUE give prior plots with black as a kernel
#'   estimate of the empirical batch effect density and red as the parametric
#' @param mean.only (Optional) FALSE If TRUE ComBat only corrects the mean of
#'   the batch effect (no scale adjustment)
#' @param ref.batch (Optional) NULL If given, will use the selected batch as a
#'   reference for batch adjustment.
#' @param BPPARAM (Optional) BiocParallelParam for parallel operation
#'
#' @return data A probe x sample genomic measure matrix, adjusted for batch
#'   effects.
#'
#' @importFrom graphics lines par
#' @importFrom stats cor density dnorm model.matrix pf ppoints prcomp predict
#' qgamma qnorm qqline qqnorm qqplot smooth.spline var
#' @importFrom utils read.delim
#'
#' @examples
#' library(bladderbatch)
#' data(bladderdata)
#' dat <- bladderEset[1:50,]
#'
#' pheno <- pData(dat)
#' edata <- exprs(dat)
#' batch <- pheno$batch
#' mod <- model.matrix(~as.factor(cancer), data=pheno)
#'
#' # parametric adjustment
#' combat_edata1 <- ComBat(dat=edata, batch=batch, mod=NULL, par.prior=TRUE, 
#'   prior.plots=FALSE)
#'
#' # non-parametric adjustment, mean-only version
#' combat_edata2 <- ComBat(dat=edata, batch=batch, mod=NULL, par.prior=FALSE,
#'   mean.only=TRUE)
#'
#' # reference-batch version, with covariates
#' combat_edata3 <- ComBat(dat=edata, batch=batch, mod=mod, par.prior=TRUE,
#'   ref.batch=3)
#'
#' @export
#'

ComBat <- function(dat, batch, mod = NULL, par.prior = TRUE,
    prior.plots = FALSE, mean.only = FALSE, ref.batch = NULL,
    BPPARAM = bpparam("SerialParam")) {
    
    # Check for only one batch and coerce dat into a matrix
    one_batch_check(batch)
    dat <- as.matrix(dat)
    
    ## find genes with zero variance in any of the batches
    zv_res <- zero_variance(dat, batch)
    
    ## Make a set of batch indicators/characteristics and create design matrix
    bd_res <- batch_design(zv_res$batch, mean.only, ref.batch, mod)
    
    ## Check if the design is confounded and for missing values
    confounding_check(bd_res$design, bd_res$n.batch, test_name = "ComBat")
    NAs <- NA_check(zv_res$dat)
    
    ## Standardize Data across genes
    B.hat <- standardize_data(NAs, bd_res$design, zv_res$dat)
    
    ## change grand.mean and var.pooled for ref batch
    ref_batch_vars <- ref_batch_variables_update(ref.batch, bd_res$ref, B.hat,
        bd_res$n.batches, bd_res$n.array, bd_res$n.batch, NAs, bd_res$batches,
        zv_res$dat, bd_res$design)
    
    ##Get regression batch effect parameters
    rx_batch_para <- regression_batch_parameters(bd_res$design, bd_res$n.batch,
        NAs, ref_batch_vars$s.data, bd_res$batches, bd_res$mean.only)
    
    ## Find priors; Plot empirical and parametric priors
    priors_res <- priors_calculations(rx_batch_para$gamma.hat,
        rx_batch_para$delta.hat, prior.plots, par.prior)
    
    ## Find EB batch adjustments
    EB_batch_adj_res <- EB_batch_adjustment(bd_res$n.batch,
        ref_batch_vars$s.data, par.prior, bd_res$mean.only,
        rx_batch_para$gamma.hat, priors_res$gamma.bar, priors_res$t2,
        bd_res$batches, rx_batch_para$delta.hat, priors_res$a.prior,
        priors_res$b.prior, bd_res$ref, ref.batch, BPPARAM)
    
    ## Normalize the Data ###
    bayesdata <- normalize_data(ref_batch_vars$s.data, bd_res$batches,
        rx_batch_para$batch.design, EB_batch_adj_res$gamma.star, 
        EB_batch_adj_res$delta.star, bd_res$n.batches,ref_batch_vars$var.pooled,
        bd_res$n.array, ref_batch_vars$stand.mean, ref.batch, bd_res$ref,
        zv_res$dat, zv_res$zero.rows, zv_res$keep.rows, zv_res$dat.orig)
    
    return(bayesdata)
}