#' Adjust for batch effects w/an empirical Bayes framework in RNA-seq raw counts
#' 
#' ComBat_seq is an improved model from ComBat using negative binomial
#' regression, which specifically targets RNA-Seq count data.
#' 
#' @param counts Raw count matrix from genomic studies (dimensions gene x
#'   sample) 
#' @param batch Vector / factor for batch
#' @param group Vector / factor for biological condition of interest 
#' @param covar_mod Model matrix for multiple covariates to include in linear
#'   model (signals from these variables are kept in data after adjustment) 
#' @param full_mod Boolean, if TRUE include condition of interest in model
#' @param shrink Boolean, whether to apply shrinkage on parameter estimation
#' @param shrink.disp Boolean, whether to apply shrinkage on dispersion
#' @param gene.subset.n Number of genes to use in empirical Bayes estimation,
#'   only useful when shrink = TRUE
#' 
#' @return data A gene x sample count matrix, adjusted for batch effects.
#' 
#' @importFrom edgeR DGEList estimateGLMCommonDisp estimateGLMTagwiseDisp glmFit
#'   glmFit.default getOffset
#' @importFrom stats dnbinom lm pnbinom qnbinom
#' @importFrom utils capture.output
#' 
#' @examples 
#' 
#' count_matrix <- matrix(rnbinom(400, size=10, prob=0.1), nrow=50, ncol=8)
#' batch <- c(rep(1, 4), rep(2, 4))
#' group <- rep(c(0,1), 4)
#' 
#' # include condition (group variable)
#' adjusted_counts <- ComBat_seq(count_matrix,
#'   batch=batch,
#'   group=group, 
#'   full_mod=TRUE)
#' 
#' # do not include condition
#' adjusted_counts <- ComBat_seq(count_matrix,
#'   batch=batch,
#'   group=NULL,
#'   full_mod=FALSE)
#' 
#' @export
#' 

ComBat_seq <- function(counts, batch, group=NULL, covar_mod=NULL, full_mod=TRUE,
    shrink=FALSE, shrink.disp=FALSE, gene.subset.n=NULL){  
    
    ########  Preparation  ########  
    prep <- preparation(batch, counts)
    
    ## Make design matrix 
    design_results <- create_design_matrix(prep$batch, group, prep$counts, 
        covar_mod, full_mod)
    
    ######## Check if the design is confounded or for missing values in counts
    confounding_check(design_results$design, prep$n_batch,
        test_name = "ComBat-Seq")
    NAs <- NA_check(prep$counts)
    
    ########  Estimate gene-wise dispersions within each batch  ########
    gene_wise_results <- estimate_gene_wise_dispersions(prep$n_batch, 
        prep$n_batches, design_results$design, design_results$batchmod,
        design_results$mod, prep$batches_ind, prep$counts, prep$batch)
    
    ########  Estimate parameters from NB GLM  ########
    para_NB_GLM_results <- est_parameters_NB_GLM(prep$dge_obj, 
        design_results$design, gene_wise_results$phi_matrix, prep$n_batch,
        prep$n_batches,prep$n_sample, prep$counts,
        gene_wise_results$genewise_disp_lst)
    
    ########  In each batch, compute posterior est. w/Monte-Carlo int. ######## 
    monte_carlo_results <- monte_carlo_est(shrink, prep$counts,
        prep$batches_ind, para_NB_GLM_results$mu_hat,
        para_NB_GLM_results$gamma_hat, shrink.disp, para_NB_GLM_results$phi_hat,
        gene.subset.n, prep$batch, prep$n_batch, prep$n_batches)
    
    ########  Adjust the data  ########  
    adjust_counts <- adjust_data(prep$counts, prep$n_batch, prep$batches_ind, 
        para_NB_GLM_results$mu_hat, para_NB_GLM_results$phi_ha,
        monte_carlo_results$mu_star, monte_carlo_results$phi_star)
    
    ## Add back genes with only 0 counts in any batch (so dims won't change)
    adjust_counts_whole <- add_0_counts(prep$countsOri, adjust_counts,
        prep$keep, prep$rm)
    
    return(adjust_counts_whole)
}
