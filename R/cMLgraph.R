#' Screen instruments and prepare LD-aware perturbation matrices
#'
#' For each pair of traits, restricts the candidate IV set to SNPs that are
#' GWAS-significant for the exposure (and not more strongly associated with
#' the outcome), and pre-computes the LD-block square-root matrices used to
#' generate joint perturbations of the GWAS estimates. Optionally inflates the
#' standard errors using LDSC-style factors.
#'
#' @param b_mat Numeric matrix of GWAS effect estimates with one column per
#'   trait. Row names must be SNP rsids.
#' @param se_mat Numeric matrix of standard errors matching `b_mat`.
#' @param n_vec Integer vector of GWAS sample sizes (one per trait).
#' @param IV_list List of length `N*(N-1)/2` giving candidate IV rsids for
#'   each trait pair, in order `(1,2), (1,3), ..., (N-1,N)`.
#' @param R_list List of LD blocks. Each element has `R` (an LD correlation
#'   matrix with rsid row/column names) and `snp` (the rsids in the block).
#' @param rho_mat `N x N` correlation matrix between GWAS estimates across
#'   traits (e.g. from bivariate LDSC).
#' @param c_vec Numeric vector of inflation factors (length `N`) applied as
#'   `sqrt(max(1, c))` to each trait's standard errors.
#' @param sig.cutoff Numeric; GWAS-significance threshold for instrument
#'   selection.
#'
#' @return A list with the per-pair screened IV indices (`IJ_snp_list`), the
#'   per-block perturbation square-root matrices (`DP_mat_list`) and the
#'   (possibly inflated) `b_mat` and `se_mat`.
#' @export
#' @importFrom stats pnorm
Graph_Screen <- function(b_mat,se_mat,n_vec,IV_list,R_list,rho_mat,c_vec=rep(1,length(n_vec)),sig.cutoff=5e-08){
  n_trait = length(n_vec)
  N_combination = n_trait * (n_trait - 1) / 2
  if(length(IV_list)!=N_combination){stop("The length of IV_list must be equal to N_combination!")}
  m_block = length(R_list) ## LD blocks
  if(sum(unlist(lapply(R_list,function(x){nrow(x$R)})))!=nrow(b_mat)){stop("LD matrix must contain all SNPs used in the analysis!")}

  IJ_snp_list = vector("list", N_combination)
  k = 1
### Screening for IVs for each pair of traits ###
  ## i:trait 1; j:trait 2
  for(i in 1:(n_trait-1)){
    for(j in (i+1):n_trait){
      SNP_used = IV_list[[k]]
      b_X = b_mat[,i]
      b_Y = b_mat[,j]
      se_X = se_mat[,i]
      se_Y = se_mat[,j]
      n_X = n_vec[i]
      n_Y = n_vec[j]
      keep_ind = intersect(which(!is.na(b_X)),which(!is.na(b_Y)))
      pvalue.X = pnorm(-abs(b_X/se_X))*2
      pvalue.Y = pnorm(-abs(b_Y/se_Y))*2
      cor_X = b_X / sqrt(b_X^2 + (n_X-2)*se_X^2)
      cor_Y = b_Y / sqrt(b_Y^2 + (n_Y-2)*se_Y^2)

      ind_X = which(pvalue.X<(sig.cutoff))
      ind_Y = which(pvalue.Y<(sig.cutoff))
      # With Screening
      intersect.ind.X.Y = intersect(ind_X,ind_Y)
      ind_X_new = setdiff(ind_X,
                          intersect.ind.X.Y[(abs(cor_X)[intersect.ind.X.Y])<
                                              (abs(cor_Y)[intersect.ind.X.Y])])
      ind_Y_new = setdiff(ind_Y,
                          intersect.ind.X.Y[(abs(cor_X)[intersect.ind.X.Y])>
                                              (abs(cor_Y)[intersect.ind.X.Y])])
      ind_X_new = intersect(ind_X_new,keep_ind)
      ind_Y_new = intersect(ind_Y_new,keep_ind)
      ind_X_new_snp = rownames(b_mat)[ind_X_new]
      ind_Y_new_snp = rownames(b_mat)[ind_Y_new]

      ind_X_new1 = ind_X_new[is.element(ind_X_new_snp,SNP_used)]
      ind_Y_new1 = ind_Y_new[is.element(ind_Y_new_snp,SNP_used)]
      IJ_snp_list[[k]]$ind_i_new = ind_X_new1
      IJ_snp_list[[k]]$ind_j_new = ind_Y_new1
      k = k + 1
    }
  }
### End of screening 1 ###
### Screening for LD of b_mat ###
  DP_mat_list = vector("list", m_block)
  for(i in 1:m_block){
      R = R_list[[i]]$R
      match_order = match(rownames(b_mat),rownames(R))
      match_order = match_order[!is.na(match_order)]
      R = R[match_order,match_order,drop=FALSE]
      snp_in_block = R_list[[i]]$snp
      PXR = kronecker(rho_mat,R)
      eigen_decomp = eigen(PXR)
      D = diag(sqrt(zapsmall(eigen_decomp$value,digits=6)))
      DP_mat_list[[i]]$V = eigen_decomp$vector %*% D
      DP_mat_list[[i]]$snp = rownames(R)
  }
### End of screening 2 ###

### Screnning for inflation factor (c) from LDSC ###
  for(i in 1:n_trait){
      se_mat[,i] = max(1,sqrt(c_vec[i])) * se_mat[,i]
  }
### End of screnning 3 ###
  out = list()
  out$IJ_snp_list = IJ_snp_list
  out$DP_mat_list = DP_mat_list
  out$b_mat = b_mat
  out$se_mat = se_mat
  return(out)
}

#' Generate one LD-aware perturbed copy of the GWAS estimates
#'
#' Adds correlated Gaussian noise to `b_mat` using the per-block square-root
#' matrices from `Graph_Screen`, preserving cross-trait correlation and
#' within-block LD.
#'
#' @param b_mat Numeric matrix of GWAS effect estimates (rows = SNPs,
#'   columns = traits).
#' @param se_mat Numeric matrix of standard errors matching `b_mat`.
#' @param n_vec Integer vector of GWAS sample sizes (one per trait).
#' @param rho_mat `N x N` correlation matrix between GWAS estimates across
#'   traits.
#' @param DP_mat_list List of LD-block perturbation matrices produced by
#'   `Graph_Screen`.
#'
#' @return A numeric matrix of the same shape as `b_mat` containing the
#'   perturbed effect estimates.
#' @export
#' @importFrom stats rnorm
Generate_Perturb <- function(b_mat,se_mat,n_vec,rho_mat,DP_mat_list){
  n_trait = length(n_vec)
  m_used = nrow(b_mat)
  e_mat_dp = matrix(0,ncol=n_trait,nrow=m_used)
  for(i in 1:length(DP_mat_list)){
      V = DP_mat_list[[i]]$V
      X = rnorm(nrow(V))
      e_vec = V %*% X
      e_mat = matrix(e_vec,ncol=n_trait,byrow=FALSE)
      snp_ind = which(is.element(rownames(b_mat),DP_mat_list[[i]]$snp))
      e_mat_dp[snp_ind,] = se_mat[snp_ind,,drop=FALSE] * e_mat
  }
  b_mat_dp = b_mat + e_mat_dp

  return(b_mat_dp)
}

#' Estimate the pairwise causal-effect graph
#'
#' For every pair of traits `(i, j)`, runs `mr_cML_O` in both directions using
#' the screened IV indices, fills in the corresponding off-diagonal entries of
#' the observed-effect graph, and applies network deconvolution to obtain the
#' direct-effect graph.
#'
#' @param b_mat Numeric matrix of GWAS effect estimates (rows = SNPs,
#'   columns = traits).
#' @param se_mat Numeric matrix of standard errors matching `b_mat`.
#' @param n_vec Integer vector of GWAS sample sizes (one per trait).
#' @param rho_mat `N x N` correlation matrix between GWAS estimates across
#'   traits.
#' @param IJ_snp_list Per-pair list of screened IV indices, as returned by
#'   `Graph_Screen`.
#' @param t Instrument-selection threshold on the z-scale, passed through to
#'   `mr_cML_O`.
#' @param random_start Integer; number of random starts per cML fit.
#'
#' @return A list with `obs_graph`, `obs_graph_se`, `obs_graph_pval` and the
#'   network-deconvolved `dir_graph`.
#' @export
Graph_Estimate <- function(b_mat,se_mat,n_vec,rho_mat,IJ_snp_list,t,random_start=10){
  n_trait = length(n_vec)
  obs_graph = matrix(0,nrow=n_trait,ncol=n_trait)
  obs_graph_pval = obs_graph_se = matrix(0,nrow=n_trait,ncol=n_trait)
  k = 1
  for(i in 1:(n_trait-1)){
    for(j in (i+1):n_trait){
      ind_i_new = IJ_snp_list[[k]]$ind_i_new
      ind_j_new = IJ_snp_list[[k]]$ind_j_new

      rho_ij = rho_mat[i,j]
      ItoJ_cML_O_res = mr_cML_O(b_exp=b_mat[ind_i_new,i],
                                 b_out=b_mat[ind_i_new,j],
                                 se_exp=se_mat[ind_i_new,i],
                                 se_out=se_mat[ind_i_new,j],
                                 n = min(n_vec[i],n_vec[j]),
                                 rho = rho_ij,
                                 random_start = random_start,t=t)
      JtoI_cML_O_res = mr_cML_O(b_exp=b_mat[ind_j_new,j],
                                     b_out=b_mat[ind_j_new,i],
                                     se_exp=se_mat[ind_j_new,j],
                                     se_out=se_mat[ind_j_new,i],
                                     n = min(n_vec[i],n_vec[j]),
                                     rho = rho_ij,
                                     random_start = random_start,t=t)

     obs_graph[i,j] = ItoJ_cML_O_res$BIC_theta
     obs_graph[j,i] = JtoI_cML_O_res$BIC_theta
     obs_graph_se[i,j] = ItoJ_cML_O_res$BIC_se
     obs_graph_se[j,i] = JtoI_cML_O_res$BIC_se
     obs_graph_pval[i,j] = ItoJ_cML_O_res$BIC_p
     obs_graph_pval[j,i] = JtoI_cML_O_res$BIC_p
      k = k + 1

    }
  }
  dir_graph = obs_graph %*% solve(diag(n_trait)+obs_graph)

  out = list()
  out$obs_graph = obs_graph
  out$obs_graph_se = obs_graph_se
  out$obs_graph_pval = obs_graph_pval
  out$dir_graph = dir_graph

  return(out)

}

#' GraphMRcML with data perturbation
#'
#' Top-level routine: screens IVs once via `Graph_Screen`, then repeatedly
#' generates LD-aware perturbed datasets and fits the pairwise causal-effect
#' graph on each one. Perturbation replicates are run in parallel via
#' [pbmcapply::pbmclapply()].
#'
#' @param b_mat Numeric matrix of GWAS effect estimates (rows = SNPs,
#'   columns = traits). Row names must be SNP rsids.
#' @param se_mat Numeric matrix of standard errors matching `b_mat`.
#' @param n_vec Integer vector of GWAS sample sizes (one per trait).
#' @param rho_mat `N x N` correlation matrix between GWAS estimates across
#'   traits.
#' @param IV_list List of length `N*(N-1)/2` giving candidate IV rsids for
#'   each trait pair.
#' @param R_list List of LD blocks (see `Graph_Screen`).
#' @param c_vec Numeric vector of LDSC inflation factors (length `N`).
#' @param sig.cutoff Numeric; GWAS-significance threshold for instrument
#'   selection.
#' @param num_pert Integer; number of data-perturbation replicates.
#' @param random_start Integer; number of random starts per cML fit.
#' @param seed Integer; RNG seed.
#' @param trait_vec Optional character vector of trait names; defaults to
#'   `colnames(b_mat)`.
#' @param curse Logical; if `TRUE`, applies the "curse-of-winner" correction
#'   by setting the cML threshold to `-qnorm(sig.cutoff/2)`.
#'
#' @return A list of per-perturbation observed graphs, their standard errors
#'   and p-values, the deconvolved direct-effect graphs, and `trait_vec`.
#' @export
#' @importFrom stats qnorm
Graph_Perturb <- function(b_mat,se_mat,n_vec,rho_mat,IV_list,R_list,c_vec=rep(1,length(n_vec)),
                          sig.cutoff=5e-08,num_pert=100,random_start=10,seed=0,trait_vec=NULL,curse=F){
    set.seed(seed)
    if(curse){
        t = -qnorm(sig.cutoff/2)
    }else{t=0}

  screen_res = Graph_Screen(b_mat=b_mat,se_mat=se_mat,
                            n_vec=n_vec,IV_list=IV_list,R_list=R_list,rho_mat=rho_mat,c_vec=c_vec,
                            sig.cutoff=sig.cutoff)
  if(is.null(trait_vec)){trait_vec=colnames(screen_res$b_mat)}
  out_list = pbmcapply::pbmclapply(1:num_pert,function(i){    b_mat_dp = Generate_Perturb(b_mat=screen_res$b_mat,
                                se_mat=screen_res$se_mat,
                                n_vec=n_vec,rho_mat=rho_mat,DP_mat_list=screen_res$DP_mat_list);
                                Graph_Estimate(b_mat=b_mat_dp,
                                     se_mat=screen_res$se_mat,
                                     n_vec=n_vec,rho_mat=rho_mat,
                                     IJ_snp_list=screen_res$IJ_snp_list,
                                     t=t,
                                     random_start=random_start)},ignore.interactive=TRUE)
  obs_graph_list = lapply(out_list,function(x){x$obs_graph})
  obs_graph_se_list = lapply(out_list,function(x){x$obs_graph_se})
  obs_graph_pval_list = lapply(out_list,function(x){x$obs_graph_pval})
  dir_graph_list = lapply(out_list,function(x){x$dir_graph})

  out = list()
  out$obs_graph_list = obs_graph_list
  out$obs_graph_se_list = obs_graph_se_list
  out$obs_graph_pval_list = obs_graph_pval_list
  out$dir_graph_list = dir_graph_list
  out$trait_vec = trait_vec
  return(out)
}






