#include <RcppArmadillo.h>
#include "progress.hpp"

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat weighted_cov(const arma::mat& x, const arma::vec& w) {
  unsigned int p = x.n_cols;

  arma::mat out(p, p, arma::fill::zeros);

  for (unsigned int k = 0; k < p; k++) {
    for (unsigned int j = k; j < p; j++) {
      out(j, k) = arma::accu(x.col(j) % x.col(k) % w);
    }
  }

  out = arma::symmatl(out);
  return out;
}

// [[Rcpp::export]]
arma::mat weighted_reg_beta(const arma::mat& X, const arma::vec& wt, const arma::mat& Y) {
  arma::mat XTWX = weighted_cov(X, wt);
  arma::mat out = arma::solve(XTWX, X.t() * (Y.each_col() % wt));
  return(out);
}

// [[Rcpp::export]]
arma::mat arma_direct2(const arma::mat& X, const arma::vec& w, const arma::mat& y) {
  arma::mat p = solve((X.each_col() % w).t() * X, X.t() * (y.each_col() % w));
  return p;
}

// [[Rcpp::depends(RcppArmadillo)]]

// # B.dim = dim(state$beta.g.h)[1:2] - c(1L, 0L)
// # state$sigma.bar.g.h.inv = lapply(1:state$G, function(k) {
// #   Sginv = solve(state$sigma.g[, , setx[k]])
// #   Sghinv = solve(state$sigma.g.h[, , sety[k]])
// #   B = state$beta.g.h[-1, , sety[k]]
// #   if (is.vector(B)) dim(B) = B.dim
// #   neg.BSghinv = -B %*% Sghinv
// #
// #   tmp = matrix(0, p, p)
// #   tmp[1:px, 1:px] = Sginv - tcrossprod(neg.BSghinv, B)
// #   tmp[1:px, px + 1:py] = neg.BSghinv
// #   tmp[px + 1:py, 1:px] = t(neg.BSghinv)
// #   tmp[px + 1:py, px + 1:py] = Sghinv
// #
// #   if (p > px + py) tmp[(px + py + 1):p, (px + py + 1):p] = diag(1/state$psi, p - px - py)
// #   return(tmp)
// # })

// [[Rcpp::export]]
arma::cube update_bars_sigma_inv(
    const unsigned int G,
    const arma::uvec& setx, const arma::uvec& sety,
    const arma::cube& Sg, const arma::cube& Sgh,
    const arma::cube& B, const arma::vec& psi,
    const unsigned int px, const unsigned int py, const unsigned int p
) {
  arma::cube out(p, p, G, arma::fill::none);

  arma::mat tmplt(p, p, arma::fill::zeros);
  if (px + py < p) {
    tmplt.submat(px + py, px + py, p - 1, p - 1) = arma::diagmat(1 / psi);
  }

  for (unsigned int k = 0; k < G; k++) {
    arma::mat tmp = tmplt;
    arma::mat Sginv = Sg.slice(setx(k) - 1).i();
    arma::mat Sghinv = Sgh.slice(sety(k) - 1).i();
    arma::mat Bslope = B.slice(sety(k) - 1);
    Bslope.shed_row(0);
    arma::mat negBSghinv = -Bslope * Sghinv;

    tmp.submat(0, 0, px - 1, px - 1) = Sginv - negBSghinv * Bslope.t();
    tmp.submat(0, px, px - 1, px + py - 1) = negBSghinv;
    tmp.submat(px, 0, px + py - 1, px - 1) = negBSghinv.t();
    tmp.submat(px, px, px + py - 1, px + py - 1) = Sghinv;

    out.slice(k) = tmp;
  }
  return(out);
}

// # A = rowSums(sapply(1:state$G, function(k) {
// #   wt.mean = rowSums(t(state$data * wt[, k] / n))
// #   muxyz = tcrossprod(wt.mean, state$mu.bar.g.h[[k]])
// #   return(tcrossprod(state$sigma.bar.g.h.inv[[k]], muxyz))
// # }))
// # dim(A) = c(p, p)
// # B = rowSums(vapply(1:ncol(wt), function(k) {
// #   invSxyz = state$sigma.bar.g.h.inv[[k]]
// #
// #   Wk = crossprod(data * wt[, k], data) / n
// #   lamk   = max(eigen(Wk, symmetric = TRUE, only.values = TRUE)$values)
// #
// #   return(tcrossprod(invSxyz, state$Gamma) %*% (Wk - diag(lamk, p)))
// # }, numeric(p * p)))
// # dim(B) = c(p, p)

// [[Rcpp::export]]
arma::mat m_step_gamma(
    const unsigned int G,
    Rcpp::List& sigma_bar_g_h_inv,
    Rcpp::List& mu_bar_g_h,
    const arma::mat& wt,
    const arma::mat& data,
    const arma::mat& Gamma,
    const unsigned int max_update = 10,
    const double epsilon = 1e-8
) {
  unsigned int p = data.n_cols;
  arma::mat wtmean_k = data.t() * wt;
  arma::cube Wk_k(p, p, G, arma::fill::none), invSxyz(p, p, G, arma::fill::none);
  arma::vec max_eigvals(G, arma::fill::none);
  arma::mat mubar(p, G, arma::fill::none);
  arma::mat kiersA(p, p, arma::fill::zeros);

  for (unsigned int k = 0; k < G; k++) {
    Wk_k.slice(k) = weighted_cov(data, wt.col(k));
    invSxyz.slice(k) = as<arma::mat>(sigma_bar_g_h_inv[k]);
    mubar.col(k) = as<arma::vec>(mu_bar_g_h[k]);
    max_eigvals(k) = arma::eig_sym(Wk_k.slice(k)).max();

    kiersA += invSxyz.slice(k) * mubar.col(k) * wtmean_k.col(k).t();
  }

  arma::mat oldgamma = Gamma, newgamma(p, p, arma::fill::none),
    P(p, p, arma::fill::none),
    Q(p, p, arma::fill::none),
    kiersF(p, p, arma::fill::none);
  arma::vec d(p);

  for (unsigned int upd_iter = 0; upd_iter < max_update; upd_iter++) {
    kiersF = kiersA;
    for (unsigned int k = 0; k < G; k++) {
      kiersF -= invSxyz.slice(k) * oldgamma.t() *
        (Wk_k.slice(k) - max_eigvals(k) * arma::eye(p, p));
    }
    svd(P, d, Q, kiersF);
    newgamma = Q * P.t();
    if (arma::accu(arma::abs(newgamma - oldgamma)) / (p * p) < epsilon) {
      // Rprintf("Early Break at %i", upd_iter);
      break;
    }
    oldgamma = newgamma;
  }
  return(newgamma);
}

// [[Rcpp::export]]
arma::mat m_step_gamma_tr(
    const unsigned int G,
    Rcpp::List& sigma_bar_g_h_inv,
    Rcpp::List& mu_bar_g_h,
    const arma::mat& wt,
    const arma::mat& data,
    const arma::mat& Gamma,
    const unsigned int max_update = 10,
    const double epsilon = 1e-8
) {
  unsigned int p = data.n_cols;
  arma::mat wtmean_k = data.t() * wt;
  arma::cube Wk(p, p, G, arma::fill::none), Sinv(p, p, G, arma::fill::none);
  arma::mat kiersA(p, p, arma::fill::zeros);

  for (unsigned int k = 0; k < G; k++) {
    Wk.slice(k) = weighted_cov(data, wt.col(k));
    Sinv.slice(k) = as<arma::mat>(sigma_bar_g_h_inv[k]);
    kiersA += wtmean_k.col(k) * as<arma::rowvec>(mu_bar_g_h[k]) * Sinv.slice(k);
    // Pre-subtract lambda
    Sinv.slice(k) -= arma::eig_sym(Sinv.slice(k)).max() * arma::eye(p, p);
  }

  arma::mat oldgamma = Gamma.t(), newgamma(p, p, arma::fill::none),
    P(p, p, arma::fill::none),
    Q(p, p, arma::fill::none),
    kiersF(p, p, arma::fill::none);
  arma::vec d(p);

  for (unsigned int upd_iter = 0; upd_iter < max_update; upd_iter++) {
    kiersF = kiersA;
    for (unsigned int k = 0; k < G; k++) kiersF -= Wk.slice(k) * oldgamma.t() * Sinv.slice(k);

    svd(P, d, Q, kiersF);
    newgamma = Q * P.t();
    if (arma::accu(arma::abs(newgamma - oldgamma)) / (p * p) < epsilon) {
      // Rprintf("Early Break at %i", upd_iter);
      break;
    }
    oldgamma = newgamma;
  }
  return(newgamma.t());
}


// [[Rcpp::export]]
arma::mat m_step_gamma_dual(
    const unsigned int G,
    Rcpp::List& sigma_bar_g_h_inv,
    Rcpp::List& mu_bar_g_h,
    const arma::mat& wt,
    const arma::mat& data,
    const arma::mat& Gamma,
    const unsigned int max_update = 10,
    const double epsilon = 1e-8
) {
  if (max_update == 0) {return(Gamma);}
  unsigned int p = data.n_cols;
  arma::mat wtmean_k = data.t() * wt;
  arma::cube
    Wk(p, p, G, arma::fill::none), Sinv(p, p, G, arma::fill::none),
    Wk_tr(p, p, G, arma::fill::none), Sinv_tr(p, p, G, arma::fill::none);
  arma::mat kiersA(p, p, arma::fill::zeros), kiersA_tr(p, p, arma::fill::zeros);

  for (unsigned int k = 0; k < G; k++) {
    Wk.slice(k)      = weighted_cov(data, wt.col(k));
    Wk_tr.slice(k)   = Wk.slice(k);
    Sinv.slice(k)    = as<arma::mat>(sigma_bar_g_h_inv[k]);
    Sinv_tr.slice(k) = Sinv.slice(k);
    kiersA_tr       += wtmean_k.col(k) * as<arma::rowvec>(mu_bar_g_h[k]) * Sinv.slice(k);
    // Pre-subtract lambda
    Sinv_tr.slice(k) -= arma::eig_sym(Sinv_tr.slice(k)).max() * arma::eye(p, p);
    Wk.slice(k) -= arma::eig_sym(Wk.slice(k)).max() * arma::eye(p, p);
  }
  kiersA = kiersA_tr.t();

  arma::mat oldgamma = Gamma.t(),
    newgamma(p, p, arma::fill::none),
    P(p, p, arma::fill::none),
    Q(p, p, arma::fill::none),
    kiersF(p, p, arma::fill::none);
  arma::vec d(p);

  for (unsigned int upd_iter = 0; upd_iter < max_update; upd_iter++) {
    // transposed
    kiersF = kiersA_tr;
    for (unsigned int k = 0; k < G; k++) kiersF -= Wk_tr.slice(k) * oldgamma.t() * Sinv_tr.slice(k);
    svd(P, d, Q, kiersF);
    newgamma = P * Q.t();

    // untransposed
    kiersF = kiersA;
    for (unsigned int k = 0; k < G; k++) kiersF -= Sinv.slice(k) * newgamma.t() * Wk.slice(k);
    svd(P, d, Q, kiersF);
    newgamma = Q * P.t();
    if (arma::accu(arma::abs(newgamma - oldgamma)) / (p * p) < epsilon) {
      return(newgamma);
    }
    oldgamma = newgamma.t();
  }
  return(newgamma);
}

// mixlogden = function(state) {
//   pxy = state$px + state$py
//   rot.data = cbind(state$X, state$Y)
//
//   zlog = vapply(1:state$G, function(k) {
//     rot.data.center = t(t(rot.data) - state$mu.bar.g.h[[k]][1:pxy])
//     siginv = state$sigma.bar.g.h.inv[[k]][1:pxy, 1:pxy]
//     mal = rowSums(rot.data.center %*% siginv * rot.data.center)
//     -0.5 * mal + 0.5 * log(det(siginv)) + log(state$tau.g.h[k])
// # constant -pxy/2 * log(2*pi) omitted
//   }, numeric(nrow(rot.data)))
//
//   return(zlog)
// }

// [[Rcpp::export]]
arma::mat mixlogden_helper(
    const unsigned int px, const unsigned int py, const unsigned int G,
    Rcpp::List& sigma_bar_g_h_inv, Rcpp::List& mu_bar_g_h,
    arma::vec& tau_g_h,
    arma::mat& rotX, arma::mat& rotY
) {
  unsigned int pxy = px + py;

  arma::mat rotdata = arma::join_horiz(rotX, rotY);
  arma::mat out(rotdata.n_rows, G, arma::fill::none);

  for (unsigned int k = 0; k < G; k++) {
    arma::mat siginv = sigma_bar_g_h_inv[k];
    siginv = siginv.submat(0, 0, pxy - 1, pxy - 1);
    arma::vec mubar = mu_bar_g_h[k];
    arma::mat rotdatacenter = rotdata.each_row() - mubar.head(pxy).t();
    arma::vec mal = arma::sum(rotdatacenter * siginv % rotdatacenter, 1);
    out.col(k) = -0.5 * mal + 0.5 * log(arma::det(siginv)) + log(tau_g_h(k));
  }

  return(out);
}

// logden = function(state) {
//   p = nrow(state$Gamma)
//   rot.data = cbind(state$X, state$Y, state$Z)
//
//   zlog = vapply(1:state$G, function(k) {
//     rot.data.center = t(t(rot.data) - state$mu.bar.g.h[[k]])
//     siginv = state$sigma.bar.g.h.inv[[k]]
//     mal = rowSums(rot.data.center %*% siginv * rot.data.center)
//     -0.5 * mal + 0.5 * log(det(siginv)) + log(state$tau.g.h[k])
//   }, numeric(nrow(rot.data)))
//
//   return(zlog - p/2*log(2*pi))
// }

// [[Rcpp::export]]
arma::mat logden_helper(
    const unsigned int G,
    Rcpp::List& sigma_bar_g_h_inv, Rcpp::List& mu_bar_g_h,
    arma::vec& tau_g_h,
    arma::mat& rotX, arma::mat& rotY, arma::mat& rotZ,
    const arma::mat& data
) {
  arma::mat out(data.n_rows, G, arma::fill::none);
  arma::mat rotdata = arma::join_horiz(rotX, rotY, rotZ);

  for (unsigned int k = 0; k < G; k++) {
    arma::mat siginv = sigma_bar_g_h_inv[k];
    arma::vec mubar = mu_bar_g_h[k];
    arma::mat rotdatacenter = rotdata.each_row() - mubar.t();
    arma::vec mal = arma::sum(rotdatacenter * siginv % rotdatacenter, 1);
    out.col(k) = -0.5 * mal + 0.5 * log(arma::det(siginv)) + log(tau_g_h(k));
  }

  out -= data.n_cols * 0.91893853320467267;

  return(out);
}

// logden.wts    <- function(z = NULL) {
//   expz = exp(z)
//   expz / rowSums(expz)
// }

// [[Rcpp::export]]
arma::mat logden_wts(arma::mat& z) {
  z.each_col() -= arma::max(z, 1);
  arma::mat out = arma::normalise(arma::exp(z), 1, 1);
  return(out);
}

// logden.loglik <- function(z = NULL) {
//   log(rowSums(exp(z)))
// }

// [[Rcpp::export]]
arma::vec logden_loglik(arma::mat& z) {
  arma::vec out = arma::log(arma::sum(arma::exp(z), 1));
  return(out);
}

// for (g in 1:state$Gx) {
//   wt.g = rowsums(z_ngh[, setx == g, drop = FALSE])
//   pi.g[g] = sum(wt.g)
//   state$mu.g[,g] = colsums(X * wt.g) / sum(wt.g)
//   XR = t(X) - state$mu.g[,g]
//   X.cov = XR %*% (t(XR) * wt.g) / sum(wt.g)
//   state$sigma.g[,,g] = X.cov
// }

// [[Rcpp::export]]
Rcpp::List m_step_primary(
    const arma::mat& z_ngh, const arma::uvec& setx,
    const unsigned int Gx, const arma::mat& X,
    const bool diagonal = false
){
  arma::vec pi_g(Gx);
  arma::mat mu_g(X.n_cols, Gx);
  arma::cube sigma_g(X.n_cols, X.n_cols, Gx);

  for (unsigned int g = 0; g < Gx; g++) {
    arma::vec wt_g = arma::sum(z_ngh.cols(arma::find(setx == (g + 1))), 1);
    pi_g(g) = arma::accu(wt_g);
    wt_g /= pi_g(g);
    mu_g.col(g) = arma::trans(arma::sum(X.each_col() % wt_g));
    arma::mat XR = X.each_row() - mu_g.col(g).t();
    if (diagonal) {
      XR = arma::square(XR);
      sigma_g.slice(g) = arma::diagmat(arma::sum(XR.each_col() % wt_g).t());
    } else {
      sigma_g.slice(g) = weighted_cov(XR, wt_g);
    }
  }

  return Rcpp::List::create(Named("pi.g") = pi_g,
                            Named("mu.g") = mu_g,
                            Named("sigma.g") = sigma_g);
}

// nh = if (state$subcluster) state$G else state$Gy[1]
// pi.g.h = numeric(nh)
//   for (h in 1:nh) {
//     h.indices = sety == h
//     wt.h = rowsums(z_ngh[, h.indices, drop = FALSE])
//     pi.g.h[h] = sum(wt.h)
//     X1.wt = X1 * wt.h
//     if (state$regression) {
//       state$beta.g.h[, , h.indices] =
//         tryCatch(solve(crossprod(X1.wt, X1), crossprod(X1.wt, Y)),
//                  error = function(e)
//                    MASS::ginv(crossprod(X1.wt, X1)) %*% crossprod(X1.wt, Y))
// # state$beta.g.h[, , h.indices] = weighted_reg_beta(X1, wt.h, Y)
//       YR = t(Y - X1 %*% state$beta.g.h[, , h])
//     } else {
//       state$beta.g.h[1, , h.indices] = colsums(Y * wt.h) / sum(wt.h)
//       YR = t(Y) - state$beta.g.h[1, , h]
//     }
//
// # Y.cov = YR %*% (t(YR) * wt.h) / pi.g.h[h]
//     Y.cov = weighted_cov(Y, wt.h) / pi.g.h[h]
//     state$sigma.g.h[, , h.indices] = Y.cov
//   }

// [[Rcpp::export]]
Rcpp::List m_step_secondary(
    const arma::mat& z_ngh, const arma::uvec& sety, const arma::uvec& Gy,
    const arma::mat& X, const arma::mat& Y,
    const bool regression, const bool subcluster,
    const bool diagonal = false
){
  unsigned int nh;
  if (subcluster) {nh = arma::accu(Gy);} else {nh = Gy(0);}

  arma::vec pi_g_h(nh, arma::fill::none);
  arma::cube beta_g_h(X.n_cols + 1, Y.n_cols, nh, arma::fill::none);
  arma::cube sigma_g_h(Y.n_cols, Y.n_cols, nh, arma::fill::none);
  arma::mat intercept(X.n_rows, 1, arma::fill::ones);
  arma::mat X1 = arma::join_horiz(intercept, X);

  for (unsigned int h = 0; h < nh; h++) {
    arma::vec wt_h = arma::sum(z_ngh.cols(arma::find(sety == h + 1)), 1);
    pi_g_h(h) = arma::accu(wt_h);
    wt_h /= pi_g_h(h);
    arma::mat YR;
    if (regression) {
      beta_g_h.slice(h) = weighted_reg_beta(X1, wt_h, Y);
      YR = Y - (X1 * beta_g_h.slice(h));
    } else {
      beta_g_h.slice(h).fill(0);
      beta_g_h.slice(h).row(0) = arma::sum(Y.each_col() % wt_h);
      YR = Y.each_row() - beta_g_h.slice(h).row(0);
    }
    if (diagonal) {
      YR = arma::square(YR);
      sigma_g_h.slice(h) = arma::diagmat(arma::sum(YR.each_col() % wt_h).t());
    } else {
      sigma_g_h.slice(h) = weighted_cov(YR, wt_h);
    }
  }

  return Rcpp::List::create(Named("pi.g.h") = pi_g_h,
                            Named("beta.g.h") = beta_g_h,
                            Named("sigma.g.h") = sigma_g_h);
}

// [[Rcpp::export]]
arma::mat e_step_helper(
    const unsigned int px, const unsigned int py, const unsigned int G,
    Rcpp::List& sigma_bar_g_h_inv, Rcpp::List& mu_bar_g_h,
    arma::vec& tau_g_h, arma::mat& rotX, arma::mat& rotY
){
  arma::mat zlog = mixlogden_helper(px, py, G, sigma_bar_g_h_inv, mu_bar_g_h, tau_g_h, rotX, rotY);
  zlog = logden_wts(zlog);
  return(zlog);
}

// In-place Version -----------------

void m_step_primary_inplace(
    arma::vec& pi_g, arma::mat& mu_g, arma::cube& sigma_g,
    const arma::mat& z_ngh, const arma::uvec& setx,
    const unsigned int Gx, const arma::mat& X, const bool diagonal
){
  arma::vec wt_g;
  arma::mat XR;
  for (unsigned int g = 0; g < Gx; g++) {
    wt_g = arma::sum(z_ngh.cols(arma::find(setx == (g + 1))), 1);
    pi_g(g) = arma::accu(wt_g);
    wt_g /= pi_g(g);
    double max_wt_g = arma::max(wt_g);
    wt_g /= max_wt_g;
    mu_g.col(g) = arma::trans(arma::sum(X.each_col() % wt_g)) * max_wt_g;
    XR = X.each_row() - mu_g.col(g).t();
    if (diagonal) {
      XR = arma::square(XR);
      sigma_g.slice(g) = arma::diagmat(arma::sum(XR.each_col() % wt_g).t()) * max_wt_g;
    } else {
      sigma_g.slice(g) = weighted_cov(XR, wt_g) * max_wt_g;
    }
  }
  return;
}

void m_step_secondary_inplace(
    arma::vec& pi_g_h, arma::cube& beta_g_h, arma::cube& sigma_g_h,
    const arma::mat& z_ngh, const arma::uvec& sety, const arma::uvec& Gy,
    const arma::mat& X, const arma::mat& Y,
    const bool regression, const bool subcluster, const bool diagonal
){
  unsigned int nh;
  if (subcluster) {nh = arma::accu(Gy);} else {nh = Gy(0);}

  arma::mat intercept(X.n_rows, 1, arma::fill::ones);
  arma::mat X1 = arma::join_horiz(intercept, X);
  arma::vec wt_h;
  arma::mat YR;

  if (!regression) {beta_g_h.fill(0);}
  for (unsigned int h = 0; h < nh; h++) {
    wt_h = arma::sum(z_ngh.cols(arma::find(sety == h + 1)), 1);
    pi_g_h(h) = arma::accu(wt_h);
    wt_h /= pi_g_h(h);
    double max_wt_h = arma::max(wt_h);
    wt_h /= max_wt_h;
    if (regression) {
      beta_g_h.slice(h) = weighted_reg_beta(X1, wt_h, Y) * max_wt_h;
      YR = Y - (X1 * beta_g_h.slice(h));
    } else {
      beta_g_h.slice(h).row(0) = arma::sum(Y.each_col() % wt_h) * max_wt_h;
      YR = Y.each_row() - beta_g_h.slice(h).row(0);
    }
    if (diagonal) {
      YR = arma::square(YR);
      sigma_g_h.slice(h) = arma::diagmat(arma::sum(YR.each_col() % wt_h).t()) * max_wt_h;
    } else {
      sigma_g_h.slice(h) = weighted_cov(YR, wt_h) * max_wt_h;
    }
  }
  return;
}

// [[Rcpp::export]]
Rcpp::List continue_em(
    Rcpp::List state,
    const unsigned int niter,
    const bool hold_z = false,
    const bool progress = true,
    const bool loglik = true,
    const bool diagonal = false,
    const unsigned int gamma_updates = 100,
    const double gamma_epsilon = 1e-8,
    const double loglik_eps = 1e-4
){

  const arma::mat data = state["data"];
  const unsigned int px = state["px"], py = state["py"], Gx = state["Gx"], G = state["G"];
  const arma::uvec Gy = state["Gy"];
  const bool regression = state["regression"], subcluster = state["subcluster"];
  arma::mat z_ngh = state["z_ngh"], X = state["X"], Y = state["Y"], Z = state["Z"];
  Rcpp::List sigma_bar_g_h_inv = state["sigma.bar.g.h.inv"], mu_bar_g_h = state["mu.bar.g.h"];
  const arma::uvec setx = state["setx"], sety = state["sety"];
  arma::vec pi_g(Gx), pi_g_h(arma::accu(Gy)), tau_g_h = state["tau.g.h"];
  arma::mat mu_g = state["mu.g"];
  arma::cube sigma_g = state["sigma.g"], beta_g_h = state["beta.g.h"], sigma_g_h = state["sigma.g.h"];
  arma::mat Gamma = state["Gamma"];
  arma::mat zlog;
  arma::vec psi = state["psi"];
  const unsigned int n_obs = data.n_rows;

  Progress p(niter, progress);
  Rcpp::NumericVector ll(niter);
  ll.fill(NA_REAL);

  for (unsigned int i = 0; i < niter; i++) {
    if (p.increment()) {
      // e.step
      if (!hold_z) {z_ngh = e_step_helper(px, py, G, sigma_bar_g_h_inv, mu_bar_g_h, tau_g_h, X, Y);}

      // m.step
      m_step_primary_inplace(pi_g, mu_g, sigma_g, z_ngh, setx, Gx, X, diagonal);
      m_step_secondary_inplace(pi_g_h, beta_g_h, sigma_g_h, z_ngh, sety, Gy, X, Y, regression, subcluster, diagonal);
      psi = arma::trans(arma::mean(Z % Z, 0));

      // update.bars
      arma::cube sigma_bar_inv =
        update_bars_sigma_inv(G, setx, sety, sigma_g, sigma_g_h, beta_g_h, psi, px, py, Gamma.n_rows);

      arma::vec tmpone(1, arma::fill::ones);
      arma::vec tmpzero(Gamma.n_rows - px - py, arma::fill::zeros);
      for (unsigned int j = 0; j < G; j++) {
        arma::vec mu_bar = arma::join_vert(
          mu_g.col(setx(j) - 1),
          (arma::join_vert(tmpone, mu_g.col(setx(j) - 1)).t() *
            beta_g_h.slice(sety(j) - 1)).t(),
            tmpzero
        );
        mu_bar_g_h[j] = mu_bar;
        sigma_bar_g_h_inv[j] = sigma_bar_inv.slice(j);
      }

      pi_g /= n_obs; pi_g_h /= n_obs;
      for (unsigned int j = 0; j < G; j++) {tau_g_h(j) = pi_g(setx(j) - 1) * pi_g_h(sety(j) - 1);}

      // m.step.gamma
      Gamma = m_step_gamma_dual(G, sigma_bar_g_h_inv, mu_bar_g_h, z_ngh, data, Gamma, gamma_updates, gamma_epsilon);

      // update.xyz
      arma::mat xyz = data * Gamma;
      X = xyz.cols(0, px - 1);
      Y = xyz.cols(px, px + py - 1);
      if (px + py < Gamma.n_rows) {Z = xyz.cols(px + py, Gamma.n_rows - 1);}

      if (loglik) {
        zlog = logden_helper(G, sigma_bar_g_h_inv, mu_bar_g_h, tau_g_h, X, Y, Z, data);
        ll(i) = arma::accu(logden_loglik(zlog));
        if (i > 1) {
          if (ll(i) - ll(i - 1) < loglik_eps) {
            break;
          }
        }
      }
    }
  }

  state["z_ngh"] = z_ngh;
  state["sigma.bar.g.h.inv"] = sigma_bar_g_h_inv;
  state["mu.bar.g.h"] = mu_bar_g_h;
  state["tau.g.h"] = tau_g_h;
  state["mu.g"] = mu_g;
  state["sigma.g"] = sigma_g;
  state["beta.g.h"] = beta_g_h;
  state["sigma.g.h"] = sigma_g_h;
  state["Gamma"] = Gamma;
  state["psi"] = psi;
  state["ll"] = ll;
  state["X"] = X;
  state["Y"] = Y;
  state["Z"] = Z;

  return(state);
}
