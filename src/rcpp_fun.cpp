#include <Rcpp.h>
#include <RcppEigen.h>
#include <queue>
#include <iostream>
#include <math.h>
using namespace Rcpp;
using namespace RcppEigen;

// [[Rcpp::depends(RcppEigen)]]

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::ArrayXd;
using Eigen::VectorXi;
using namespace std;

int isnan(double x) { return x != x; }
int isinf(double x) { return !isnan(x) && isnan(x - x); }

MatrixXd AtA(const MatrixXd & A) { // transpose of A times itself
  int n(A.cols());
  return MatrixXd(n,n).setZero().selfadjointView<Eigen::Lower>().rankUpdate(A.adjoint());
}

MatrixXd AAt (const MatrixXd & A) { // A times its transpose
  int m(A.rows());
  return MatrixXd(m, m).setZero().selfadjointView<Eigen::Lower>().rankUpdate(A);
}

struct line_search_par {
  bool line_search;
  double ls_beta;
  double ls_eps;
  int ls_maxit;
};

struct TISP_fun_grad_inp {
  MatrixXd X_affinity;
  VectorXd y_affinity;
  double y_sc;
};

struct optim_par {
  int maxit;
  double eps;
  double step_size;
  bool accelerate;
};

struct thresh_par {
  int q;
  double lambda;
  double eta;
  double lambda_sc;
};

//// thresholding functions
vector<int> top_i_pq(const VectorXd & v, const int & n) {
  typedef pair<double, int> Elt;
  priority_queue< Elt, vector<Elt>, greater<Elt> > pq;
  vector<int> result;
  
  for (int i = 0; i != v.size(); ++i) {
    if (pq.size() < n)
      pq.push(Elt(v(i), i));
    else {
      Elt elt = Elt(v(i), i);
      if (pq.top() < elt) {
        pq.pop();
        pq.push(elt);
      }
    }
  }
  
  result.reserve(pq.size());
  while (!pq.empty()) {
    result.push_back(pq.top().second);
    pq.pop();
  }
  
  return result;
}


VectorXd quantile_thresh(const VectorXd & inp, const thresh_par & thresh_inp) {
  int n = inp.size();
  VectorXd zs = VectorXd::Zero(n);
  if (thresh_inp.q > 0) {
    VectorXd inp_abs = inp.cwiseAbs();
    vector<int> order_idx = top_i_pq(inp_abs, thresh_inp.q);
    int curr_idx;
    for(int i = 0; i < thresh_inp.q; i++) {
      curr_idx = order_idx[i];
      zs(curr_idx) = inp(curr_idx);
    }
  } else {
    zs = inp;
  }
  return zs;
}

VectorXd quantile_ridge_thresh(const VectorXd & inp, 
                               const thresh_par & thresh_inp) {
  VectorXd thresh = quantile_thresh(inp, thresh_inp);
  return thresh / (1.0 + thresh_inp.eta);
}

VectorXd hard_thresh(const VectorXd & inp,
                     const thresh_par & thresh_inp) {
  int n = inp.size();
  VectorXd outp = VectorXd::Zero(n);
  for(int i = 0; i < n; i++) {
    if (abs(inp[i]) > thresh_inp.lambda_sc) {
      outp[i] = inp[i] - thresh_inp.lambda_sc;
    }
  }
  return outp;
}

VectorXd hard_ridge_thresh(const VectorXd & inp,
                           const thresh_par & thresh_inp) {
  return hard_thresh(inp / (1.0 + thresh_inp.eta), thresh_inp);
}

double sign_func(double x)
{
  if (x > 0)
    return +1.0;
  else if (x == 0)
    return 0.0;
  else
    return -1.0;
}

VectorXd soft_thresh(const VectorXd & inp,
                     const thresh_par & thresh_inp) {
  int n = inp.size();
  VectorXd outp = VectorXd::Zero(n);
  for(int i = 0; i < n; i++) {
    if (inp[i] > thresh_inp.lambda_sc) {
      outp[i] = inp[i] - thresh_inp.lambda_sc;
    } else if (inp[i] < -thresh_inp.lambda_sc) {
      outp[i] = inp[i] + thresh_inp.lambda_sc;
    }
  }
  return outp;
}
////

double fun_TISP_cpp(const VectorXd & beta,
                    const TISP_fun_grad_inp & extra_par) {
  return (extra_par.X_affinity * beta).dot(beta) / 2.0 - extra_par.y_affinity.dot(beta);
}

VectorXd grad_TISP_cpp(const VectorXd & beta,
                     const TISP_fun_grad_inp & extra_par) {
  return extra_par.X_affinity * beta - extra_par.y_affinity;
}

double surrogate_prox_cpp(const VectorXd & beta_new,
                         const VectorXd & beta_curr,
                         const VectorXd & grad_curr,
                         const double & fval_curr,
                         const double & step_size) {
  return fval_curr + grad_curr.dot(beta_new - beta_curr) + 1.0 / (2.0 * step_size) * (beta_new - beta_curr).squaredNorm();
}


void TISP_line_search(VectorXd & beta_new, 
                      const VectorXd & beta_curr, 
                      const VectorXd & grad_curr, 
                      const TISP_fun_grad_inp & extra_par,
                      double & step_size,
                      thresh_par & thresh_inp,
                      VectorXd (*thresh_fun)(const VectorXd &, const thresh_par &),
                      const line_search_par & ls_inp) {
  double fval_curr = fun_TISP_cpp(beta_curr, extra_par);
  double fun_diff = fun_TISP_cpp(beta_new, extra_par) - surrogate_prox_cpp(beta_new, beta_curr, grad_curr, fval_curr, step_size);
  int ls_iter = 0;
  while (fun_diff > 0 && abs(fun_diff) > ls_inp.ls_eps && ls_iter < ls_inp.ls_maxit) {
    ls_iter += 1;
    step_size *= ls_inp.ls_beta;
    thresh_inp.lambda_sc = thresh_inp.lambda * step_size;
    beta_new = (*thresh_fun)(beta_curr - step_size * grad_curr, thresh_inp);
    fun_diff = fun_TISP_cpp(beta_new, extra_par) - surrogate_prox_cpp(beta_new, beta_curr, grad_curr, fval_curr, step_size);
  }
}

void TISP_fwd_line_search(VectorXd & beta_new, 
                          const VectorXd & beta_curr, 
                          const VectorXd & grad_curr, 
                          const TISP_fun_grad_inp & extra_par,
                          double & step_size,
                          thresh_par & thresh_inp,
                          VectorXd (*thresh_fun)(const VectorXd &, const thresh_par &),
                          const line_search_par & ls_inp) {
  double fval_curr = fun_TISP_cpp(beta_curr, extra_par);
  double fun_diff = fun_TISP_cpp(beta_new, extra_par) - surrogate_prox_cpp(beta_new, beta_curr, grad_curr, fval_curr, step_size);
  int ls_iter = 0;
  VectorXd beta_ls_old = beta_new;
  while (ls_iter < ls_inp.ls_maxit) {
    ls_iter += 1;
    step_size *= ls_inp.ls_beta;
    thresh_inp.lambda_sc = thresh_inp.lambda * step_size;
    beta_new = (*thresh_fun)(beta_curr - step_size * grad_curr, thresh_inp);
    fun_diff = fun_TISP_cpp(beta_new, extra_par) - surrogate_prox_cpp(beta_new, beta_curr, grad_curr, fval_curr, step_size);
    if (fun_diff > 0 || abs(fun_diff) > ls_inp.ls_eps) {
      beta_new = beta_ls_old;
      break;
    } else {
      beta_ls_old = beta_new;
    }
  }
}

List TISP_loop(const VectorXd & beta_init,
               const TISP_fun_grad_inp & extra_par,
               thresh_par & thresh_inp,
               const line_search_par & ls_inp,
               const optim_par & optim_inp,
               VectorXd (*thresh_fun)(const VectorXd &, const thresh_par &),
               const std::string & ls_type,
               const std::string & tol_type) {
  VectorXd beta_curr = beta_init;
  VectorXd grad_val = beta_curr;
  VectorXd beta_new = beta_curr;
  VectorXd update_term = beta_curr;
  VectorXd accel_update = beta_curr;
  double omega = 2.0;
  double one_min_omega = (1.0 - omega);
  double loss_curr = fun_TISP_cpp(beta_curr, extra_par);
  double loss_new;
  double tol;
  double step_size_sc = optim_inp.step_size;
  VectorXd loss_vals = VectorXd::Zero(optim_inp.maxit);
  int j = 0;
  bool converged = false;
  while (!converged && j < optim_inp.maxit) {
    if(j % 100 == 0)
      Rcpp::checkUserInterrupt();
    j += 1;
    grad_val = grad_TISP_cpp(beta_curr, extra_par);
    update_term = beta_curr - optim_inp.step_size * grad_val;
    if (optim_inp.accelerate) {
      accel_update = one_min_omega * accel_update + omega * update_term;
      beta_new = (*thresh_fun)(accel_update, thresh_inp);
    } else {
      beta_new = (*thresh_fun)(update_term, thresh_inp);
      if (ls_inp.line_search) {
        if (!ls_type.compare("forward")) {
          step_size_sc = optim_inp.step_size;
          thresh_inp.lambda_sc = thresh_inp.lambda * step_size_sc;
          TISP_fwd_line_search(beta_new, beta_curr, grad_val, extra_par, step_size_sc, thresh_inp, thresh_fun, ls_inp);
        } else if (!ls_type.compare("backward")) {
          step_size_sc = 1.0;
          thresh_inp.lambda_sc = thresh_inp.lambda * step_size_sc;
          TISP_line_search(beta_new, beta_curr, grad_val, extra_par, step_size_sc, thresh_inp, thresh_fun, ls_inp);
        }
      }
    }
    loss_new = fun_TISP_cpp(beta_new, extra_par);
    if (!tol_type.compare("loss")) {
      tol = abs(loss_new - loss_curr);
    } else if (!tol_type.compare("beta")) {
      tol = (beta_new - beta_curr).array().abs().maxCoeff();
    }
    loss_curr = loss_new;
    loss_vals(j) = loss_curr;
    beta_curr = beta_new;
    if (tol < optim_inp.eps) {
      converged = true;
    }
  }
  loss_vals = loss_vals.head(j);
  double loss_final = loss_curr + extra_par.y_sc / 2.0;
  List ret;
  ret["beta_opt"] = beta_curr;
  ret["iter"] = j;
  ret["step_size"] = step_size_sc;
  ret["eta"] = thresh_inp.eta;
  ret["tol"] = tol;
  ret["loss_final"] = loss_final;
  ret["loss_vals"] = loss_vals;
  ret["converged"] = converged;
  return ret;
}

// [[Rcpp::export]]
List TISP_inner(const Eigen::Map<Eigen::MatrixXd> & X_affinity,
                const Eigen::Map<Eigen::VectorXd> & y_affinity,
                const Eigen::Map<Eigen::VectorXd> & beta_init,
                const double & y_sc,
                const double & lambda,
                const int & q,
                double & step_size,
                const bool & accelerate = true,
                const bool & line_search = true,
                const double & ls_beta = 1.5,
                const double & ls_eps = 1e-10,
                const int & ls_maxit = 100,
                const double & eta = 1e-4,
                const std::string & thresh_type = "quantile_ridge",
                const std::string & tol_type = "beta",
                const int & maxit = 2000,
                const double & eps = 1e-5,
                const std::string & ls_type = "forward") {
  line_search_par ls_inp = {line_search, ls_beta, ls_eps, ls_maxit};
  TISP_fun_grad_inp extra_par = {X_affinity, y_affinity, y_sc};
  thresh_par thresh_inp = {q, lambda, eta, lambda * step_size};
  optim_par optim_inp = {maxit, eps, step_size, accelerate};
  List ret;
  if (!thresh_type.compare("quantile_ridge")) {
    ret = TISP_loop(beta_init, extra_par, thresh_inp, ls_inp, optim_inp, &quantile_ridge_thresh, ls_type, tol_type);
  } else if (!thresh_type.compare("quantile")) {
    ret = TISP_loop(beta_init, extra_par, thresh_inp, ls_inp, optim_inp, &quantile_thresh, ls_type, tol_type);
  } else if (!thresh_type.compare("soft")) {
    ret = TISP_loop(beta_init, extra_par, thresh_inp, ls_inp, optim_inp, &soft_thresh, ls_type, tol_type);
  } else if (!thresh_type.compare("hard_ridge")) {
    ret = TISP_loop(beta_init, extra_par, thresh_inp, ls_inp, optim_inp, &hard_ridge_thresh, ls_type, tol_type);
  } else if (!thresh_type.compare("hard")) {
    ret = TISP_loop(beta_init, extra_par, thresh_inp, ls_inp, optim_inp, &hard_thresh, ls_type, tol_type);
  } else {
    stop("thresholding type not supported");
  }
  return ret;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//
