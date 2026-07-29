// lnm_gibbs_core_sym.cpp — FULLY SYMMETRIC (baseline-free) logistic-normal
// multinomial Gibbs (now with country RE). All K categories treated identically:
//  * symmetric PG over all K (c_ij = logsumexp of the OTHER K-1, no baseline 0)
//  * B per-category with zero-sum (CLR/Msym) horseshoe prior applied IN the draw
//    (pulls B[v,j] toward mean(B[v,-j])) then centered -> exact zero-sum, order-invariant
//  * country RE per (group, cat), centered to zero-sum; u / Sigma_od K x K, centered
// PG via large-h normal approx (levels). Baseline-free everywhere (validated by a
// category-permutation invariance test).
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
inline double pgm(double b,double c){double a=std::fabs(c);if(a<1e-6)return b*0.25;return (b/(2.0*c))*std::tanh(c/2.0);}
inline double pgv(double b,double c){double a=std::fabs(c);if(a<1e-6)return b/24.0;double cc=std::min(a,30.0)*(c<0?-1:1);double ch=std::cosh(cc/2.0);return std::max(b*(std::sinh(cc)-cc)/(4.0*cc*cc*cc*ch*ch),1e-12);}

// [[Rcpp::export]]
List lnm_gibbs_core_sym(const arma::mat& X, const arma::mat& Y, const arma::vec& N,
                        const arma::mat& Xre, const arma::uvec& gstart, const arma::uvec& gsize,
                        const arma::uvec& hs_idx0, int niter, int nburn,
                        double A0, double slab_s2, double iw_df, const arma::mat& iw_scale,
                        double re_a, double re_b, bool use_re, bool use_hs, int thin, unsigned int seed,
                        const arma::uvec& re_pin, const arma::vec& re_support,
                        const arma::vec& fe_support, bool re_asis, bool use_tempering = false, double tempering_T0 = 0.05,
                        double init_jitter = 0.0) {
  arma::arma_rng::set_seed(seed);
  unsigned int n=X.n_rows, p=X.n_cols, K=Y.n_cols, G=gstart.n_elem, q=Xre.n_cols;
  arma::mat kappa(n,K); for(unsigned j=0;j<K;++j) kappa.col(j)=Y.col(j)-N/2.0;
  arma::mat B(p,K,arma::fill::zeros), u(n,K,arma::fill::zeros), bre(std::max(G*q,1u),K,arma::fill::zeros);
  if(init_jitter>0) B += init_jitter*(2.0*arma::randu<arma::mat>(p,K)-1.0);   // per-chain overdispersed init (uniform +/- jitter; arma RNG seeded per chain)
  arma::vec tau_re2(K,arma::fill::ones);
  arma::mat Sigma=arma::eye(K,K), SigInv=arma::eye(K,K);
  arma::vec lambda2(p,arma::fill::ones), nuhs(p,arma::fill::ones), prior_prec(p);
  prior_prec.fill(1.0/(A0*A0)); double tau2=1,xi=1,c2=slab_s2;
  arma::mat postB_sum(p,K,arma::fill::zeros), postSig_sum(K,K,arma::fill::zeros), bre_sum(std::max(G*q,1u),K,arma::fill::zeros);
  arma::vec maxB(niter), od(niter);
  int n_keep=(niter-nburn)/thin; if(n_keep<1)n_keep=1;
  arma::mat Bdraws(p*K,n_keep,arma::fill::zeros); int kept=0;
  auto re_fit=[&](arma::mat& RF){RF.zeros(); if(use_re) for(unsigned g=0;g<G;++g){unsigned r0=gstart[g],r1=gstart[g]+gsize[g]-1; RF.rows(r0,r1)=Xre.rows(r0,r1)*bre.rows(g*q,g*q+q-1);}};

  for(int it=0; it<niter; ++it){
    double temp_iter = 1.0;
    if (use_tempering && it <= nburn / 4 && nburn > 0) {
      temp_iter = tempering_T0 + (1.0 - tempering_T0) * ((double)it / (nburn / 4));
    }
    arma::vec N_iter = N * temp_iter;
    arma::mat kappa_iter = kappa * temp_iter;
    arma::mat RF(n,K); re_fit(RF);
    arma::mat eta = X*B + RF + u;
    arma::mat expe = arma::exp(arma::clamp(eta,-30.0,30.0));
    arma::vec tot = arma::sum(expe,1);
    arma::mat om(n,K), z(n,K);
    for(unsigned j=0;j<K;++j){
      arma::vec cj = arma::log(arma::clamp(tot-expe.col(j),1e-12,arma::datum::inf));
      for(unsigned i=0;i<n;++i){double o=R::rnorm(pgm(N_iter[i],eta(i,j)-cj[i]),std::sqrt(pgv(N_iter[i],eta(i,j)-cj[i])));om(i,j)=std::max(o,1e-9);}
      z.col(j)=kappa_iter.col(j)/om.col(j)+cj;
    }
    if(use_hs) for(arma::uword idx:hs_idx0){double ve=c2*tau2*lambda2[idx]/(c2+tau2*lambda2[idx]);double fs=(fe_support.n_elem>idx?fe_support[idx]:1.0);prior_prec[idx]=fs/std::max(ve,1e-12);}  // fe_support: global participation-ratio FE shrinkage
    // B | .  per category with zero-sum (Msym) prior
    arma::mat rB = z - RF - u; double km=(double)K/(double)(K-1);
    for(unsigned j=0;j<K;++j){
      arma::vec wj=om.col(j); arma::mat P=X.t()*(X.each_col()%wj);
      arma::vec pm=(arma::sum(B,1)-B.col(j))/(double)(K-1); arma::vec pp=prior_prec*km; P.diag()+=pp;
      arma::vec b=X.t()*(wj%rB.col(j))+pp%pm; arma::mat L=arma::chol(P,"lower");
      arma::vec mu=arma::solve(arma::trimatu(L.t()),arma::solve(arma::trimatl(L),b));
      B.col(j)=mu+arma::solve(arma::trimatu(L.t()),arma::randn<arma::vec>(p));
    }
    B.each_col()-=arma::mean(B,1);
    // country RE | .  per (group,cat), centered to zero-sum
    if(use_re){
      arma::mat rr=z-X*B-u;
      for(unsigned g=0;g<G;++g){unsigned r0=gstart[g],r1=gstart[g]+gsize[g]-1; arma::mat Xg=Xre.rows(r0,r1);
        arma::vec rsg(q); for(unsigned r=0;r<q;++r) rsg[r]=(re_support.n_elem>g*q+r?re_support[g*q+r]:1.0);  // within-group participation-ratio support factor (>=1, shrinks info-sparse RE slopes)
        for(unsigned j=0;j<K;++j){arma::vec wj=om.col(j).subvec(r0,r1); arma::mat P=Xg.t()*(Xg.each_col()%wj); P.diag()+=rsg/tau_re2[j]+1e-6;
          arma::vec b=Xg.t()*(wj%rr.col(j).subvec(r0,r1)); arma::mat L=arma::chol(P,"lower");
          arma::vec mu=arma::solve(arma::trimatu(L.t()),arma::solve(arma::trimatl(L),b));
          bre.submat(g*q,j,g*q+q-1,j)=mu+arma::solve(arma::trimatu(L.t()),arma::randn<arma::vec>(q));}
        arma::mat blk=bre.rows(g*q,g*q+q-1); blk.each_col()-=arma::mean(blk,1); bre.rows(g*q,g*q+q-1)=blk;
        for(unsigned r=0;r<q;++r) if(re_pin.n_elem>g*q+r && re_pin[g*q+r]) bre.row(g*q+r).zeros();}  // pin screened (unidentified) RE slopes
      for(unsigned j=0;j<K;++j){double ss=arma::dot(bre.col(j),bre.col(j)); tau_re2[j]=1.0/R::rgamma(re_a+G*q/2.0,1.0/(re_b+ss/2.0));}
      re_fit(RF);
      // RE-ASIS interweave: redraw each per-category RE scale tau_re[j] in the non-centered
      // parameterization (bre.col(j) = s * btil, btil fixed) via an MH step. The PG-weighted
      // likelihood proposal N(C/A,1/sqrt(A)) cancels, leaving the IG-prior ratio + s^2->s Jacobian.
      // Decouples tau from the REs -> better joint mixing. Default off; on==off in mean at re_idx=1.
      if(re_asis){
        for(unsigned j=0;j<K;++j){
          double s_old=std::sqrt(tau_re2[j]); if(!std::isfinite(s_old)||s_old<1e-8) continue;
          arma::vec dcol=RF.col(j)/s_old;                                   // NCP design (RF.col(j) = s_old * dcol)
          double A=arma::dot(om.col(j)%dcol,dcol), C=arma::dot(om.col(j)%dcol,rr.col(j));
          if(!std::isfinite(A)||A<=0) continue;
          double s_prop=R::rnorm(C/A,1.0/std::sqrt(A));
          // Guard: reject (stay) on non-finite or extreme per-step scale moves. Rejecting is a
          // valid MH move (keeps the target); without it an accepted huge f blows eta into the
          // Polya-Gamma overflow zone -> hard crash at full random slopes (LNM RE fragility).
          if(!std::isfinite(s_prop) || s_prop<=0.1*s_old || s_prop>=10.0*s_old) continue;
          double la=(-(re_a+1.0)*std::log(s_prop*s_prop)-re_b/(s_prop*s_prop)+std::log(s_prop))
                   -(-(re_a+1.0)*std::log(tau_re2[j])-re_b/tau_re2[j]+std::log(s_old));
          if(std::isfinite(la) && std::log(R::runif(0,1))<la){ bre.col(j)*=(s_prop/s_old); tau_re2[j]=s_prop*s_prop; }
        }
        re_fit(RF);
      }
    }
    // u | .  per pixel, centered to zero-sum
    arma::mat rb=z-X*B-RF;
    for(unsigned i=0;i<n;++i){arma::mat Pu=SigInv;Pu.diag()+=om.row(i).t();arma::mat L=arma::chol(Pu,"lower");
      arma::vec bi=(om.row(i)%rb.row(i)).t();arma::vec mu=arma::solve(arma::trimatu(L.t()),arma::solve(arma::trimatl(L),bi));
      u.row(i)=(mu+arma::solve(arma::trimatu(L.t()),arma::randn<arma::vec>(K))).t();}
    u.each_col()-=arma::mean(u,1);
    Sigma=arma::iwishrnd(iw_scale+u.t()*u+1e-6*arma::eye(K,K), iw_df+n); SigInv=arma::inv_sympd(Sigma);
    if(use_hs){
      auto ss=[&](arma::uword idx){arma::rowvec br=B.row(idx);double mf=arma::accu(br)/(double)K;return arma::accu(arma::square(br-mf));};
      for(arma::uword idx:hs_idx0){double s=ss(idx);nuhs[idx]=1.0/R::rgamma(1.0,1.0/(1.0+1.0/lambda2[idx]));
        lambda2[idx]=std::min(std::max(1.0/R::rgamma((K+1)/2.0,1.0/(1.0/nuhs[idx]+s/(2.0*tau2))),1e-10),1e6);}
      double ssa=0;for(arma::uword idx:hs_idx0)ssa+=ss(idx)/lambda2[idx];
      xi=1.0/R::rgamma(1.0,1.0/(1.0+1.0/tau2));
      tau2=std::min(std::max(1.0/R::rgamma((hs_idx0.n_elem*K+1)/2.0,1.0/(1.0/xi+ssa/2.0)),1e-10),1e6);
    }
    maxB[it]=arma::abs(B).max(); od[it]=arma::mean(Sigma.diag());
    if(it>=nburn){postB_sum+=B;postSig_sum+=Sigma;bre_sum+=bre;
      if(((it-nburn)%thin==0)&&kept<n_keep){Bdraws.col(kept)=arma::vectorise(B);++kept;}}
  }
  return List::create(_["B_sum"]=postB_sum,_["Sig_sum"]=postSig_sum,_["bre_mean"]=bre_sum/(niter-nburn),
                      _["nret"]=niter-nburn,_["maxB"]=maxB,_["od"]=od,_["Bdraws"]=Bdraws.cols(0,kept-1));
}
