#include "BM.hh"
#include "config.hh"

int main(int argc, char **argv)
{
  int D=3, P=1;
  int K=10, M=1000, N=1000, T=1000;
  Layer F(D+P,D+P);

  char prefix[32]="g2", odir[32] = "dir-dat";
  double connection_scale = 1.0;
  g_energy_idx  = 0;
  g_denergy_idx = 0;
  g_prop_idx    = 0;
  g_annealing_schedule_idx = 1;
  g_init_layer_idx = 0;

  std::fstream ofs_F, ofs_t_gibbs, ofs_cfg;
  cook_args(argc, argv, D, P, K, M, N, T,
	    connection_scale, prefix, odir);

  init_plot_table(odir, prefix, 
		  ofs_F, ofs_t_gibbs, ofs_cfg,
		  std::ios_base::out);

  banded_view Fb(F,0);
  SubLayer L(F ,ub::range(0,D  ),ub::range(0,D  ));
  SubLayer W(F ,ub::range(0,D  ),ub::range(D,D+P));
  SubLayer J(F ,ub::range(D,D+P),ub::range(D,D+P));
  Layer Ld(D,D), Lm(D,D);
  Layer Wd(D,P), Wm(D,P);
  Layer Jd(P,P), Jm(P,P);
  symmetric_view Fs(F);
  symmetricsub_view Ls(L), Js(J);

  Probability prob_Q (0x01<<D,0); // environment
  Probability accm_Q (prob_Q.size(),0);
  Probability boltz_P(prob_Q.size(),0);
  Probability gibbs_P(D,0); // can't print KLd. so alternate vf mean
  double alpha = 1.0;	    // learning rate
  double temperature = 1.0; // temperature
  double kl;	// learning rate, store the kl divergence
  

  // all state
  Patterns all_v (0x01<<D,Pattern(D,0));
  Patterns all_h (0x01<<P,Pattern(P,0));
  Patterns all_mu(all_v.size(),Pattern(P,0));

  //PatternTbl bton_tbl;
  //  BtonTbl npattern_to_digit;
  DtoB_table dtob_tbl;
  //  DtoB_value ptod_tbl[D];

  // fantasy particle
  Patterns all_vf(M,Pattern(D,0));// fantasy particle
  Patterns all_hf(M,Pattern(P,0));// fantasy particle
    
  if (D == 3){
    prob_Q[0]=0.10; prob_Q[1]=0.10; prob_Q[2]=0.05; prob_Q[3]=0.05;
    prob_Q[4]=0.10; prob_Q[5]=0.10; prob_Q[6]=0.40; prob_Q[7]=0.10;
  } 
  else if (D == 4){
    size_t p = 0;
    double q = 0x01<<(D-3);
    prob_Q[p+0]=0.10/q; prob_Q[p+1]=0.10/q; prob_Q[p+2]=0.05/q; prob_Q[p+3]=0.05/q;
    prob_Q[p+4]=0.10/q; prob_Q[p+5]=0.10/q; prob_Q[p+6]=0.40/q; prob_Q[p+7]=0.10/q;
    p = 0x01 << 3;
    prob_Q[p+0]=0.10/q; prob_Q[p+1]=0.10/q; prob_Q[p+2]=0.05/q; prob_Q[p+3]=0.05/q;
    prob_Q[p+4]=0.10/q; prob_Q[p+5]=0.10/q; prob_Q[p+6]=0.40/q; prob_Q[p+7]=0.10/q;
  }
  else if (D == 5){
    size_t p = 0;
    double q = 0x01<<(D-3);
    prob_Q[p+0]=0.10/q; prob_Q[p+1]=0.10/q; prob_Q[p+2]=0.05/q; prob_Q[p+3]=0.05/q;
    prob_Q[p+4]=0.10/q; prob_Q[p+5]=0.10/q; prob_Q[p+6]=0.40/q; prob_Q[p+7]=0.10/q;
    p = 0x01 << 3;
    prob_Q[p+0]=0.10/q; prob_Q[p+1]=0.10/q; prob_Q[p+2]=0.05/q; prob_Q[p+3]=0.05/q;
    prob_Q[p+4]=0.10/q; prob_Q[p+5]=0.10/q; prob_Q[p+6]=0.40/q; prob_Q[p+7]=0.10/q;
    p += 0x01 << 3;
    prob_Q[p+0]=0.10/q; prob_Q[p+1]=0.10/q; prob_Q[p+2]=0.05/q; prob_Q[p+3]=0.05/q;
    prob_Q[p+4]=0.10/q; prob_Q[p+5]=0.10/q; prob_Q[p+6]=0.40/q; prob_Q[p+7]=0.10/q;
    p += 0x01 << 3;
    prob_Q[p+0]=0.10/q; prob_Q[p+1]=0.10/q; prob_Q[p+2]=0.05/q; prob_Q[p+3]=0.05/q;
    prob_Q[p+4]=0.10/q; prob_Q[p+5]=0.10/q; prob_Q[p+6]=0.40/q; prob_Q[p+7]=0.10/q;
  }

  std::partial_sum(prob_Q.begin(),prob_Q.end(),accm_Q.begin());

  std::cout << _eq(prob_Q) << std::endl;
  std::cout << _eq(accm_Q) << std::endl;

  assert(accm_Q[accm_Q.size()-1] > 0.99 && accm_Q[accm_Q.size()-1] < 1.01);
  accm_Q[accm_Q.size()-1];
  for (size_t i=0; i<all_v.size(); i++){
    all_v[i] = dtob(i,D);
    dtob_tbl[btod(all_v[i])] = make_pair(i,all_v[i]);
  }
  for (size_t j=0; j<all_h.size(); j++)
    all_h[j] = dtob(j,P);

  //1. Randomly initialize mu and
  F = init_layer(D+P,connection_scale);
  for(size_t m=0; m<M; m++){
    for(size_t i=0; i<D; i++)
      all_vf[m][i] = drand48()>0.5?ONBIT:OFFBIT;
    for(size_t i=0; i<P; i++)
      all_hf[m][i] = drand48()>0.5?ONBIT:OFFBIT;
  }

  {// plot configuration
    symcfg( D, P, K, M, N, T,
	    prefix,
	    F,L,W,J,
	    //all_vf,all_hf,
	    prob_Q,accm_Q,
	    alpha,
	    all_v,all_h,all_mu,
	    connection_scale,
	    ofs_cfg);
  }

  symcfg( D, P, K, M, N, T,
	  prefix,F,L,W,J,
	  prob_Q,accm_Q,alpha,
	  all_v,all_h,all_mu,
	  connection_scale,std::cout);
  //finalize_plot_table(); assert(false);

  const Layer L0 = zerom(D,D);
  const Layer W0 = zerom(D,P);
  const Layer J0 = zerom(P,P);
  const Probability gibbs_P0 = zerov(gibbs_P.size());
  for(size_t t=0; t<T; t++){
    // zero initialize.
    Ld = Lm = L0;
    Wd = Wm = W0;
    //(a) n=1 to N, data dependent expection
    for(size_t n=0; n<N; n++){
      size_t k;
      double r = drand48();
      for (k=0; k<accm_Q.size(); k++)
	if (r < accm_Q[k])
	  break;
      Pattern &vn = all_v[k], &mu = all_mu[k];
      // mean-field updates until convergence
      // this values convergence about 10-20
      for(size_t i=0; i < 32; i++){
	size_t j = rand()%P;
	mu[j] = sigma(denergy_up(j,vn,mu,L,J,W));
      }
      Ld += ub::outer_prod(vn,vn);
      Wd += ub::outer_prod(vn,mu);
      Jd += ub::outer_prod(mu,mu);
    }
    Ld /= N;
    Wd /= N;
    Jd /= N;

    //(b) m=1 to M
    gibbs_P = gibbs_P0;
    for(size_t m=0; m<M; m++){
      Pattern &vf = all_vf[m], &hf = all_hf[m];
      // k-step Gibbs sampler using 4,5
      for(size_t k=0; k<K; k++){
	size_t j = rand()%P;
	size_t i = rand()%D;
	hf[j] = (sigma(denergy_up(j,vf,hf,L,J,W))>drand48()) ? ONBIT : OFFBIT;
	vf[i] = (sigma(denergy_dw(i,vf,hf,L,J,W))>drand48()) ? ONBIT : OFFBIT;
      }
      gibbs_P += vf;
      Lm += ub::outer_prod(vf,vf);
      Wm += ub::outer_prod(vf,hf);
      Jm += ub::outer_prod(hf,hf);
    }
    Lm /= M;
    Wm /= M;
    Jm /= M;
    gibbs_P /= M;

    {//(c) update
      L = L + alpha*(Ld - Lm);
      W = W + alpha*(Wd - Wm);
      J = J + alpha*(Jd - Jm);
    }
    for (size_t i=0; i<all_v.size(); i++)
      boltz_P[i] = prob_v(all_v[i],L,J,W,all_v,all_h);
    
    //alpha = 0.01;
    alpha = annealing_schedule(t);
    ofs_F << t << "\t" << F << std::endl;
    ofs_t_gibbs << t << "\t" << gibbs_P << std::endl;
  }

  finalize_plot_table();

  return 0;
}
