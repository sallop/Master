#include "FBM.hh"


// matrix is symmetrix, and w[i,i] is a mean values of x[i]
int main(int argc, char **argv)
{
  size_t D=3, P=3;
  size_t K=5, M=1000, N=1000, T=10000;
  char prefix[32]="K3", odir[32] = "dir-dat";
  g_energy_idx  =  8;
  g_denergy_idx =  4;
  g_prop_idx    = -1;
  g_annealing_schedule_idx = 0;
// [kurata's matrix _energy3, _energy8] : [date's matrix _energy0]
// [kurata's matrix _denergy1, _energy4] : [date's matrix _denergy0]

  std::fstream ofs_L, ofs_t_gibbs, ofs_cfg;
  cook_args(argc, argv, D, P, K, M, N, T, prefix, odir);
  init_plot_table(odir, prefix, ofs_L, ofs_t_gibbs, ofs_cfg, std::ios_base::out);
  _p(std::cout, D);_p(std::cout, P);_p(std::cout, K);
  _p(std::cout, M);_p(std::cout, N);_p(std::cout, T);
  _p(std::cout, prefix);

  


  Layer L(D,D), data (D,D), model(D,D); // symmetrize
  symmetric_view sL(L);
  symmetric_view sdata(data);
  symmetric_view smodel(model);
  banded_view bL(L,0);		// exist for test(), but not use.
  banded_view bdata(data,0);	// exist for test(), but not use.
  banded_view bmomdel(model,0);	// exist for test(), but not use.

  Pattern v(D,0), vf(D,0);
  Probability boltz_P(0x01<<D, 0);
  Probability gibbs_P(0x01<<D, 0);
  Probability prob_Q (0x01<<D, 0);
  Probability accm_Q (0x01<<D, 0);
  double alpha = 1.0;
  double kl;
  
  prob_Q[0]=0.10;prob_Q[1]=0.10;prob_Q[2]=0.05;prob_Q[3]=0.05;
  prob_Q[4]=0.10;prob_Q[5]=0.10;prob_Q[6]=0.40;prob_Q[7]=0.10;

  std::partial_sum(prob_Q.begin(),prob_Q.end(),accm_Q.begin());
  accm_Q[accm_Q.size()-1] = 1.0;

  // modifies later. now degree is D, and nothing a special on bit
  Patterns all_v(0x01<<D, Pattern(D));
  for(size_t n=0; n < 0x01<<D; ++n)
    all_v[n] = dtob(n,D);

  //1. Randomly initialize
  for(size_t i=0;i<D;i++)
    for(size_t j=0;j<D;j++)
      L(i,j)=drand48();
  // symmetrization
  L = sL;
  {// cfg-closed
    symcfg(D,P,K,M,N,T,
	   prefix,
	   L, data, model,sL,sdata,smodel,bL,bdata,bmomdel,
	   v,vf,boltz_P,gibbs_P,prob_Q,accm_Q,
	   alpha, kl, all_v, ofs_cfg);
  }

  for(size_t t=1; t<T; t++){
    std::cout << _eq(t) << std::endl;
    data  = zerom(D,D);
    model = zerom(D,D);
    // data dependent expection
    for(size_t n=0; n<N; n++){
      size_t k;
      double r = drand48();
      for (k=0; k < accm_Q.size(); k++)
	if (r < accm_Q[k])
	  break;
      data += ub::outer_prod(all_v[k],all_v[k]);
    }
    data /= N;

    // model dependent expection
    gibbs_P = zerov(0x01<<D);
    vf = zerov(D);
    for(size_t m=0; m<M; m++){
      kstep_gibbs(1,vf,L,alpha);
      model += outer_prod(vf,vf);
    }
    model /= M;
    gibbs_P /= M;

    // kurata's impliment is symmetry
    L = L + alpha*(data - model);
    // but diagonal elements store mean values
    // L = sL - bL;
    // annealing schedule
    alpha = annealing_schedule(t);
    // output data
    ofs_L       << t << "\t" << L       << std::endl;
    ofs_t_gibbs << t << "\t" << gibbs_P << std::endl;
  }

  finalize_plot_table();

  for (size_t i=0; i<all_v.size(); i++)
    std::cout << "v^{"<<_eq(i)<<"}=" << all_v[i] << std::endl;

  return 0;
}
