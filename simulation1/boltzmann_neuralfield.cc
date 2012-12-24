#include "BM.hh"
#include "util.hh"
// matrix is symmetrix, and w[i,i] is a mean values of x[i]
int main(int argc, char **argv)
{
  int D=4, P=4;
  int K=5, M=1000, N=1000, T=100;
  char prefix[32]="bnf", odir[32] = "dir-dat";
  double connection_scale = 1.0;
  g_energy_idx  = 0;
  g_denergy_idx = 0;
  g_prop_idx    = 0;
  g_annealing_schedule_idx = 0;
  g_init_layer_idx = 0;

// [kurata's matrix _energy3, _energy8] : [date's matrix _energy0]
// [kurata's matrix _denergy1, _energy4] : [date's matrix _denergy0]
  std::fstream ofs_W, ofs_t_gibbs, ofs_cfg;

  cook_args(argc, argv, D, P, K, M, N, T, 
	    connection_scale, prefix, odir);

  init_plot_table(odir, prefix, ofs_W, ofs_t_gibbs, ofs_cfg,
		  std::ios_base::out);

  _p(std::cout, D);
  _p(std::cout, P);
  _p(std::cout, K);
  _p(std::cout, M);
  _p(std::cout, N);
  _p(std::cout, T);
  _p(std::cout, prefix);

  Layer W(D,P), data (D,P), model(D,P); // symmetrize
  symmetric_view sW(W)  , sdata(data)  , smodel (model)  ;
  banded_view    bW(W,0), bdata(data,0), bmomdel(model,0);

  Pattern v(D,0), vf(D,0);
  Probability boltz_P(0x01<<D, 0);
  Probability gibbs_P(0x01<<D, 0);
  Probability prob_Q (0x01<<D, 0);
  Probability accm_Q (0x01<<D, 0);
  double alpha = 1.0;
  double temperature = 1.0;

  double kl;
  
  prob_Q[0]=0.10; prob_Q[1]=0.10; prob_Q[2]=0.05; prob_Q[3]=0.05;
  prob_Q[4]=0.10; prob_Q[5]=0.10; prob_Q[6]=0.40; prob_Q[7]=0.10;

  std::partial_sum(prob_Q.begin(),prob_Q.end(),accm_Q.begin());
  accm_Q[accm_Q.size()-1] = 1.0;

  // minor change.
  // only adjacent 2 cell is file
  Patterns all_v(D, Pattern(D));
  //for(size_t n=0; n < 0x01<<D; ++n)
  //all_v[n] = dtob(n,D);
  

  for (size_t i=0; i<all_v.size(); i++)
    _p(std::cout, all_v[i]);
  assert(false);
  //1. Randomly initialize
  // need to minor change.
  // boltzmann neural field.
  W = init_layer(D,connection_scale);
  
  {// cfg-closed
    symcfg(D, P, K, M, N, T, 
	   prefix,
	   W,  data,  model,
	   sW, sdata, smodel,
	   bW, bdata, bmomdel,
	   v, vf,
	   boltz_P, gibbs_P,
	   prob_Q , accm_Q,
	   alpha, kl,
	   all_v,
	   connection_scale,
	   ofs_cfg);
  }

  for(size_t t=1; t<T; t++){
    _p(std::cout, t);
    data  = zerom(D,P);
    model = zerom(D,P);
    // data dependent expection
    // need to minor change
    // Input Layer procedure
    // hebb learning
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
    // need to minor change
    // Output Layer procedure. 
    // anti hebb learning
    gibbs_P = zerov(0x01<<D);
    for(size_t m=0; m<M; m++){
      vf = zerov(D);
      kstep_gibbs(K,vf,W,temperature);
      model += outer_prod(vf,vf);
      gibbs_P[btod(vf)]++;
    }
    model /= M;
    gibbs_P /= M;

    // learning update 
    W = W + alpha*(data - model);
    

    alpha = annealing_schedule(t);
    
    // output data
    ofs_W       << t << "\t" << W       << std::endl;
    ofs_t_gibbs << t << "\t" << gibbs_P << std::endl;
  }

  finalize_plot_table();

  std::cout << "check all_v about side effect" << std::endl;
  for (size_t i=0; i<all_v.size(); i++)
    std::cout << "v^{"<<_eq(i)<<"}=" << all_v[i] << std::endl;

  return 0;
}
