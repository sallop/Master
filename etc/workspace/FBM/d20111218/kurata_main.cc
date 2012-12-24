#include "FBM.hh"

size_t g_energy_idx  = 10;
size_t g_denergy_idx =  4;
size_t g_prop_idx    = -1;

double energy(Pattern &v, Layer &L)
{
  static energy_ptr f = g_energy_tbl[g_energy_idx];
  return f(v, L);
}

double denergy(size_t i, Pattern &v, Layer &L)
{
  static denergy_ptr f = g_denergy_tbl[g_denergy_idx];
  return f(i,v,L);
}

double prop(int i, Pattern& v, Layer& L, double temperature)
{
  return sigma(denergy(i,v,L)/temperature);
}

// matrix is symmetrix, and w[i,i] is a mean values of x[i]

int main(int argc, char **argv)
{
  int D=3, P=3;
  size_t K=5, M=1000, N=1000, T=10000;
  char prefix[32]="K", odir[32] = "dir-dat";
  bool f_test = false;

  int result;
  while ((result = getopt(argc,argv, "D:P:K:M:N:T:d:e:u:p:")) != -1){
    switch (result){
    case 'D': D = atoi(optarg); break;
    case 'P': P = atoi(optarg); break;
    case 'K': K = atoi(optarg); break;
    case 'M': M = atoi(optarg); break;
    case 'N': N = atoi(optarg); break;
    case 'T': T = atoi(optarg); break;
    case 'd': strcpy(odir  ,optarg); break;
    case 'p': strcpy(prefix,optarg); break;
    case 'e': g_energy_idx = atoi(optarg); break;
    case 'u': g_denergy_idx = atoi(optarg); break;
    case '?': std::cout << "optarg=" << optarg << std::endl;
    }
  }

  _p(std::cout, D);_p(std::cout, P);_p(std::cout, K);
  _p(std::cout, M);_p(std::cout, N);_p(std::cout, T);
  _p(std::cout, prefix);_p(std::cout, f_test);

  std::ofstream ofs_L, ofs_t_gibbs, ofs_cfg;

  struct {/* output stream & identifies */
    std::ofstream *os; const char *identify;
  } *idx,ftbl[] = {
    {&ofs_L      ,"-mat.dat"},
    {&ofs_t_gibbs,"-t_gibbs.dat"},
    {&ofs_cfg    ,"-cfg.dat"},
    {NULL, NULL},
  };
  for (idx=&ftbl[0]; idx->identify != NULL; idx++){
    char fname[32];
    sprintf(fname,"%s/%s%s",odir, prefix, idx->identify);
    idx->os->open(fname);
  }

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
    test(D,P,K,M,N,T,
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
    for(size_t m=0; m<M; m++){
      vf = zerov(D);
      kstep_gibbs(K,vf,L,alpha);
      model += outer_prod(vf,vf);
      gibbs_P[btod(vf)]++;
    }
    model /= M;
    gibbs_P /= M;

    // kurata's impliment is symmetry
    L = L + alpha*(data - model);
    // but diagonal elements store mean values
    // L = sL - bL;

    // annealing schedule
#ifdef ANNEALING == 0
    alpha = 0.01;
#elif ANNEALING == 1
    // c being greater than or equal to the largest energy barrier
    // in the problem.
    // d is usually set equal to one.
    #define __c 5.0
    #define __d 1.0
    alpha = 1.0/(log(t+__d));
    #undef __c
    #undef __d
#elif ANNEALING == 2
    #define __D D
    alpha = (__D - 1)/log(t);
    #undef
#else
    alpha = 1.0/(t+1);
#endif
    // output data
    ofs_L       << t << "\t" << L       << std::endl;
    ofs_t_gibbs << t << "\t" << gibbs_P << std::endl;
  }

  for (idx=&ftbl[0]; idx->identify != NULL; idx++)
    idx->os->close();

  for (size_t i=0; i<all_v.size(); i++)
    std::cout << "v^{"<<_eq(i)<<"}=" << all_v[i] << std::endl;

  return 0;
}
