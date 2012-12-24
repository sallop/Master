#include "FBM.hh"

size_t g_energy_idx  = 0;
size_t g_denergy_idx = 0;
size_t g_prop_idx    = -1;	// not use

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

// matrix is symmetrix, and w[0,:], w[:,0] is a mean values of x[i]
// but now not implemented yet.

int main(int argc, char **argv)
{
  int D=3, P=3;
  size_t K=5, M=1000, N=1000, T=10000;
  double alpha = 0.01;
  char prefix[32]="D", odir[32] = "dir-dat";

  int result;
  while ((result = getopt(argc,argv, "D:P:K:M:N:T:a:e:u:d:p:")) != -1){
    switch (result){
    case 'D': D = atoi(optarg); break;
    case 'P': P = atoi(optarg); break;
    case 'K': K = atoi(optarg); break;
    case 'M': M = atoi(optarg); break;
    case 'N': N = atoi(optarg); break;
    case 'T': T = atoi(optarg); break;
    case 'a': alpha = atof(optarg); break;
    case 'e': g_energy_idx  = atoi(optarg); break;
    case 'u': g_denergy_idx = atoi(optarg); break;
    case 'd': strcpy(odir  ,optarg); break;
    case 'p': strcpy(prefix,optarg); break;
    case '?': std::cout << "optarg=" << optarg << std::endl;
    }
  }

  _p(std::cout, D);_p(std::cout, P);_p(std::cout, K);
  _p(std::cout, M);_p(std::cout, N);_p(std::cout, T);
  _p(std::cout, g_energy_idx);
  _p(std::cout, g_denergy_idx);
  _p(std::cout, prefix);

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
    char fname[32]; sprintf(fname,"%s/%s%s",odir,prefix,idx->identify);
    idx->os->open(fname);
  }

  Layer L(D+1,D+1), data (D+1,D+1), model(D+1,D+1); // symmetrize
  symmetric_view sL(L), sdata(data), smodel(model);
  banded_view bL(L,0), bdata(data,0), bmomdel(model,0);
  Pattern v(D+1,0), vf(D+1,0);
  Probability boltz_P(0x01<<D, 0), gibbs_P(0x01<<D, 0);
  Probability prob_Q (0x01<<D, 0), accm_Q (0x01<<D, 0);

  double kl;

  prob_Q[0]=0.10;prob_Q[1]=0.10;prob_Q[2]=0.05;prob_Q[3]=0.05;
  prob_Q[4]=0.10;prob_Q[5]=0.10;prob_Q[6]=0.40;prob_Q[7]=0.10;

  std::partial_sum(prob_Q.begin(),prob_Q.end(),accm_Q.begin());
  accm_Q[accm_Q.size()-1] = 1.0;

  // modifies later. now degree is D, and nothing a special on bit
  Patterns all_v(0x01<<D, Pattern(D+1));
  for(size_t n=0; n < 0x01<<D; ++n){
    all_v[n] = dtob(n,D+1);
    all_v[n][0]=1;
  }
  //  assert(false);
  //1. Randomly initialize
  for(size_t i=0;i<D+1;i++)
    for(size_t j=0;j<D+1;j++)
      L(i,j)=drand48();
  // symmetrization
  L = sL - bL;
  {// print out simulation config
    test(D,P,K,M,N,T,
	 prefix,
	 L, data, model,sL,sdata,smodel,bL,bdata,bmomdel,
	 v,vf,boltz_P,gibbs_P,prob_Q,accm_Q,
	 alpha, kl, all_v, ofs_cfg);
  }

  for(size_t t=1; t<T; t++){
    std::cout << _eq(t) << std::endl;
    data  = zerom(D+1,D+1);
    model = zerom(D+1,D+1);
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
      vf = zerov(D+1);
      vf[0]=1;
      //kstep_gibbs(K,vf,L,alpha);
      kstep_gibbs(K,vf,L,1.0);
      model += outer_prod(vf,vf);
      gibbs_P[btod(vf)]++;

    }
    model /= M;
    gibbs_P /= M;

    // date's impliment is symmetry
    L = L + alpha*(data - model);
    L = sL - bL;
    //    alpha = 1.0/(t+1);
    
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
