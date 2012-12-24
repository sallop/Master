#include "util.hh"
#include "BM.hh"
//#include "FBM.hh" FBM.cc include energy_tbl[] impliments, so you don't include this
//#include "GBM.hh" if you include GBM, or FBM, then may be duplicate error.

PlotTable g_plot_table;

void
init_plot_table(char *dir, char *prefix,
		std::fstream& fs_mat,
		std::fstream& fs_t_gibbs,
		std::fstream& fs_cfg,
		std::ios_base::openmode mode)
{
  g_plot_table.insert(PlotRecord("-mat.dat"    ,&fs_mat));
  g_plot_table.insert(PlotRecord("-t_gibbs.dat",&fs_t_gibbs));
  g_plot_table.insert(PlotRecord("-cfg.dat",&fs_cfg));

  for(PlotTable::iterator it=g_plot_table.begin();
      it != g_plot_table.end();
      it++)
    {
      char fname[128];
      sprintf(fname, "%s/%s%s", dir, prefix, it->first.c_str());
      printf("%s\n",fname);
      it->second->open(fname,mode);
    }
}

void
finalize_plot_table()
{
  for (PlotTable::iterator it=g_plot_table.begin();
       it != g_plot_table.end();
       it++)
    it->second->close();
}


void symcfg(int D,int P,int K,int M,int N,int T,
	    char *prefix,
	    Layer &F,
	    Probability &prob_Q,
	    Probability &accm_Q,
	    double &alpha,
	    Patterns &all_v,
	    double &connection_scale,
	    std::ostream& os = std::cout)
{
  using namespace std;
  for (size_t i=0; i<all_v.size(); i++)
    os << _eq(i) <<" "<< _eq(all_v[i]) << endl;
  _p(os,prob_Q);
  _p(os,accm_Q);
  _p(os,connection_scale);
  _p(os,D);
  _p(os,P);
  _p(os,K);
  _p(os,M);
  _p(os,N);
  _p(os,T);
  _p(os,g_energy_idx);
  _p(os,g_denergy_idx);
  _p(os,g_prop_idx);
  _p(os,g_annealing_schedule_idx);
  _p(os,g_init_layer_idx);
  _p(os,alpha);
  _p(os,F);
}
// for (size_t i=0; i<all_v.size(); i++)
//     os << _eq(i) <<" "<< _eq(all_v[i]) << endl;
//   gibbs_P = zerov(0x01<<D);
//   _p(os, prob_Q);
//   // plot a data dependent expection 
//   for (size_t n=0; n < N; n++){
//     size_t k;
//     for (k=0; k<accm_Q.size(); k++)
//       if (n < N*accm_Q[k])
// 	break;
//     data += ub::outer_prod(all_v[k], all_v[k]);
//   }
//   _p(os,data);


void symcfg(int D,int P,int K,int M,int N,int T,
	    char *prefix,
	    Layer &F,
	    SubLayer &L,
	    SubLayer &W,
	    SubLayer &J,
	    Patterns &vf,
	    Patterns &hf,
	    Probability &prob_Q,
	    Probability &accm_Q,
	    double &alpha,
	    Patterns &all_v,
	    Patterns &all_h,
	    Patterns &all_mu,
	    double &connection_scale,
	    std::ostream& os = std::cout)
{
  using namespace std;
  symcfg(D, P, K, M, N, T,
	 prefix,
	 F,
	 sF,
	 prob_Q,
	 accm_Q,
	 alpha,
	 all_v,
	 connection_scale,
	 os = std::cout);
  _p(os,L);
  _p(os,W);
  _p(os,J);
  for (size_t i=0; i<all_h.size(); i++)
    _p(os,all_h[i]);
  for (size_t i=0; i<all_mu.size(); i++)
    _p(os,all_mu[i]);
}

void
cook_args(int argc, char **argv,
	  int& D, int& P, int& K, int& M, int& N, int& T,
	  double& connection_scale,
	  char prefix[], char odir[])
{
  int result;
  char *tp;
  while ((result = getopt(argc,argv, "D:P:K:M:N:T:a:e:t:u:d:p:m:"))!= -1){
    switch (result){
    case 'D':
      D = atoi(optarg);
      break;
    case 'P':
      P = atoi(optarg);
      break;
    case 'K':
      K = atoi(optarg);
      break;
    case 'M':
      M = atoi(optarg);
      break;
    case 'N':
      N = atoi(optarg);
      break;
    case 'T':
      T = atoi(optarg);
      break;
    case 'a':
      g_annealing_schedule_idx = atoi(optarg);
      break;
    case 'e':
      g_energy_idx = atoi(optarg);
      break;
    case 'u':
      g_denergy_idx = atoi(optarg);
      break;
    case 't':
      g_prop_idx = atoi(optarg);
      break;
    case 'd':
      strcpy(odir  ,optarg);
      break;
    case 'p':
      strcpy(prefix,optarg);
      break;
    case 'm': 
      tp=strtok(optarg,",");
      g_init_layer_idx = atoi(tp);
      tp=strtok(NULL,",");
      connection_scale = atof(tp);
      break;
    case '?':
      std::cout << "optarg=" << optarg << std::endl;
      assert(false);
    }
  }

  _p(std::cout,D);
  _p(std::cout,P);
  _p(std::cout,K);
  _p(std::cout,M);
  _p(std::cout,N);
  _p(std::cout,T);
  _p(std::cout,odir);
  _p(std::cout,g_annealing_schedule_idx);
  _p(std::cout,prefix);
  _p(std::cout,g_energy_idx);
  _p(std::cout,g_denergy_idx);
}
