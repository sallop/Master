#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <utility>
#include <iostream>
#include <boost/numeric/ublas/io.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/regex.hpp>
#include "Layer.hh"
#include "BM.hh"
#include "config.hh"

extern PlotTable g_plot_table;

using namespace std;
using namespace boost::algorithm;
using namespace boost;

#define BUFSZ 256

string&
trim(string& s, const char* t=" \t\n\r\f\v")
{
  s.erase(0,s.find_first_not_of(t));
  s.erase(s.find_last_not_of(t) + 1);
  return s;
}

string
read_param(const string& line, string param)
{
  size_t s = line.find(param), e;
  if (string::npos == s) return "";
  s = s + param.size();
  e = line.size();
  string subline = line.substr(s,e);
  return trim(subline);
}

// this function need to call before init_plot_table();
void
readcfg(size_t &D,size_t &P,size_t &K,size_t &M,size_t &N,size_t &T,
	Probability &prob_Q,
	Probability &accm_Q,
	Probability &boltzP,
	double &connection_scale,
	Patterns &all_v,
	Patterns &all_h,
	Patterns &all_mu,
	double &alpha,
	Layer &F)
{
  typedef vector<string> strings;
  typedef map<string,strings> config_table;
  config_table cfg_tbl;
#define SET(x) cfg_tbl.insert(make_pair((x), strings()))
  SET("all_v[i]=");
  SET("prob_Q=");
  SET("accm_Q=");
  SET("connection_scale=");
  SET("D=");
  SET("P=");
  SET("K=");
  SET("M=");
  SET("N=");
  SET("T=");
  SET("g_energy_idx=");
  SET("g_denergy_idx=");
  SET("g_prop_idx=");
  SET("g_annealing_schedule_idx=");
  SET("g_init_layer_idx=");
  SET("alpha=");
  SET("F=");
  SET("all_h[i]=");
  SET("all_mu[i]=");
#undef SET
  
#define PLOTTBL(){					\
    for (config_table::iterator it=cfg_tbl.begin();	\
	 it != cfg_tbl.end(); it++){			\
      cout << "(key = " << it->first << ", value = ";	\
      for(strings::iterator it2 = it->second.begin();	\
	  it2 != it->second.end(); it2++){		\
	cout << *it2 << ",";				\
      }							\
      cout << ")" << endl;				\
    }							\
  }
  
  cout << "-------before-------" << endl;
  PLOTTBL();
  fstream &ifs = *g_plot_table["-cfg.dat"];
  string line, tmp;
  istringstream iss;
  for (size_t nline=0; getline(ifs,line); nline++){
    for (config_table::iterator it = cfg_tbl.begin();
	 it != cfg_tbl.end(); it++){
      tmp = read_param(line,it->first);
      if (tmp != ""){
	it->second.push_back(tmp);
	break;
      }
    }
  }
  cout << "-------after-------" << endl;
  PLOTTBL();
#undef PLOTTBL
#define GET(x) cfg_tbl[(x)][0].c_str()
  D = atoi(GET("D="));
  P = atoi(GET("P="));
  K = atoi(GET("K="));
  M = atoi(GET("M="));
  N = atoi(GET("N="));
  T = atoi(GET("T="));
  g_energy_idx  = atoi(GET("g_energy_idx=")); // -> error
  g_denergy_idx = atoi(GET("g_denergy_idx="));
  g_prop_idx    = atoi(GET("g_prop_idx="));
  g_annealing_schedule_idx = atoi(GET("g_annealing_schedule_idx="));
  iss.str(GET("prob_Q="));
  iss >> prob_Q;
  
  D = atoi(GET("D="));
  P = atoi(GET("P="));
  K = atoi(GET("K="));
  M = atoi(GET("M="));
  N = atoi(GET("N="));
  T = atoi(GET("T="));
  g_energy_idx		   = atoi(GET("g_energy_idx="));
  g_denergy_idx		   = atoi(GET("g_denergy_idx="));
  g_prop_idx		   = atoi(GET("g_prop_idx="));
  g_annealing_schedule_idx = atoi(GET("g_annealing_schedule_idx="));
  g_init_layer_idx	   = atoi(GET("g_init_layer_idx="));

  alpha = atof(GET("alpha="));
  connection_scale = atof(GET("connection_scale="));

  iss.str(GET("F="))     ; iss >> F;
  iss.str(GET("prob_Q=")); iss >> prob_Q;
  iss.str(GET("accm_Q=")); iss >> accm_Q;

  for (size_t i=0; i < cfg_tbl["all_v[i]="].size(); i++){
    Pattern v;
    iss.str(cfg_tbl["all_v[i]="][i]); iss >> v;
    all_v.push_back(v);
  }
  
  for (size_t i=0; i < cfg_tbl["all_h[i]="].size(); i++){
    Pattern h;
    iss.str(cfg_tbl["all_h[i]="][i]); iss >> h;
    all_h.push_back(h);
  }

  for (size_t i=0; i < cfg_tbl["all_mu[i]="].size(); i++){
    Pattern mu;
    iss.str(cfg_tbl["all_mu[i]="][i]); iss >> mu;
    all_mu.push_back(mu);
  }
  boltzP.resize(all_v.size());
#undef GET
}

void
read_t_gibbs(char *buf, size_t size, double &t, Probability &prob)
{
  istringstream iss(buf, istringstream::in);
  iss >> t >> prob;
}

void
read_t_gibbs(const string& str, double &t, Probability &prob)
{
  istringstream iss(str, istringstream::in);
  iss >> t >> prob;
}

void
read_t_mat(const string& str, double &t, Layer &L)
{
  istringstream iss(str, istringstream::in);
  iss >> t >> L;
}

int
main(int argc, char *argv[])
{
  char dir[BUFSZ]="dir-dat", prefix[BUFSZ]="g1";
  bool output_flg = false;
  fstream fs_mat, fs_t_gibbs, fs_cfg;
  double pause_time = 4.2;
  int result;
  while ((result = getopt(argc,argv,"d:p:o")) != -1){
    char *tp;
    switch (result){
    case 'd': strcpy(dir, optarg)   ; break;
    case 'p': strcpy(prefix, optarg); break;
    case 'o':
      output_flg = true;
      pause_time = 0.0;
      break;
    case '?':
      cout << "optarg =" << optarg << endl;
      assert(false);
    }
  }

  Patterns all_v, all_h;
  Probability gibbsP, boltzP;
  Probability probQ, accmQ;
  double connection_scale;
  size_t D, P, K, M, N, T;
  double alpha;
  Layer F;//,L,W,J;
  Patterns all_mu;

  init_plot_table(dir, prefix, fs_mat, fs_t_gibbs, fs_cfg, ios_base::in);
  
  for (PlotTable::iterator it=g_plot_table.begin();
       it != g_plot_table.end(); it++){
    cout << "(key=" << it->first << ",value=" << it->second << ")" << endl;
  }

  readcfg(D, P, K, M, N, T,
	  probQ, accmQ, boltzP,
	  connection_scale,
	  all_v,all_h,all_mu,
	  alpha,
	  F);

  SubLayer L(F,ub::range(0,D)  ,ub::range(0,D));
  SubLayer W(F,ub::range(0,D)  ,ub::range(D,D+P));
  SubLayer J(F,ub::range(D,D+P),ub::range(D,D+P));

#define TQ	"Q"
#define TPg	"P_g"
#define TPb	"P_b"
#define TKL(P) str(KL(P||Q))
#define TE(n)	 str(energy##n)
#define str(s)   #s
#define EACHA(ary) for(size_t i=0; i<sizeof((ary))/sizeof((ary)[0]); i++)
#define EACHV(vec) for(size_t i=0; i<(vec).size(); i++)
#define RANGE(start,end) for(size_t i=(start); i<(end); i++)
#define GPPLOT(slot,...) fprintf(gtbl[(slot)].fp, __VA_ARGS__)
#define FPPLOT(slot,...) fprintf(ftbl[(slot)].fp, __VA_ARGS__)
  struct {
    FILE *fp; const char *id;
  } gtbl[] = {
    { NULL, "",},	                      // 0
    { popen("gnuplot","w"), "boltzP_t_kl",},  // 1
    { NULL, "",},                             // 2
    { popen("gnuplot","w"), "probQ_boltzP",}, // 3
  },ftbl[] = {
    { NULL, "t_boltzP",},	     // 0
    { NULL, "t_energy",},	     // 1
    { NULL, "gibbsP_t_kl",},  // 2 add 2011/12/26
    { NULL, "boltzP_t_kl",}, 	     // 3 add 2011/12/26
  };

  EACHA(ftbl){
    char filename[BUFSZ];
    sprintf(filename, "%s/%s-%s.dat", dir, prefix, ftbl[i].id);
    ftbl[i].fp = fopen(filename,"w");
  }

  EACHA(gtbl){
    if(gtbl[i].fp != NULL){
      if (output_flg){
	GPPLOT(i,"set terminal postscript enhanced color eps\n");
	GPPLOT(i,"set output '%s/%s-%s.eps'\n", dir, prefix, gtbl[i].id);
      }

      GPPLOT(i,"set style line 1 lt 1 lw 2 lc 1\n"); // KL(gibbsP||Q)
      GPPLOT(i,"set style line 2 lt 1 lw 2 lc 3\n"); // KL(boltzP||Q)
      GPPLOT(i,"set style line 3 lt 1 lw 2 lc 4\n"); // prob_Q
      GPPLOT(i,"set style line 4 lt 1 lw 2 lc 5\n"); // gibbsP
      GPPLOT(i,"set style line 5 lt 1 lw 2 lc 8\n"); // boltzP
      GPPLOT(i,"set grid\n");
    }
  }  

  //  GPPLOT(1, "set yrange[0:.3]\n");
  GPPLOT(1, "plot '-' using 1:2 with lines title '%s' ls %d\n", TKL(P_b), 2);
  std::string line;
  
  // print KL divergence loop
  while (getline(fs_mat,line)){
    double t, kl_boltzP_probQ;
    read_t_mat  (line, t, F);
    double sum_energy = 0.0;
    kl_boltzP_probQ = KLd(boltzP, probQ);
    GPPLOT(1, "%lf %lf\n", t, kl_boltzP_probQ);// t_KL(boltzP||probQ)

    FPPLOT(0, "%lf ", t);// t_boltzP
    FPPLOT(1, "%lf ", t);// t_energy
    FPPLOT(3, "%lf %lf\n", t, kl_boltzP_probQ);
    //    fprintf(stdout, "%lf %lf\n", t, kl_boltzP_probQ);
    for(size_t i=0; i<all_v.size(); i++){
      Pattern &v = all_v[i];
      boltzP[i] = prob_v(v,L,J,W,all_v,all_h);
      FPPLOT(0, "%lf ", boltzP[i]);
      for (size_t j=0; j<all_h.size(); j++){
	Pattern &h = all_h[j];
	double j_energy = energy(v,h,L,J,W);
	sum_energy += j_energy;
	FPPLOT(1, "%lf ", j_energy);
      }
    }
    FPPLOT(0, "%lf\n", std::accumulate(boltzP.begin(),boltzP.end(),0.0));
    FPPLOT(1, "%lf\n", sum_energy);
  }
  GPPLOT(1,"e\n");
  GPPLOT(1,"pause %lf\n",pause_time);
  //  GPPLOT(3, "set yrange[0.0:1.0]\n");
  GPPLOT(3, "plot ");
  GPPLOT(3,
	 "'-' using ($1):2:(0.2) with boxes fs pattern 2 title '%s' ls %d,"
	 ,TQ ,3);
  GPPLOT(3,
	 "'-' using ($1+0.33):2:(0.2) with boxes fs pattern 2 title '%s' ls %d\n"
	 ,TPb,5);
  EACHV(probQ ){ GPPLOT(3,"%d %lf\n", i,  probQ[i]);} GPPLOT(3, "e\n");
  EACHV(boltzP){ GPPLOT(3,"%d %lf\n", i, boltzP[i]);} GPPLOT(3, "e\n");
  GPPLOT(3, "pause %lf\n",pause_time);
  pclose(gtbl[1].fp);
  pclose(gtbl[3].fp);
  EACHA(ftbl){fclose(ftbl[i].fp);}

  char ifname[sizeof(ftbl)/sizeof(ftbl[0])][BUFSZ];
  char ofname[sizeof(ftbl)/sizeof(ftbl[0])][BUFSZ];
  RANGE(0,2){
    ftbl[i].fp = popen("gnuplot","w");
    sprintf(ifname[i], "%s/%s-%s.dat", dir, prefix, ftbl[i].id);
    sprintf(ofname[i], "%s/%s-%s.eps", dir, prefix, ftbl[i].id);
    
    FPPLOT(i,"set style line 1 lt 1 lw 10 lc 1\n"); // KL(gibbsP||Q)
    FPPLOT(i,"set style line 2 lt 1 lw 10 lc 3\n"); // KL(boltzP||Q)
    FPPLOT(i,"set style line 3 lt 1 lw 10 lc 4\n"); // prob_Q
    FPPLOT(i,"set style line 4 lt 1 lw 10 lc 5\n"); // gibbsP
    FPPLOT(i,"set style line 5 lt 1 lw 10 lc 8\n"); // boltzP
    FPPLOT(i,"set style line 6 lt 3 lw 2 lc 3\n");  // sum Energy
    FPPLOT(i,"set grid\n");
    if (output_flg){
      FPPLOT(i,"set terminal postscript enhanced color eps\n");
      FPPLOT(i,"set output '%s'\n", ofname[i]);
    }
    FPPLOT(i,"plot ");
    switch (i){
      size_t j,k;
    case 0:
      for (j=0; j<all_v.size(); j++){
	FPPLOT(i,"'%s' u 1:%d with lines linewidth 2 title 'P(v^%d)',",
	       ifname[i], j+2, j);
      }
      FPPLOT(i,"'%s' u 1:%d with lines  linewidth 2 title '\\sum_n P(v^n)'\n",
	     ifname[i], j+2);
      break;
    case 1:
      for (j=0; j<all_v.size(); j++){
	for (k=0; k<all_h.size(); k++){
	  FPPLOT(i, "'%s' u 1:%d with lines linewidth 2 title 'En(v^%d,h^%d)',",
		 ifname[i],
		 j*all_h.size() + k + 2,
		 j,
		 k);
	}
      }
      FPPLOT(i,"'%s' u 1:%d with lines ls 6 title '\\sum_n En(v^n)'\n",
	     ifname[i], (j-1)*all_h.size() + k + 2);
      break;
    }
    FPPLOT(i,"pause %lf\n",pause_time);
    pclose(ftbl[i].fp);
  }

  finalize_plot_table();
  return 0;
}
