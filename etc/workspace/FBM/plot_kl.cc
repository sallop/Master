#include "FBM.hh"
#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <utility>
#include <iostream>
#include <boost/numeric/ublas/io.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/regex.hpp>

using namespace std;
using namespace boost::algorithm;
using namespace boost;

string& trim(string& s, const char* t=" \t\n\r\f\v")
{
  s.erase(0,s.find_first_not_of(t));
  s.erase(s.find_last_not_of(t) + 1);
  return s;
}

string read_param(const string& line, string param)
{
  size_t s = line.find(param), e;
  if (string::npos == s) return "";
  s = s + param.size();
  e = line.size();
  string subline = line.substr(s,e);
  return trim(subline);
}

// this function need to call before init_plot_table();
void readcfg(size_t& D,size_t& P,size_t& K,size_t& M,size_t& N,size_t& T,
	     Probability &boltzP, Probability &gibbsP, Probability &prob_Q,
	     double &alpha, double &kl, Patterns &all_v)
{
  typedef vector<string> strings;
  typedef map<string,strings> config_table;
  config_table cfg_tbl;
#define SET(x) cfg_tbl.insert(make_pair((x), strings()))
  // readable parameter
  SET("D=");SET("P=");SET("K=");SET("M=");SET("N=");SET("T=");
  SET("alpha="); SET("all_v[i]="); SET("prob_Q=");
  SET("g_energy_idx="); SET("g_denergy_idx=");
  SET("g_prop_idx="); SET("g_annealing_schedule_idx=");
#undef SET
#define PLOTTBL() {					\
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
  for (size_t nline=0; nline < 32 && getline(ifs,line); nline++){
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
  D = atoi(GET("D=")); P = atoi(GET("P=")); K = atoi(GET("K="));
  M = atoi(GET("M=")); N = atoi(GET("N=")); T = atoi(GET("T="));
  g_energy_idx  = atoi(GET("g_energy_idx=")); // -> error
  g_denergy_idx = atoi(GET("g_denergy_idx="));
  g_prop_idx    = atoi(GET("g_prop_idx="));
  g_annealing_schedule_idx = atoi(GET("g_annealing_schedule_idx="));
  iss.str(GET("prob_Q="));
  iss >> prob_Q;
  gibbsP.resize(prob_Q.size());
  boltzP.resize(prob_Q.size());
#undef GET
  Pattern v;
  for (size_t i=0; i < cfg_tbl["all_v[i]="].size(); i++){
    iss.str(cfg_tbl["all_v[i]="][i]);
    iss >> v;
    cout << _eq(v) << endl; 
    all_v.push_back(v);
  }
}

void read_t_gibbs(char *buf, size_t size, double &t, Probability &prob)
{
  istringstream iss(buf, istringstream::in);
  iss >> t >> prob;
}

void read_t_gibbs(const string& str, double &t, Probability &prob)
{
  istringstream iss(str, istringstream::in);
  iss >> t >> prob;
}

void read_t_mat(const string& str, double &t, Layer &L)
{
  istringstream iss(str, istringstream::in);
  iss >> t >> L;
}

int main(int argc, char *argv[])
{
  char dir[32]="dir-dat", prefix[32]="K";
  bool output_flg = false;
  fstream fs_mat, fs_t_gibbs, fs_cfg;
  double pause_time = 4.2;
  int result;
  while ((result = getopt(argc,argv,"d:p:o")) != -1){
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
  size_t D, P, K, M, N, T;
  Probability gibbsP, boltzP, probQ;
  double alpha, kl;
  double sum_energy;
  Patterns all_v;
  Layer L;

  init_plot_table(dir, prefix, fs_mat, fs_t_gibbs, fs_cfg, ios_base::in);
  
  for (PlotTable::iterator it=g_plot_table.begin();
       it != g_plot_table.end(); it++){
    cout << "(key=" << it->first << ",value=" << it->second << ")" << endl;
  }

  readcfg(D, P, K, M, N, T, boltzP, gibbsP, probQ, alpha, kl, all_v);
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
    { popen("gnuplot","w"), "gibbsP_t_kl",},	     // 0
    { popen("gnuplot","w"), "boltzP_t_kl",}, 	     // 1
    { popen("gnuplot","w"), "gibbsP_boltzP_t_kl",},  // 2
    { popen("gnuplot","w"), "probQ_gibbsP_boltzP",},	// 3
  },ftbl[] = {
    { NULL, "t_boltzP",},	     	// 0
    { NULL, "t_energy",},	     	// 1
  };

  EACHA(ftbl){
    char filename[32];
    sprintf(filename, "%s/%s-%s.dat", dir, prefix, ftbl[i].id);
    ftbl[i].fp = fopen(filename,"w");
  }

  EACHA(gtbl){
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

  GPPLOT(0, "plot '-' using 1:2 with lines title '%s' ls %d\n", TKL(P_g), 1);
  GPPLOT(2, "plot ");
  GPPLOT(2, "'-' using 1:2 with lines title '%s' ls %d, ", TKL(P_g), 1);
  GPPLOT(2, "'-' using 1:2 with lines title '%s' ls %d\n", TKL(P_b), 2);
  string line;
  while (getline(fs_t_gibbs,line)){
    double t, kl_gibbsP_probQ;
    read_t_gibbs(line, t, gibbsP);
    kl_gibbsP_probQ = KLd(gibbsP, probQ);
    GPPLOT(0, "%lf %lf\n", t, kl_gibbsP_probQ);
    GPPLOT(2, "%lf %lf\n", t, kl_gibbsP_probQ);
  }
  GPPLOT(0, "e\n");
  GPPLOT(2, "e\n");

  GPPLOT(1, "plot '-' using 1:2 with lines title '%s' ls %d\n", TKL(P_b), 2);
  while (getline(fs_mat,line)){
    double t, kl_boltzP_probQ;
    read_t_mat  (line, t, L);
    sum_energy = 0.0;
    kl_boltzP_probQ = KLd(boltzP, probQ);
    GPPLOT(1, "%lf %lf\n", t, kl_boltzP_probQ);
    GPPLOT(2, "%lf %lf\n", t, kl_boltzP_probQ);
    FPPLOT(0, "%lf ", t);
    FPPLOT(1, "%lf ", t);
    EACHV(all_v){
      double i_energy = energy(all_v[i],L);
      sum_energy += i_energy;
      boltzP[i] = boltzmann_distribution(all_v[i],L,all_v);
      FPPLOT(0, "%lf ", boltzP[i]);
      FPPLOT(1, "%lf ", i_energy );
    }
    FPPLOT(0, "%lf\n", std::accumulate(boltzP.begin(),boltzP.end(),0.0));
    FPPLOT(1, "%lf\n", sum_energy);
  }
  GPPLOT(1, "e\n");
  GPPLOT(2, "e\n");
  EACHA(gtbl){ GPPLOT(i,"pause %lf\n",pause_time);}

  //  GPPLOT(3, "set xrange[-0.5:%lf]\n",probQ.size()-1);
  GPPLOT(3, "set yrange[ 0.0:1.0]\n");
  GPPLOT(3, "plot ");
  GPPLOT(3, "'-' using ($1-0.33):2:(0.2) with boxes fs pattern 2 title '%s' ls %d," ,TQ ,3);
  GPPLOT(3, "'-' using ($1+0.00):2:(0.2) with boxes fs pattern 2 title '%s' ls %d," ,TPg,4);
  GPPLOT(3, "'-' using ($1+0.33):2:(0.2) with boxes fs pattern 2 title '%s' ls %d\n",TPb,5);
  EACHV(probQ ){ GPPLOT(3,"%d %lf\n", i,  probQ[i]);} GPPLOT(3, "e\n");
  EACHV(gibbsP){ GPPLOT(3,"%d %lf\n", i, gibbsP[i]);} GPPLOT(3, "e\n");
  EACHV(boltzP){ GPPLOT(3,"%d %lf\n", i, boltzP[i]);} GPPLOT(3, "e\n");
  GPPLOT(3, "pause %lf\n",pause_time);
  EACHA(gtbl){ pclose(gtbl[i].fp);}
  EACHA(ftbl){ fclose(ftbl[i].fp);}

  char ifname[sizeof(ftbl)/sizeof(ftbl[0])][64];
  char ofname[sizeof(ftbl)/sizeof(ftbl[0])][64];
  EACHA(ftbl){
    ftbl[i].fp = popen("gnuplot","w");
    sprintf(ifname[i], "%s/%s-%s.dat", dir, prefix, ftbl[i].id);
    sprintf(ofname[i], "%s/%s-%s.eps", dir, prefix, ftbl[i].id);
    
    FPPLOT(i,"set style line 1 lt 1 lw 10 lc 1\n"); // KL(gibbsP||Q)
    FPPLOT(i,"set style line 2 lt 1 lw 10 lc 3\n"); // KL(boltzP||Q)
    FPPLOT(i,"set style line 3 lt 1 lw 10 lc 4\n"); // prob_Q
    FPPLOT(i,"set style line 4 lt 1 lw 10 lc 5\n"); // gibbsP
    FPPLOT(i,"set style line 5 lt 1 lw 10 lc 8\n"); // boltzP
    FPPLOT(i,"set grid\n");
    if (output_flg){
      FPPLOT(i,"set terminal postscript enhanced color eps\n");
      FPPLOT(i,"set output '%s'\n", ofname[i]);
    }
    FPPLOT(i,"plot ");
    switch (i){
      size_t j;
    case 0:
      for (j=0; j<all_v.size(); j++)
	FPPLOT(i, "'%s' u 1:%d with lines title 'P(v^%d)',",ifname[i], j+2, j);
      FPPLOT(i, "'%s' u 1:%d with lines title '\\sum_n P(v^n)'\n",ifname[i], j+2);
      break;
    case 1:
      for (j=0; j<all_v.size(); j++)
	FPPLOT(i, "'%s' u 1:%d with lines title 'En(v^%d)',", ifname[i], j+2, j);
      FPPLOT(i, "'%s' u 1:%d with lines title '\\sum_n En(v^n)'\n", ifname[i], j+2, j);
      break;
    }
    FPPLOT(i,"pause %lf\n",pause_time);
    pclose(ftbl[i].fp);
  }

  finalize_plot_table();
  return 0;
}
