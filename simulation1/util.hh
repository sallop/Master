#ifndef UTIL_H_
#define UTIL_H_

#include "BM.hh"
#include <iostream>
#include <unistd.h>
#include <cmath>
#include <cstring>
#include <cfloat>
#include <cstdlib>
#include <cassert>
#include <time.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <utility>
#include <vector>
#include <map>
#include <utility>
#include <functional>
#include <numeric>

#define SETW std::setw(16)
#define _eq(val) std::setiosflags(std::ios::left) << #val << "="	\
  << std::setw(16) << std::setiosflags(std::ios::left) << val
#define _p(os,val) os << _eq(val) << std::endl

typedef std::pair<std::string, std::fstream*> PlotRecord;
typedef std::map <std::string, std::fstream*> PlotTable ;

extern PlotTable g_plot_table;

void cook_args(int argc, char **argv,
	       int &D, int &P, int &K, int &M, int &N, int &T,
	       double &connection_scale,
	       char *prefix, char *odir);

void symcfg(int D, int P,int K,int M,int N,int T,
	    char *prefix,
	    Layer &L, 
	    Layer &data,
	    Layer &model,
	    symmetric_view &sL,
	    symmetric_view &sdata,
	    symmetric_view &smodel,
	    banded_view &bL,
	    banded_view &bdata,
	    banded_view &bmomdel,
	    Pattern &v,
	    Pattern &vf,
	    Probability &boltz_P,
	    Probability &gibbs_P,
	    Probability &prob_Q,
	    Probability &accm_Q,
	    double &alpha,
	    double &kl,
	    Patterns &all_v,
	    double &connection_scale,
	    std::ostream &os);

// using GBM
void symcfg(int D,int P,int K,int M,int N,int T,
	    char *prefix,
	    Layer &F,
	    Layer &L,Layer &Ld,Layer &Lm,
	    Layer &W,Layer &Wd,Layer &Wm,
	    Layer &J,Layer &Jd,Layer &Jm,
	    Pattern &v, Pattern &vf,
	    Pattern &h, Pattern &hf,
	    Probability &boltz_P,
	    Probability &gibbs_P,
	    Probability &prob_Q,
	    Probability &accm_Q,
	    double &alpha,
	    double &kl,
	    Patterns &all_v,
	    Patterns &mu,
	    double &connection_scale,
	    std::ostream& os);

void init_plot_table(char *odir,
		     char *prefix,
		     std::fstream &fs_mat,
		     std::fstream &fs_t_gibbs,
		     std::fstream &fs_cfg,
		     std::ios_base::openmode mode);
		     //= std::ios_base::in|std::ios_base::out);

void finalize_plot_table();

#endif	// UTIL_H_
