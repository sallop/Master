#include "FBM.hh"

size_t g_energy_idx ;
size_t g_denergy_idx;
size_t g_prop_idx   ;
size_t g_annealing_schedule_idx;
size_t g_init_layer_idx;

PlotTable g_plot_table;
// [kurata's matrix _energy3, _energy8] : [date's matrix _energy0]
// [kurata's matrix _denergy1, _energy4] : [date's matrix _denergy0]
energy_ptr g_energy_tbl[]={
  _energy0,// threshold is L(0,i) and L(i,0). date's impliments
  _energy1,// simple inner prod and then 0.5 times
  _energy2,// for loop version inner prod
  _energy3,// correct code at Kurata's matrix. diag is threshold.
  _energy4,// upper triangle elements
  _energy5,// (all elements - diagonal elements)*0.5
  _energy6,// (inner_prod - diagonal)*0.5
  _energy7,// (all elements - diagonal_const)*0.5
  _energy8,// correct code upper triangle matrix. same way energy3
  _energy9,// threshold is const in a Kurata's report [boltzmann machine]. 
  _energy10,// threshold is one side judge.
};

denergy_ptr g_denergy_tbl[]= {
  _denergy0,// linear sum + threshold_const. date's impliments.
  _denergy1,// collect linear sum - self input + threshold_const
  _denergy2,// linear sum - self inpu
  _denergy3,// linear sum
  _denergy4,// collect linear sum - self input + threshold_const
};

prop_ptr g_prop_tbl[] = {_prop0,_prop1,_prop2,_prop3,_prop4,};

annealing_ptr g_annealing_schedule_tbl[]={
  _annealing_schedule0,		// const value 0.01
  _annealing_schedule1,		// 1./(time + 1.)
  _annealing_schedule2,		// __C/(time + __D)
  _annealing_schedule3,		// _D/log(time)
};

init_layer_ptr g_init_layer_tbl[]={
  _init_layer0,			// matrix elements is 0 initialize
  _init_layer1,			// random*scale random domain is [0,1)
  _init_layer2,			// random*scale random. diagonal is zero 
  _init_layer3,			// const value init. diagonal is const
  _init_layer4,			// const value init. diagonal is zero
  _init_layer5,			// random*scale random domain is[-0.5,0.5)
  _init_layer6,			// random*scale diag is zero.
};

#define UDITER(W)				\
  for (size_t i=0; i<(W).size1(); i++)		\
    for (size_t j=i; j<(W).size2(); j++)

#define UITER(W)				\
  for (size_t i=0; i<(W).size1(); i++)		\
    for (size_t j=i+1; j<(W).size2(); j++)

Layer _init_layer0(int D, double scale)
{
  Layer L(D,D);
  UDITER(L){ L(i,j) = L(j,i) = 0.0; }
  return L;
}

Layer _init_layer1(int D, double scale)
{
  Layer L(D,D);
  UDITER(L){L(i,j) = L(j,i) = scale*drand48();}
  return L;
}

Layer _init_layer2(int D, double scale)
{
  Layer L = _init_layer2(D,scale);
  for (size_t i=0; i<D; i++) L(i,i)=0.0;
  return L;
}

Layer _init_layer3(int D, double scale)
{
  Layer L(D,D);
  UDITER(L){L(i,j) = L(j,i) = scale;}
  return L;
}

Layer _init_layer4(int D, double scale)
{
  Layer L = _init_layer3(D,scale);
  for (size_t i=0; i<D; i++) L(i,i) = 0.0;
  return L;
}

Layer _init_layer5(int D, double scale)
{
  Layer L(D,D);
  UDITER(L){ L(i,j) = L(j,i) = scale*(drand48()-0.5);}
  return L;
}

Layer _init_layer6(int D, double scale)
{
  Layer L = _init_layer5(D,scale);
  for (size_t i=0; i<D; i++){ L(i,i) = 0.0; }
  return L;
}

Layer init_layer(int D, double scale)
{
  static init_layer_ptr fn = g_init_layer_tbl[g_init_layer_idx];
  return fn(D,scale);
}

double energy(Pattern& v, Layer& L)
{
  static energy_ptr fn = g_energy_tbl[g_energy_idx];
  return fn(v,L);
}

double denergy(size_t i, Pattern &v, Layer &L)
{
  static denergy_ptr fn = g_denergy_tbl[g_denergy_idx];
  return fn(i,v,L);
}

double prop(int i, Pattern &v, Layer &L, double temperature)
{
  static prop_ptr fn = g_prop_tbl[g_prop_idx];
  return fn(i,v,L,temperature);
}

double annealing_schedule(double time)
{
  static annealing_ptr fn = g_annealing_schedule_tbl[g_annealing_schedule_idx];
  return fn(time);
}

void pvec(Pattern v, std::ostream& os = std::cout)
{
  using namespace std;
  for (size_t i=0; i<v.size(); i++)
    os << setprecision(8) << SETW << v[i] << " ";
  os << "\n";
}

inline double trace(Layer& W)
{
  using namespace boost::numeric::ublas;
  size_t D = W.size1();
  matrix_vector_range<Layer> diag(W,range(0,D),range(0,D));
  return sum(diag);
}

void pprobability(Probability& Pr, std::ostream& os=std::cout)
{
  for(size_t i=0; i<Pr.size(); i++)
    os << Pr[i] << " ";
  os << std::accumulate(Pr.begin(), Pr.end(), 0.0) << std::endl;
}

Pattern dtob(int n, size_t size)
{
  Pattern ret(size, 0);
#ifdef DATEIMP
  for (size_t i=0; i < size-1; i++)
    ret[size-1-i] = (n>>i) & 0x01;
    //ret[size-i] = (n>>(i-1)) & 0x01;
  //  std::cout << "DATEIMPLIMENT" << std::endl;
#else
  for (size_t i=0; i < size; i++)
    ret[size-1-i] = (n>>i) & 0x01;
  //  std::cout << "MACRO NOT INFARRENCE" << std::endl;
#endif
  return ret;
}

int btod(const Pattern& x)
{
  size_t i, ret=0, size=x.size();
#ifdef DATEIMP
  for(i=0; i < size-1; i++)
    ret += x[size-1-i]*(0x01<<i);
    //ret += x[size-i]*(0x01<<(i-1));
#else
  for(i=0; i < size; i++)
    ret += x[size-1-i]*(0x01 <<i);
#endif
  return ret;
}

double sigma(double x){return 1./(1+exp(-x));}

double _energy0(Pattern& v, Layer& L)
{
  double ene=0.0;
  for (size_t i=1; i<L.size1(); i++){
    ene += v[i]*L(0,i);
    for (size_t j=i+1; j<L.size2(); j++)
      ene += v[i]*L(i,j)*v[j];
  }
    
  return -ene;
}

double _energy1(Pattern& v, Layer& L)
{
  using namespace boost::numeric::ublas;
  double term = inner_prod(prod(v,L),v);
  return -0.5*term;
}

double _energy2(Pattern& v, Layer& L)
{
  double val2 = 0.0;
  for (size_t i=0; i<L.size1(); i++)
    for (size_t j=0; j<L.size2(); j++)
      val2 += v[i]*L(i,j)*v[j];
  return -0.5*val2;
}

// collect! for Kurata's model. result is "energy_cmp.cc"
double _energy3(Pattern& v, Layer& L)
{
  double val3 = 0.0;
  for (size_t i=0; i<L.size1(); i++)
    for (size_t j=i; j<L.size2(); j++)
      val3 += v[i]*L(i,j)*v[j];
  return -val3;
}

double _energy4(Pattern& v, Layer& L)
{
  double val4 = 0.0;
  for (size_t i=0; i<L.size1(); i++)
    for (size_t j=i+1; j<L.size2(); j++)
      val4 += v[i]*L(i,j)*v[j];
  return -val4;
}

double _energy5(Pattern& v, Layer& L)
{
  double val5 = 0.0;
  for (size_t i=0; i<L.size1(); i++){
    for (size_t j=0; j<L.size2(); j++){
      val5 += v[i]*L(i,j)*v[j];
    }
    val5 -= v[i]*L(i,i)*v[i];
  }
  return -0.5*val5;
}

double _energy6(Pattern& v, Layer& L)
{
  double val6 = 0.0;
  val6 = inner_prod(prod(v,L),v);
  val6 -= trace(L);
  return -0.5*val6;
}

double _energy7(Pattern& v, Layer& L)
{
  double val7 = 0.0;
  for (size_t i=0; i<L.size1(); i++){
    for (size_t j=0; j<L.size2(); j++){
      val7 += v[i]*L(i,j)*v[j];
    }
    val7 -= L(i,i);
  }
  return -0.5*val7;
}

// collect! for Kurata's model. result is "energy_cmp.cc"
double _energy8(Pattern& v, Layer& L)
{
  double ret = 0.0;
  for (size_t i=0; i<L.size1(); i++)
    for (size_t j=i; j<L.size2(); j++)
      ret += v[i]*L(i,j)*v[j];
  return -ret;
}

// I think this is hold water.
// Kurata's Report 「ボルツマン・マシン」 describe this
double _energy9(Pattern& v, Layer& L)
{
  double ret=0.0;
  for (size_t i=0; i<L.size1(); i++){
    for (size_t j=i+1; j<L.size2(); j++){
      ret += v[i]*L(i,j)*v[j];
    }
    ret += L(i,i);
  }
  return -ret;
}

// I think this is hold water 2. 
// 「ニューロンコンピューティングの基礎理論」
double _energy10(Pattern& v, Layer& L)
{
  double ret=0.0;
  for (size_t i=0; i<L.size1(); i++){
    ret += v[i]*L(i,i);
    for (size_t j=i+1; j<L.size2(); j++)
      ret += v[i]*L(i,j)*v[j];
  }
  return -ret;
}

double energy_test(Pattern& v, Layer& L, std::ostream &os=std::cout)
{
  _p(os,_energy0(v,L)); _p(os,_energy1(v,L)); _p(os,_energy2(v,L));
  _p(os,_energy3(v,L)); _p(os,_energy4(v,L)); _p(os,_energy5(v,L));
  _p(os,_energy6(v,L)); _p(os,_energy7(v,L)); _p(os,_energy8(v,L));
  _p(os,_energy9(v,L)); _p(os,_energy10(v,L));
  return energy(v,L);
}

double _denergy0(size_t i, Pattern& v, Layer& L)
{
  double u_i=0.0;
  for (size_t j=1; j<v.size(); j++)
    u_i += L(i,j)*v[j];
  return u_i + L(0,i);
}


// collect! result is "energy_cmp.cc"
double _denergy1(size_t i, Pattern& v, Layer& L)
{
  double u_i=0.0;
  for (size_t j=0; j<v.size(); j++)
    u_i += L(i,j)*v[j];
  return u_i - L(i,i)*v[i] + L(i,i);
}

double _denergy2(size_t i, Pattern& v, Layer& L)
{
  double u_i=0.0;
  for (size_t j=0; j<v.size(); j++)
    u_i += L(i,j)*v[j];
  return u_i - L(i,i)*v[i];
}

double _denergy3(size_t i, Pattern& v, Layer& L)
{
  double u_i=0.0;
  for (size_t j=0; j<v.size(); j++)
    u_i += L(i,j)*v[j];
  return u_i;
}

// collect! for Kurata's model. result is "energy_cmp.cc"
double _denergy4(size_t i, Pattern& v, Layer& L)
{
  double u_i=0.0;
  for (size_t j=0; j<i; j++)
    u_i += L(i,j)*v[j];
  for (size_t j=i+1; j<v.size(); j++)
    u_i += L(i,j)*v[j];
  u_i += L(i,i);// bias term
  return u_i;
}

double denergy_test(size_t i, Pattern& v, Layer& L,
		    std::ostream &os=std::cout)
{
  _p(os,_denergy0(i,v,L));
  _p(os,_denergy1(i,v,L));
  return denergy(i,v,L);
}

double _prop0(int i, Pattern& v, Layer& L, double temperature)
{
  return sigma(_denergy0(i,v,L)/temperature);
}

double _prop1(int i, Pattern& v, Layer& L, double temperature)
{
  return sigma(_denergy1(i,v,L)/temperature);
}

double _prop2(int i, Pattern& v, Layer& L, double temperature)
{
  return sigma(_denergy2(i,v,L)/temperature);
}

double _prop3(int i, Pattern& v, Layer& L, double temperature)
{
  return sigma(_denergy3(i,v,L)/temperature);
}

double _prop4(int i, Pattern& v, Layer& L, double temperature)
{
  return sigma(_denergy4(i,v,L)/temperature);
}

double _annealing_schedule0(double time)
{				// const
#ifdef TEMPERATURE
  return TEMPERATURE;
#else
  return 0.01;
#endif
}

double _annealing_schedule1(double time)
{
  return 1.0/(time + 1);
}

double _annealing_schedule2(double time)
{
#ifndef __C			// Largest Energy Barrier
#define __C 5.0
#endif
#ifndef __D			// usually set equal to one
#define __D 1.0
#endif
  return __C/(time + __D);
#undef __C
#undef __D
}

double _annealing_schedule3(double time)
{
#ifndef _D
#define _D (0.1)
#endif
  return _D/log(time);
#undef _D
}

double prop_test(int i, Pattern& v, Layer& L,
		 double temperature,std::ostream& os=std::cout)
{
  _p(os,_prop0(i,v,L,temperature));
  _p(os,_prop1(i,v,L,temperature));
  _p(os,_prop2(i,v,L,temperature));
  _p(os,_prop3(i,v,L,temperature));
  return prop(i,v,L,temperature);
}

// Z(\theta)
double partition_function(Layer& L, Patterns& all_v)
{
  double ret=0.0;
  for(Patterns::iterator v = all_v.begin(); v != all_v.end(); ++v){
    ret += exp(-energy(*v,L));
  }
  return ret;
}

//double boltzmann_distribution(Pattern& v,
double boltzmann_distribution(Pattern& v,Layer& L, Patterns& all_v)
{
  double numer = exp(-energy(v,L));
  double denom = partition_function(L,all_v);
  return numer/denom;
}

double KLd(Probability& p, Probability& q)
{
  size_t i, size = (p.size() == q.size()) ? p.size() : 0;
  double d=0.0;
  double term1, term2;
  //for (i=0; i < size; i++){d += p[i]*log(p[i]/q[i]);}
  for (i=0; i < size; i++){
    term1 = p[i]*log(p[i]);
    term2 = p[i]*log(q[i]);
    if (!isnan(term1)){ d += term1;}
    if (!isnan(term2)){ d -= term2;}
    //d += p[i]*log(p[i]) - p[i]*log(q[i]);
  }
  return d;
}

inline void kstep_gibbs(size_t K, Pattern& v, Layer& L, double temperature)
{
  size_t D = v.size();
  for(size_t k=0; k < K; k++){
    size_t i;
#ifdef DATEIMP
    i = (rand()%(D-1)) + 1;
#else
    i = rand()%D;
#endif
    v[i] = (prop(i, v, L, temperature) > drand48()) ? ONBIT : OFFBIT;
  }
}



void init_plot_table(char *dir, char *prefix,
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
      char fname[32];
      sprintf(fname, "%s/%s%s", dir, prefix, it->first.c_str());
      printf("%s\n",fname);
      it->second->open(fname,mode);
    }
}

void finalize_plot_table()
{
  for (PlotTable::iterator it=g_plot_table.begin();
       it != g_plot_table.end();
       it++)
    {
      it->second->close();
    }
}


void symcfg(int D,int P,int K,int M,int N,int T,
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
	    std::ostream& os = std::cout)
{
  using namespace std;
  double temperature = 1.0;

  _p(os,D);_p(os,P);_p(os,K);_p(os,M);_p(os,N);_p(os,T);
  _p(os,g_energy_idx);
  _p(os,g_denergy_idx);
  _p(os,g_prop_idx);
  _p(os,g_annealing_schedule_idx);
  _p(os,g_init_layer_tbl);
  _p(os,alpha);
  _p(os,L);_p(os,sL);_p(os,bL);
  _p(os,prob_Q);_p(os,accm_Q);
  _p(os,connection_scale);

  for (size_t i=0; i<all_v.size(); i++)
    os << _eq(i) <<" "<< _eq(all_v[i]) << endl;

  for (int n=0; n < 0x01<<D; ++n){
    Pattern& vn = all_v[n];
    boltz_P[n] = boltzmann_distribution(vn, L, all_v);
    os << _eq(n)  << " "
       << _eq(vn) << " "
       << _eq(energy(vn,L)) << " "
       << _eq(boltz_P[n])   << "\n";
  }

  _p(os, partition_function(L, all_v));

  gibbs_P = zerov(0x01<<D);
  _p(os, prob_Q);
  for (int n=0; n<0x01<<D; n++){
    vf = all_v[n];
    os << _eq(vf) << " ";
    for (size_t m=0; m<M; m++){
      kstep_gibbs(K,vf,L,temperature);
      gibbs_P[btod(vf)]++;
    }
    gibbs_P /= M;
    os << _eq(gibbs_P) << "\n";
  }

  for (int n=0; n<0x01<<D; n++){
    vf = all_v[n];
    os << _eq(n) << " " << _eq(vf) << endl;
    for (size_t i=0; i<D; i++){
      os << _eq(i) << " "
	 << _eq(prop(i,vf,L,temperature))
	 << std::endl;
    }
    os << endl;
  }
  
  {//plot kl divergence
    os << "KL-divergence gibbs-sampling" << endl;
    _p(os, KLd(gibbs_P,gibbs_P));
    _p(os, KLd(prob_Q , prob_Q));
    _p(os, KLd(gibbs_P, prob_Q));
    _p(os, KLd(prob_Q ,gibbs_P));
    os << "KL-divergence boltzmann-dist" << endl;
    _p(os, KLd(boltz_P,boltz_P));
    _p(os, KLd(prob_Q , prob_Q));
    _p(os, KLd(boltz_P, prob_Q));
    _p(os, KLd(prob_Q ,boltz_P));
  }

  {// data dependent expection
    for (size_t n=0; n < N; n++){
      size_t k;
      for (k=0; k<accm_Q.size(); k++)
	if (n < N*accm_Q[k])
	  break;
      data += ub::outer_prod(all_v[k], all_v[k]);
    }
    _p(os,data);
  }
}

void cook_args(int argc, char **argv,
	       int& D, int& P, int& K, int& M, int& N, int& T,
	       double& connection_scale,
	       char prefix[], char odir[])
{
  int result;
  while ((result = getopt(argc,argv, "D:P:K:M:N:T:a:e:u:d:p:m:"))
	 != -1){
    switch (result){
    case 'D': D = atoi(optarg); break;
    case 'P': P = atoi(optarg); break;
    case 'K': K = atoi(optarg); break;
    case 'M': M = atoi(optarg); break;
    case 'N': N = atoi(optarg); break;
    case 'T': T = atoi(optarg); break;
    case 'a': g_annealing_schedule_idx = atoi(optarg);
    case 'e': g_energy_idx = atoi(optarg); break;
    case 'u': g_denergy_idx = atoi(optarg); break;
    case 'd': strcpy(odir  ,optarg); break;
    case 'p': strcpy(prefix,optarg); break;
    case 'm':
      g_init_layer_idx = optarg[0] - '0';
      connection_scale = atof(&optarg[1]);
      std::cout << _eq(g_init_layer_idx) << std::endl;
      std::cout << _eq(connection_scale) << std::endl;
      assert(false);
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
