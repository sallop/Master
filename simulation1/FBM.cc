#include "BM.hh"
#include "FBM.hh"
//#include "util.hh"
#include "config.hh"

//PlotTable g_plot_table;
// k[urata's matrix _energy3, _energy8] : [date's matrix _energy0]
// [kurata's matrix _denergy1, _energy4] : [date's matrix _denergy0]
energy_ptr g_energy_tbl[]={
  _energy0 ,// threshold is L(0,i) and L(i,0). date's impliments
  _energy1 ,// simple inner prod and then 0.5 times
  _energy2 ,// for loop version inner prod
  _energy3 ,// correct code at Kurata's matrix. diag is threshold.
  _energy4 ,// upper triangle elements
  _energy5 ,// (all elements - diagonal elements)*0.5
  _energy6 ,// (inner_prod - diagonal)*0.5
  _energy7 ,// (all elements - diagonal_const)*0.5
  _energy8 ,// correct code upper triangle matrix. same way energy3
  _energy9 ,// threshold is const in a Kurata's report [boltzmann machine]. 
  _energy10,// threshold is one side judge.
};

denergy_ptr g_denergy_tbl[]= {
  _denergy0,// linear sum + threshold_const. date's impliments.
  _denergy1,// collect linear sum - self input + threshold_const
  _denergy2,// linear sum - self inpu
  _denergy3,// linear sum
  _denergy4,// collect linear sum - self input + threshold_const
};

prop_ptr g_prop_tbl[] = {
  _prop0,			// sigma( u / T )
  _prop1,			// sigma( u * T )
  _prop2,
  _prop3,
  _prop4,
};

// annealing_ptr g_annealing_schedule_tbl[]={
//   _annealing_schedule0,		// const value 0.01
//   _annealing_schedule1,		// 1./(time + 1.)
//   _annealing_schedule2,		// __C/(time + __D)
//   _annealing_schedule3,		// _D/log(time)
//   _annealing_schedule4,		// 1/time
// };

// init_layer_ptr g_init_layer_tbl[]={
//   _init_layer0,			// matrix elements is 0 initialize
//   _init_layer1,			// random*scale random domain is [0,1)
//   _init_layer2,			// random*scale random. diagonal is zero 
//   _init_layer3,			// const value init. diagonal is const
//   _init_layer4,			// const value init. diagonal is zero
//   _init_layer5,			// random*scale random domain is[-0.5,0.5)
//   _init_layer6,			// random*scale diag is zero.
//};

// Layer _init_layer0(int D, double scale)
// {
//   Layer L(D,D);
//   UDITER(L){ L(i,j) = L(j,i) = 0.0; }
//   return L;
// }

// Layer _init_layer1(int D, double scale)
// {
//   Layer L(D,D);
//   UDITER(L){L(i,j) = L(j,i) = scale*drand48();}
//   return L;
// }

// Layer _init_layer2(int D, double scale)
// {
//   Layer L = _init_layer2(D,scale);
//   for (size_t i=0; i<D; i++) L(i,i)=0.0;
//   return L;
// }

// Layer _init_layer3(int D, double scale)
// {
//   Layer L(D,D);
//   UDITER(L){L(i,j) = L(j,i) = scale;}
//   return L;
// }

// Layer _init_layer4(int D, double scale)
// {
//   Layer L = _init_layer3(D,scale);
//   for (size_t i=0; i<D; i++) L(i,i) = 0.0;
//   return L;
// }

// Layer _init_layer5(int D, double scale)
// {
//   Layer L(D,D);
//   UDITER(L){ L(i,j) = L(j,i) = scale*(drand48()-0.5);}
//   return L;
// }

// Layer _init_layer6(int D, double scale)
// {
//   Layer L = _init_layer5(D,scale);
//   for (size_t i=0; i<D; i++){ L(i,i) = 0.0; }
//   return L;
//}

double
_energy0(Pattern& v, Layer& L)
{
  double ene=0.0;
  for (size_t i=1; i<L.size1(); i++){
    ene += v[i]*L(0,i);
    for (size_t j=i+1; j<L.size2(); j++)
      ene += v[i]*L(i,j)*v[j];
  }
    
  return -ene;
}

double
_energy1(Pattern& v, Layer& L)
{
  using namespace boost::numeric::ublas;
  double term = inner_prod(prod(v,L),v);
  return -0.5*term;
}

double
_energy2(Pattern& v, Layer& L)
{
  double val2 = 0.0;
  for (size_t i=0; i<L.size1(); i++)
    for (size_t j=0; j<L.size2(); j++)
      val2 += v[i]*L(i,j)*v[j];
  return -0.5*val2;
}

// collect! for Kurata's model. result is "energy_cmp.cc"
double
_energy3(Pattern& v, Layer& L)
{
  double val3 = 0.0;
  for (size_t i=0; i<L.size1(); i++)
    for (size_t j=i; j<L.size2(); j++)
      val3 += v[i]*L(i,j)*v[j];
  return -val3;
}

double
_energy4(Pattern& v, Layer& L)
{
  double val4 = 0.0;
  for (size_t i=0; i<L.size1(); i++)
    for (size_t j=i+1; j<L.size2(); j++)
      val4 += v[i]*L(i,j)*v[j];
  return -val4;
}

double
_energy5(Pattern& v, Layer& L)
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

double
_energy6(Pattern& v, Layer& L)
{
  double val6 = 0.0;
  val6 = inner_prod(prod(v,L),v);
  val6 -= trace(L);
  return -0.5*val6;
}

double
_energy7(Pattern& v, Layer& L)
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
double
_energy8(Pattern& v, Layer& L)
{
  double ret = 0.0;
  for (size_t i=0; i<L.size1(); i++)
    for (size_t j=i; j<L.size2(); j++)
      ret += v[i]*L(i,j)*v[j];
  return -ret;
}

// I think this is hold water.
// Kurata's Report 「ボルツマン・マシン」 describe this
double
_energy9(Pattern& v, Layer& L)
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
double
_energy10(Pattern& v, Layer& L)
{
  double ret=0.0;
  for (size_t i=0; i<L.size1(); i++){
    ret += v[i]*L(i,i);
    for (size_t j=i+1; j<L.size2(); j++)
      ret += v[i]*L(i,j)*v[j];
  }
  return -ret;
}

// double
// energy_test(Pattern& v, Layer& L, std::ostream &os=std::cout)
// {
//   _p(os,_energy0(v,L));
//   _p(os,_energy1(v,L));
//   _p(os,_energy2(v,L));
//   _p(os,_energy3(v,L));
//   _p(os,_energy4(v,L));
//   _p(os,_energy5(v,L));
//   _p(os,_energy6(v,L));
//   _p(os,_energy7(v,L));
//   _p(os,_energy8(v,L));
//   _p(os,_energy9(v,L));
//   _p(os,_energy10(v,L));
//   return energy(v,L);
// }

double
_denergy0(size_t i, Pattern& v, Layer& L)
{
  double u_i=0.0;
  for (size_t j=1; j<v.size(); j++)
    u_i += L(i,j)*v[j];
  return u_i + L(0,i);
}

// collect! result is "energy_cmp.cc"
double
_denergy1(size_t i, Pattern& v, Layer& L)
{
  double u_i=0.0;
  for (size_t j=0; j<v.size(); j++)
    u_i += L(i,j)*v[j];
  return u_i - L(i,i)*v[i] + L(i,i);
}

double
_denergy2(size_t i, Pattern& v, Layer& L)
{
  double u_i=0.0;
  for (size_t j=0; j<v.size(); j++)
    u_i += L(i,j)*v[j];
  return u_i - L(i,i)*v[i];
}

double
_denergy3(size_t i, Pattern& v, Layer& L)
{
  double u_i=0.0;
  for (size_t j=0; j<v.size(); j++)
    u_i += L(i,j)*v[j];
  return u_i;
}

// collect! for Kurata's model. result is "energy_cmp.cc"
double
_denergy4(size_t i, Pattern& v, Layer& L)
{
  double u_i=0.0;
  for (size_t j=0; j<i; j++)
    u_i += L(i,j)*v[j];
  for (size_t j=i+1; j<v.size(); j++)
    u_i += L(i,j)*v[j];
  u_i += L(i,i);// bias term
  return u_i;
}

// double
// denergy_test(size_t i, Pattern& v, Layer& L, std::ostream &os=std::cout)
// {
//   _p(os,_denergy0(i,v,L));
//   _p(os,_denergy1(i,v,L));
//   return denergy(i,v,L);
// }

double
_prop0(int i, Pattern& v, Layer& L, double temperature)
{
  return sigma(denergy(i,v,L)/temperature);
}

double
_prop1(int i, Pattern& v, Layer& L, double temperature)
{
  return sigma(denergy(i,v,L)*temperature);
}

double _prop2(int i, Pattern& v, Layer& L, double temperature)
{
  return 0.0;
}

double _prop3(int i, Pattern& v, Layer& L, double temperature)
{
  return 0.0;
}

double _prop4(int i, Pattern& v, Layer& L, double temperature)
{
  return 0.0;
}

// double _annealing_schedule0(double time)
// {				// const
// #ifdef TEMPERATURE
//   return TEMPERATURE;
// #else
//   return 0.01;
// #endif
// }

// double _annealing_schedule1(double time)
// {
//   return 1.0/(time + 1);
// }

// double _annealing_schedule2(double time)
// {
// #ifndef __C			// Largest Energy Barrier
// #define __C 5.0
// #endif
// #ifndef __D			// usually set equal to one
// #define __D 1.0
// #endif
//   return __C/(time + __D);
// #undef __C
// #undef __D
// }

// double _annealing_schedule3(double time)
// {
// #ifndef _D
// #define _D (0.1)
// #endif
//   return _D/log(time);
// #undef _D
// }

// double _annealing_schedule4(double time)
// {
//   return 1./time;
//}

// double prop_test(int i, Pattern& v, Layer& L,
// 		 double temperature,std::ostream& os=std::cout)
// {
//   _p(os,_prop0(i,v,L,temperature));
//   _p(os,_prop1(i,v,L,temperature));
//   _p(os,_prop2(i,v,L,temperature));
//   _p(os,_prop3(i,v,L,temperature));
//   return prop(i,v,L,temperature);
// }

// Layer
// init_layer(int D, double scale)
// {
//   static init_layer_ptr fn = g_init_layer_tbl[g_init_layer_idx];
//   return fn(D,scale);
// }

double
energy(Pattern& v, Layer& L)
{
  static energy_ptr fn = g_energy_tbl[g_energy_idx];
  return fn(v,L);
}

double
denergy(size_t i, Pattern &v, Layer &L)
{
  static denergy_ptr fn = g_denergy_tbl[g_denergy_idx];
  return fn(i,v,L);
}

double
prop(int i, Pattern &v, Layer &L, double temperature)
{
  static prop_ptr fn = g_prop_tbl[g_prop_idx];
  return fn(i,v,L,temperature);
}

// double
// annealing_schedule(double time)
// {
//   static annealing_ptr fn = g_annealing_schedule_tbl[g_annealing_schedule_idx];
//   return fn(time);
// }

double
partition_function(Layer& L, Patterns& all_v)
{
  double ret=0.0;
  for(Patterns::iterator v = all_v.begin(); v != all_v.end(); ++v){
    ret += exp(-energy(*v,L));
  }
  return ret;
}

double
boltzmann_distribution(Pattern& v,Layer& L, Patterns& all_v)
{
  double numer = exp(-energy(v,L));
  double denom = partition_function(L,all_v);
  return numer/denom;
}

void
kstep_gibbs(size_t K, Pattern& v, Layer& L, double temperature)
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

void
kstep_gibbs_count_onbit(size_t K, Probability &vf, Layer &L,
			double temperature, Probability &gibbs_P)
{
  size_t D = vf.size();
  Pattern v(D,0);
  for(size_t k=0; k < K; k++){
    size_t i;
#ifdef DATEIMP
    i = (rand()%(D-1)) + 1;
#else
    i = rand()%D;
#endif
    v[i] = (prop(i, v, L, temperature) > drand48()) ? ONBIT : OFFBIT;
    gibbs_P[btod(v)]++;
    vf += v;
  }
  vf /= K;
}

void
gibbs_sampling(Pattern& v,  Layer& L, double temperature,
	       Probability &onbit_count)
{
  size_t D = v.size();
  size_t i;
#ifdef DATEIMP
  i = (rand()%(D-1)) + 1;
#else
  i = rand()%D;
#endif
  v[i] = (prop(i, v, L, temperature) > drand48()) ? ONBIT : OFFBIT;
  onbit_count[btod(v)]++;
}
