#include "BM.hh"

double g_init_alpha;		     // append 2012/01/21
size_t g_annealing_schedule_idx;
size_t g_boltzmann_neural_field_idx; // append 2012/01/21
size_t g_denergy_dw_idx;
size_t g_denergy_idx;
size_t g_denergy_up_idx;
size_t g_energy_idx;
size_t g_init_layer_idx;
size_t g_neighbor_k;		// append 2012/01/21
size_t g_prop_idx;
unsigned g_seed;

// using FBM,GBM,RBM at all BM.
annealing_ptr g_annealing_schedule_tbl[]={
  _annealing_schedule0,		// const value 0.01
  _annealing_schedule1,		// 1./(time + 1.)
  _annealing_schedule2,		// __C/(time + __D)
  _annealing_schedule3,		// _D/log(time)
  _annealing_schedule4,		// 1/time
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

// using FBM energy functions
// energy_ptr g_energy_tbl[]={
//   _energy0 ,// threshold is L(0,i) and L(i,0). date's impliments
//   _energy1 ,// simple inner prod and then 0.5 times
//   _energy2 ,// for loop version inner prod
//   _energy3 ,// correct code at Kurata's matrix. diag is threshold.
//   _energy4 ,// upper triangle elements
//   _energy5 ,// (all elements - diagonal elements)*0.5
//   _energy6 ,// (inner_prod - diagonal)*0.5
//   _energy7 ,// (all elements - diagonal_const)*0.5
//   _energy8 ,// correct code upper triangle matrix. same way energy3
//   _energy9 ,// threshold is const in a Kurata's report [boltzmann machine]. 
//   _energy10,// threshold is one side judge.
// };

// using FBM energy functions
// denergy_ptr g_denergy_tbl[]={
//   _denergy0,// linear sum + threshold_const. date's impliments.
//   _denergy1,// collect linear sum - self input + threshold_const
//   _denergy2,// linear sum - self inpu
//   _denergy3,// linear sum
//   _denergy4,// collect linear sum - self input + threshold_const
// };

// using FBM energy functions
// prop_ptr g_prop_tbl[] = {
//   _prop0,			// sigma( u / T )
//   _prop1,			// sigma( u * T )
//   _prop2,
//   _prop3,
//   _prop4,
// };

// energy_gbm_ptr g_energy_gbm_tbl[]={
//   _energy0,
// };
// denergy_gbm_ptr g_denergy_up_gbm_tbl[]={
//   _denergy_up0,
// };
// denergy_gbm_ptr g_denergy_dw_gbm_tbl[]={
//   _denergy_dw0,
// };

// General purpose Functions
Layer
init_layer(int D, double scale)
{
  static init_layer_ptr fn = g_init_layer_tbl[g_init_layer_idx];
  return fn(D,scale);
}

double
annealing_schedule(double time)
{
  static annealing_ptr f = g_annealing_schedule_tbl[g_annealing_schedule_idx];
  return f(time);
}

Pattern
dtob(int n, size_t size)
{
  Pattern ret(size, 0);
#ifdef DATEIMP
  for (size_t i=0; i < size-1; i++)
    ret[size-1-i] = (n>>i) & 0x01;
#else
  for (size_t i=0; i < size; i++)
    ret[size-1-i] = (n>>i) & 0x01;
#endif
  return ret;
}

int
btod(const Pattern& x)
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

double
sigma(double x)
{
  return 1./(1.+exp(-x));
}

double
KLd(Probability& p, Probability& q)
{
  size_t i, size = (p.size() == q.size()) ? p.size() : 0;
  double d=0.0;
  for (i=0; i < size; i++){
    double term1, term2;
    if (!isnan(term1 = p[i]*log(p[i]) ))
      d += term1;
    if (!isnan(term2 = p[i]*log(q[i]) ))
      d -= term2;
  }
  return d;
}

double
trace(Layer& W)
{
  using namespace boost::numeric::ublas;
  size_t D = W.size1();
  matrix_vector_range<Layer> diag(W,range(0,D),range(0,D));
  return sum(diag);
}

// Sub Functions
// annealing schedule dependent time. this value is const.
double _annealing_schedule0(double immutable_alpha)
{
#ifdef TEMPERATURE
  return TEMPERATURE;
#else
  return immutable_alpha;
#endif
}

double _annealing_schedule1(double time)
{
  return 1.0/(time + 1);
}

double _annealing_schedule2(double time)
{
#ifndef __C			// Largest Energy Barrier
#define __C 10
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
  return _D/log(time+1);
#undef _D
}

double _annealing_schedule4(double time)
{
  return 1./time;
}

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
