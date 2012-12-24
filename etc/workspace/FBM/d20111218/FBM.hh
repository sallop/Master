#ifndef FBM_H_
#define FBM_H_

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
#include <functional>
#include <numeric>
#include <execinfo.h>
#include <exception>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/banded.hpp>

#define ONBIT  1
#define OFFBIT 0

#define SETW setw(16)

#define _eq(val) std::setiosflags(std::ios::left) << #val << "="	\
  << std::setw(16) << std::setiosflags(std::ios::left) << val

#define _p(os,val) os << _eq(val) << std::endl

namespace ub = boost::numeric::ublas;
typedef double UnitType;
typedef ub::vector<UnitType> Pattern;
typedef ub::zero_vector<UnitType> zerov;
typedef ub::zero_matrix<UnitType> zerom;
typedef ub::vector<Pattern>  Patterns;
typedef ub::matrix<UnitType> Layer  ;
typedef ub::matrix_row   <Layer> LayerRow;
typedef ub::matrix_column<Layer> LayerCol;
typedef ub::vector<double> Threshold;
typedef ub::vector<double> Probability;
typedef ub::symmetric_adaptor<Layer,ub::lower> symmetric_view;
typedef ub::banded_adaptor<Layer> banded_view;

typedef double (*energy_ptr)(Pattern& s, Layer& W);
typedef double (*denergy_ptr)(size_t i, Pattern& s, Layer& W);
typedef double (*prop_ptr)(int i, Pattern& s, Layer& W, double temperature);

extern energy_ptr g_energy_tbl[];
extern denergy_ptr g_denergy_tbl[];
extern prop_ptr g_prop_tbl[];
extern size_t g_energy_idx;
extern size_t g_denergy_idx;
extern size_t g_prop_idx;

template<typename T>
void pmat(T& m, std::ostream& os = std::cout)
{
  using namespace std;
  for (size_t i=0; i<m.size1(); i++){
    ub::matrix_row<T> mr(m,i);
    os << mr << "\n";
  }
}

void pvec(Pattern v, std::ostream& os);
inline double trace(Layer& W);
void pprobability(Probability& Pr, std::ostream& os);
Pattern dtob(int n, size_t size);
int btod(const Pattern& x);
double sigma(double x);

double energy_test(Pattern& v, Layer& L,std::ostream &os);
double denergy_test(size_t i, Pattern& v, Layer& L,std::ostream &os);
double prop_test(int i, Pattern& v, Layer& L,double temperature,std::ostream& os);

double _energy0(Pattern& v, Layer& L);
double _energy1(Pattern& v, Layer& L);
double _energy2(Pattern& v, Layer& L);
double _energy3(Pattern& v, Layer& L);
double _energy4(Pattern& v, Layer& L);
double _energy5(Pattern& v, Layer& L);
double _energy6(Pattern& v, Layer& L);
double _energy7(Pattern& v, Layer& L);
double _energy8(Pattern& v, Layer& L); // collect code candidate
double _energy9(Pattern& v, Layer& L); // hold water
double _energy10(Pattern& v, Layer& L);// maybe same as _energy9
double energy(Pattern& v, Layer& L);

double _denergy0(size_t i, Pattern& v, Layer& L);
double _denergy1(size_t i, Pattern& v, Layer& L);
double _denergy2(size_t i, Pattern& v, Layer& L);// correct code candidate
double _denergy3(size_t i, Pattern& v, Layer& L);
double _denergy4(size_t i, Pattern& v, Layer& L);// I think this is hold water
double denergy(size_t i, Pattern& v, Layer& L);

double _prop0(int i, Pattern& v, Layer& L, double temperature);
double _prop1(int i, Pattern& v, Layer& L, double temperature);
double _prop2(int i, Pattern& v, Layer& L, double temperature);
double _prop3(int i, Pattern& v, Layer& L, double temperature);
double _prop4(int i, Pattern& v, Layer& L, double temperature);
double prop(int i, Pattern& v, Layer& L, double temperature);

double partition_function(Layer& L, Patterns& all_v);
double boltzmann_distribution(Pattern& v,Layer& L, Patterns& all_v);
double KLd(Probability& p, Probability& q);
void kstep_gibbs(size_t K, Pattern& v, Layer& L,double temperature);

void test(int D,int P,int K,int M,int N,int T,char prefix[32],
	  Layer &L,Layer &data,Layer &model,
	  symmetric_view &sL,
	  symmetric_view &sdata,
	  symmetric_view &smodel,
	  banded_view &bL,
	  banded_view &bdata,
	  banded_view &bmomdel,
	  Pattern &v,Pattern &vf,
	  Probability &boltz_P,
	  Probability &gibbs_P,
	  Probability &prob_Q,
	  Probability &accm_Q,
	  double alpha,double kl,Patterns& all_v,
	  std::ostream& os);

#endif	// FBM_H_
