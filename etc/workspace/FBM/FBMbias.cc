#include <iostream>
#include <cmath>
#include <cmath>
#include <cfloat>
#include <cstdlib>
#include <cassert>
#include <time.h>
#include <iostream>
#include <iomanip>
#include <fstream>
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

typedef double UnitType;
typedef boost::numeric::ublas::vector<UnitType> Pattern;
typedef boost::numeric::ublas::zero_vector<UnitType> zerov;
typedef boost::numeric::ublas::zero_matrix<UnitType> zerom;
typedef boost::numeric::ublas::vector<Pattern>  Patterns;
typedef boost::numeric::ublas::matrix<UnitType> Layer  ;
typedef boost::numeric::ublas::matrix_row   <Layer> LayerRow;
typedef boost::numeric::ublas::matrix_column<Layer> LayerCol;
typedef boost::numeric::ublas::vector<double> Threshold;
typedef boost::numeric::ublas::vector<double> Probability;

namespace ub = boost::numeric::ublas;

Pattern dtob(int n, size_t size)
{
  Pattern ret(size, 0);
  for (size_t i=0; i < size; i++){
    ret[size-1-i] = (n>>i) & 0x01;
  }
  return ret;
}

int btod(const Pattern& x)
{
  size_t i, ret=0, size=x.size();
  for(i=0; i < size; i++){
    ret += x[size-1-i]*(0x01 << i);
  }
  return ret;
}

double sigma(double x)
{
  return 1./(1+exp(-x));
}

double energy(Pattern& v, Layer& L)
{
  using namespace boost::numeric::ublas;
  double term = inner_prod(prod(v,L),v);
  return -0.5*term;
}


inline double prop(int i, Pattern& v, Layer& L, double temperature)
{
  //ub::matrix_row<Layer> Li(L,i);
  LayerRow Li(L,i);
  double term = ub::inner_prod(v,Li);
  return sigma(term/temperature);
  //  return sigma(term - Li[i]*v[i]);
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
  for (i=0; i < size; i++){d += p[i]*log(p[i]/q[i]);}
  return d;
}

void print(Layer& W)
{
  size_t i, j;
  for(i=0; i < W.size1(); i++){
    ub::matrix_row<Layer> wr(W,i);
    std::cout << wr << std::endl;
  }
}

inline void kstep_gibbs(size_t K, Pattern& v, Layer& L, double temperature)
{
  size_t D = v.size();
  for(size_t k=0; k < K; k++){
    size_t i = rand()%D;
    v[i] = (prop(i, v, L, temperature) > drand48()) ? ONBIT : OFFBIT;
  }
}

template<typename T>
void pmat(T& m)
{
  for (size_t i=0; i < m.size1(); i++){
    ub::matrix_row<T> mr(m,i);
    std::cout << mr << std::endl;
  }
  std::cout << std::endl;
}

int main(int argc, char **argv)
{
  //  const int D=3, P=2;
  //  const int D=3, P=2;
  const int D=3, P=3;
  const size_t K=5, M=1000, N=1000, T=10000;
  //  const size_t K=50, M=2000, N=2000, T=1000;
  //const size_t K=50, M=3000, N=3000, T=1000;
  //  const size_t K=50, M=3000, N=1000, T=1000;
  //  const size_t K=50, M=4000, N=4000, T=5000;
  //  const size_t K=100, M=4000, N=4000, T=50000;
  //const size_t K=10, M=4000, N=4000, T=50000;
  Layer L(D,D), data (D,D), model(D,D); // symmetrize
  ub::symmetric_adaptor<Layer,ub::lower> sL(L);
  ub::symmetric_adaptor<Layer,ub::lower> sdata(data);
  ub::symmetric_adaptor<Layer,ub::lower> smodel(model);
  ub::banded_adaptor<Layer> bL(L,0);
  ub::banded_adaptor<Layer> bdata(data,0);
  ub::banded_adaptor<Layer> bmomdel(model,0);

  Pattern v(D,0), vf(D,0);
  Probability boltz_P   (0x01<<D, 0);
  Probability gibbs_P   (0x01<<D, 0);
  Probability gibbs_temp(0x01<<D, 0);
  Probability prob_Q    (0x01<<D, 0);
  Probability accm_Q    (0x01<<D, 0);
  double alpha = 1.0;
  double kl;

  prob_Q[0]=0.10;
  prob_Q[1]=0.10;
  prob_Q[2]=0.05;
  prob_Q[3]=0.05;
  prob_Q[4]=0.10;
  prob_Q[5]=0.10;
  prob_Q[6]=0.40;
  prob_Q[7]=0.10;

  std::partial_sum(prob_Q.begin(),
		   prob_Q.end()  ,
		   accm_Q.begin());
  accm_Q[accm_Q.size()-1] = 1.0;
  // all state
  Patterns all_v(0x01<<D, Pattern(D));
  for(size_t n=0; n < 0x01<<D; ++n){
    all_v[n] = dtob(n,D);
  }
  
  //1. Randomly initialize mu and
  for(size_t i=0;i<D;i++)
    for(size_t j=0;j<D;j++)
      L(i,j)=0.01*drand48();

  std::cout << "L" << std::endl;
  pmat<Layer>(L);
  std::cout << "sL" << std::endl;
  pmat<ub::symmetric_adaptor<Layer,ub::lower> >(sL);
  std::cout << "bL" << std::endl;
  pmat<ub::banded_adaptor<Layer> >(bL);

  
  std::ofstream ofs_L("L3x3.dat");
  std::ofstream ofs_t_gibbs("t_gibbs3x3.dat");
  for(size_t t=1; t<T; t++){
    data  = zerom(D,D);
    model = zerom(D,D);
    // data dependent expection
    for(size_t n=0; n<N; n++){
      size_t k;
      double r = drand48();
      for (k=0; k < accm_Q.size(); k++){
	if (r < accm_Q[k])
	  break;
      }
      data += ub::outer_prod(all_v[k],all_v[k]);
    }
    data /= N;

    // model dependent expection
    gibbs_temp = zerov(0x01<<D);
    for(size_t m=0; m<M; m++){
      vf = zerov(D);
      //kstep_gibbs(K,vf,L,t);
      kstep_gibbs(K,vf,L,alpha);
      model += outer_prod(vf,vf);
      int d = btod(vf);
      gibbs_temp[d]++;
    }
    model /= M;
    gibbs_temp /= M;

    L = L + alpha*(data - model);
    L = sL - bL;
    alpha = 1.0/(t+1);

    ofs_L       << t << "\t" << L          << std::endl;
    ofs_t_gibbs << t << "\t" << gibbs_temp << std::endl;
    // for (size_t i=0; i < 0x01<<D; i++){
    //       std::cout << prob_Q[i]     << "\t"
    // 		<< gibbs_temp[i] << "\t"
    // 		<< boltzmann_distribution(all_v[i],L,all_v)  << "\t"
    // 		<< std::endl;
    //     }
    //    std::cout << std::endl;
    //    pmat<Layer>(L);
    //    assert(false);
  }

  ofs_L.close();
  ofs_t_gibbs.close();




  return 0;
}
