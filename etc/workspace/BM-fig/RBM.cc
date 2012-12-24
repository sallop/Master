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
#include <boost/tuple/tuple.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#define ONBIT  1
#define OFFBIT 0

namespace ub = boost::numeric::ublas;

typedef double UnitType;
typedef ub::vector<UnitType>      Pattern    ;
typedef std::vector<Pattern>      Patterns   ;
typedef ub::matrix<UnitType>      Layer      ;
typedef ub::zero_vector<UnitType> ZeroPattern;
typedef ub::zero_matrix<UnitType> ZeroLayer  ;
typedef std::vector<double> Threshold;
typedef std::vector<double> Probability;

class TrainData
{
  Patterns    &patterns_ ;
  Probability &prob_     ;
  Probability cumulative_;
public:
  TrainData(Patterns &patterns, Probability &prob) :
    patterns_(patterns), prob_(prob), cumulative_(prob.size(), 0.0)
  {
    int i;
    double sum;
    std::cout << "patterns[0].size() = " << patterns[0].size() << std::endl;
    std::cout << "prob.size()        = " << prob.size()        << std::endl;
    std::cout << "patterns.size()    = " << patterns.size()    << std::endl;
    std::cout << "cumulative_.size() = " << cumulative_.size() << std::endl;
    assert(prob_.size() == cumulative_.size());
    // std::partial_sum(prob_.begin(), prob_.end(), cumulative_.begin());
    // 0.2     , 0.25    , 0.3     , 0.4     , 0.5     , 0.9     , 1
    for(i=0, sum = 0; i < cumulative_.size(); i++){
      sum += prob_[i];
      cumulative_[i] = sum;
    }

    std::cout << "prob_[i]" << "\t" << "cumulative_[i]" << std::endl;
    for(i=0; i < cumulative_.size(); i++){
      std::cout << prob_[i] << "\t" << cumulative_[i] << std::endl;
    }
    //assert(false);
  }

  double get_cumulative(size_t idx){ return cumulative_[idx]; }
  double n_traindata(){ return patterns_.size(); }
  double d_pattern(){ return patterns_[0].size(); }
};

template<typename T> inline
void symmetrization_U(ub::matrix<T>& W)
{
  size_t i, j, size = W.size1();
  assert(W.size1() == W.size2());
  for(i=0; i < size; i++){
    for(j=i; j < size; j++)
      W(j,i)=W(i,j);
    W(i,i)=0.0;
  }
}

// Pattern dtob(int n, size_t size)
// {
//   Pattern ret(size, 0);
//   for (size_t i=0; i < size; i++){ ret[i] = (n>>i) & 0x01;}
//   return ret;
// }
// int btod(const Pattern& x)
// {
//   size_t i, ret=0, size=x.size();
//   for(i=0; i < size; i++){ ret += x[i]*(0x01 << i);}
//   return ret;
//}

Pattern dtob(int n, size_t size)
{
  Pattern ret(size, 0);
  for (size_t i=0; i < size; i++){
    ret[size-i-1] = (n>>i) & 0x01;
  }
  return ret;
}

int btod(const Pattern& x)
{
  size_t i, ret = 0, size = x.size();
  for(i=0; i < size; i++){ ret += x[size-i-1]*(0x01 << i);}
  return ret;
}

template<typename Iterator>
void dtob(int d, Iterator p0, Iterator pN)
{
  int n = pN - p0 + 1;
  while (n--){*pN-- = d & 0x01; d >>= 1;}
}

template<typename Iterator>
int btod(Iterator p0, Iterator pN)
{
  int n = 0, MSB = 0x01 << (pN-p0);
  while (MSB){ n += (*p0++)*MSB; MSB >>= 1;}
  return n;
}

double sigma(double x)
{
  return 1./(1+exp(-x));
}

double energy(Pattern& v, Pattern& h, Layer& W, Threshold& a, Threshold& b)
{
  size_t i, j, k, m, D=v.size(), P=h.size();
  double term1=0.0,term2=0.0,term3=0.0;
  for(i=0;i<D;i++)
    term1 += a[i]*v[i];
  for(i=0;i<P;i++)
    term2 += b[i]*h[i];
  for(i=0;i<D;i++)
    for(j=0;j<P;j++)
      term3 += v(i)*W(i,j)*h(j);
  return - term1 - term2 - term3;
}

inline double mean_field_update(int j, Layer& W, Pattern& v)
{
  size_t i, m, D = W.size1(), P = W.size2();
  double term1 = 0, term2 = 0;
  for(i=0; i<D; i++) term1 += W(i,j)* v(i);
  return sigma(term1);
}

inline double propup  (int j, Layer& W, Pattern& v, Threshold& b)
{
  size_t i, m, D = W.size1(), P=W.size2();
  double term1=0;
  for(i=0; i<D; i++) term1 += W(i,j)*v(i);
  return sigma(b[j] + term1);
}

inline double propdown(int i, Layer& W, Pattern& h, Threshold& a)
{
  size_t j, k, D = W.size1(), P = W.size2();
  double term1=0;
  for(j=0; j<P; j++) term1 += W(i,j)*h(j);
  return sigma(a[i] + term1);  
}

// Z(\theta)
double partition_function(Layer& W, Threshold& a, Threshold& b, Patterns& all_v, Patterns& all_h)
{
  size_t i, j, nV = all_v.size(), nH = all_h.size();
  double ret=0.0;
  assert(W.size1() == all_v[0].size());
  assert(W.size2() == all_h[0].size());
  for(i=0; i < nV; i++)
    for(j=0; j < nH; j++){
      Pattern& v = all_v[i], &h = all_h[j];
      ret += exp(-energy(v,h,W,a,b));
    }
  return ret;
}

//double boltzmann_distribution(Pattern& v,
double prob_model_assign_v(Pattern&      v,
			   Layer&        W,
			   Threshold&    a,
			   Threshold&    b,
			   Patterns& all_v,
			   Patterns& all_h)
{
  size_t i, j, k, nV=all_v.size(), nH=all_h.size();
  double numer=0.0, denom=0.0;
  for(i = 0; i < nH; i++)
    numer += exp(-energy( v, all_h[i], W, a, b));
  denom = partition_function(W, a, b, all_v, all_h);
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
    for(j=0; j < W.size2(); j++)
      printf("%8.6lf\t", W(i,j));
    printf("\n");
  }
}

int main(int argc, char **argv)
{
  size_t d, i, j, k, m, n, p, r, t;
  const int    D = 3, P=2;
  const size_t K = 50, M=1000, N=1000, T=1000;
  //const size_t K=100, M=10, N=10, T=100;
  //const int D=4, P=2; //const int D=4, P=4; //const int D=8, P=4;
  //const size_t K=5, M=10, N=10, T=100; //const size_t K=100, M=10, N=10, T=100;
  Layer     W(D,P), W_data(D,P), W_model(D,P);
  ZeroLayer W_zero (D,P);
  Pattern  v(D,0), vf(D,0), h(P,0), hf(P,0), mu(P,0);
  ZeroPattern v_zero(D), h_zero(D), mu_zero(P);

  Threshold a(D,0), b(P,0);

  Probability prob_P(0x01<<D, 0), prob_Q(0x01<<D, 0), free_running(0x01<<D, 0);
  double alpha = 1.0, kl, rnd;	// learning rate, store the kl divergence
  // all state
  Patterns all_v(0x01<<D, Pattern(D,0));
  Patterns all_h(0x01<<P, Pattern(P,0));
  for(i=0; i < (0x01<<D); i++){ all_v[i] = dtob(i,D);}
  for(i=0; i < (0x01<<P); i++){ all_h[i] = dtob(i,P);}
  FILE* fp;
  // circit's all_state
  prob_Q[0] = 0.10;
  prob_Q[1] = 0.10;
  prob_Q[2] = 0.05;
  prob_Q[3] = 0.05;
  prob_Q[4] = 0.10;
  prob_Q[5] = 0.10;
  prob_Q[6] = 0.40;
  prob_Q[7] = 0.10;

  // setting a train data
  TrainData train_data(all_v, prob_Q);

  //1. Randomly initialize mu and
  for(i=0;i<D;i++){for(j=0;j<P;j++){W(i,j)=0.01*(drand48()-0.5);}}

  fp = fopen("dir-dat/t-kl.dat","w");
  // 1. Randomly initialize parameters \theta^0 and M fantasy particles
  for(t=0; t<T; t++){
    W_data = W_zero;

    //(a) n=1 to N, get a data dependent expection
    for(n=0; n<N; n++){
      rnd = drand48();
      for(k=0; k < train_data.n_traindata(); k++)
	if(rnd < train_data.get_cumulative(k))
	  break;
      // mean-field updates until convergence.
      // this values convergence about 10-20
      Pattern& vn = all_v[k];//&mu = mu_N[k];
      for(i=0; i <  P ; i++){ mu[i] = 0.01*(drand48()-0.5);}
      for(i=0; i < 100; i++){
	j = rand()%P;	
	mu[j] = mean_field_update(j, W, vn);
      }
      W_data += ub::outer_prod(vn, mu);
    }

    //(b) m=1 to M, get a model dependent expection
    W_model = W_zero;
    // randomly initialize mu
    for(m=0; m<M; m++){
      // k-step Gibbs sampler using 4,5
      for(k=0; k<K; k++){
	j = rand()%P;
	i = rand()%D;
	hf[j] = (propup  (j,W,vf,b) > drand48()) ? ONBIT : OFFBIT;
	vf[i] = (propdown(i,W,hf,a) > drand48()) ? ONBIT : OFFBIT;
      }// end k
      W_model += ub::outer_prod(vf,hf);
    }// end m

    //(c) update
    {
      W = W + alpha*(W_data/N - W_model/M);
    }
    alpha = 1.0/(t+1);

    for(int i=0; i < all_v.size(); i++){
      Pattern &v = all_v[i];
      prob_P[i] = prob_model_assign_v(v, W, a, b, all_v, all_h);
    }
    fprintf(fp, "%d %lf\n", t, KLd(prob_P,prob_Q));
  }
  fclose(fp);
  
  fp=fopen("dir-dat/Layer.dat","w");
  fprintf(fp, "W=\n");
  for(i=0; i < D; i++){
    for(j=0; j < P; j++)
      fprintf(fp, "%2.1lf\t", W(i,j));
    fprintf(fp,"\n");
  }
  fprintf(fp,"\n");
  fclose(fp);

  //  k-step Gibbs sampler using 4,5
  for(t=0; t<10000; t++){
    for(k=0; k<50; k++){
      j = rand()%P;
      i = rand()%D;
      h[j] = (propup  (j,W,v,b) > drand48()) ? ONBIT : OFFBIT;
      v[i] = (propdown(i,W,h,a) > drand48()) ? ONBIT : OFFBIT;
    }// end k
    //free_running[btod<Pattern::iterator>(v.begin(), v.end())]++;
    free_running[btod(v)]++;
  }
  double sum = std::accumulate(free_running.begin(),
			       free_running.end()  , 0.0);
  std::cout << "state\t\t"
	    << "Q(x)\t"
	    << "run\t"
	    << "p(v^n)"
	    << std::endl;
  for(i=0; i < free_running.size(); i++){
    std::cout << all_v[i]  << "\t"
   	      << prob_Q[i] << "\t"
   	      << free_running[i]/sum << "\t"
   	      << prob_model_assign_v(all_v[i], W, a, b, all_v, all_h) << "\t"
   	      << std::endl;
  }
  return 0;
}
