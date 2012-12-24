#include <cstdio>
#include <unistd.h>
#include <time.h>
#include <iostream>
#include <cmath>
#include <cfloat>
#include <cstdlib>
#include <cassert>
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

typedef double				UnitType   ;
typedef ub::vector<UnitType>		Pattern    ;
typedef std::vector<Pattern>		Patterns   ;
typedef ub::matrix<UnitType>		Layer      ;
typedef ub::zero_vector<UnitType>	ZeroPattern;
typedef ub::zero_matrix<UnitType>	ZeroLayer  ;
typedef Pattern				Threshold  ;
typedef Pattern				Probability;

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

Pattern dtob(int n, size_t size)
{
  Pattern ret(size, 0);
  for (size_t i=0; i < size; i++)
    ret[size-i-1] = (n>>i) & 0x01;
  return ret;
}

int btod(const Pattern& x)
{
  size_t i, ret = 0, size = x.size();
  for(i=0; i < size; i++) 
    ret += x[size-i-1]*(0x01 << i);
  return ret;
}

template<typename Iterator>
void dtob(int d, Iterator p0, Iterator pN)
{
  int n = pN - p0 + 1;
  while (n--){
    *pN-- = d & 0x01;
    d >>= 1;
  }
}

template<typename Iterator>
int btod(Iterator p0, Iterator pN)
{
  int n = 0, MSB = 0x01 << (pN-p0);
  while (MSB){
    n += (*p0++)*MSB;
    MSB >>= 1;
  }
  return n;
}

double KLd(Probability& p, Probability& q)
{
  size_t i, size = (p.size() == q.size()) ? p.size() : 0;
  double d=0.0;
  for (i=0; i < size; i++){
    d += p[i]*log(p[i]/q[i]);
  }
  return d;
}

class BM
{
protected:
  std::string name;
  size_t	D_;		// visible dimension
  size_t	P_;		// hidden vector dimension
  size_t	K_;		// k-step gibbs sampling
  size_t	M_;		// Fantasy particle
  size_t	N_;		// number of trainig data
  size_t	T_;		// training time
  double	alpha_;		// Learning rate
  Probability   prob_P_;
  Probability   prob_Q_;	// Environment
  Probability   cumulative_Q_;	// Environment
public:
  double Sigma(double x){ return 1./(1.+exp(-x));}
  virtual double PartitionFunction()=0;
  virtual double ProbModelAssignV(Pattern& v)=0;
  virtual void Learn() = 0;
  BM(size_t D, size_t P, size_t K, size_t M, size_t N, size_t T, size_t alpha,
     Probability environment)
    : prob_P_(0x01<<D,0), prob_Q_(environment), cumulative_Q_(0x01<<D, 0)
  {
    size_t i,j,k; 
    name = "BM";
    D_ = D;
    P_ = P;
    K_ = K;
    M_ = M;
    N_ = N;
    T_ = T;
    alpha_ = alpha;

    cumulative_Q_[0] = prob_Q_[0];
    for (i=1; i<prob_Q_.size(); i++)
      cumulative_Q_[i] = cumulative_Q_[i-1] + prob_Q_[i];
  }

  std::string ExperimentCondition()
  {
    std::stringstream ss;
    ss << "# ";
    ss << "D=" << D_ << ",";
    ss << "P=" << P_ << ",";
    ss << "K=" << K_ << ",";
    ss << "M=" << M_ << ",";
    ss << "N=" << N_ << ",";
    ss << "T=" << T_ << ",";
    ss << "alpha=" << alpha_;
    return ss.str();
  }
  
  virtual void PritProbability()
  {
    Pattern v(D_,0);
    for (size_t i=0; i < 0x01<<D_; i++){
      v = dtob(i,D_);
      printf("P[%d]=%lf\t", i, ProbModelAssignV(v));
    }
    printf("\n");
  }
};

//完全結合のボルツマンマシン
class FBM : public BM
{
protected:
  Layer		L_;
  Threshold     a_;		// bias
  Patterns &all_v_;
public:
  std::string ExperimentCondition()
  {
    std::stringstream ss;
    ss << BM::ExperimentCondition();
    ss << ",";
    ss << "FBM";
    return ss.str();
  }

  FBM(size_t D, size_t K, size_t M, size_t N, size_t T, size_t alpha,
      Probability environment, Patterns &all_v)
    : BM(D,0,K,M,N,T,alpha,environment),
      L_(D,D), a_(D,0), all_v_(all_v)
  {
    size_t i,j,k;
    name = "FBM";
    for (i=0; i<D_; i++){
      a_[i] = 0.01*(drand48()-0.5);
      for (j=0; j<D_; j++){
	L_(i,j) = 0.01*(drand48()-0.5);
      }
    }
    symmetrization_U<double>(L_);
  }

  virtual double Energy(Pattern& v) = 0;
  virtual double LinearSum(int i, Pattern& v) = 0;

  double PartitionFunction()
  {
    size_t i,j,k;
    double z=0.0;
    for (i=0; i<all_v_.size(); i++)
      z += exp(-Energy( all_v_[i] ));
    return z;
  }

  double ProbModelAssignV(Pattern& v)
  {
    size_t i,j,k;
    double numer=0.0, denom=0.0;
    //for (i=0; i<all_v_.size(); i++)
    numer = exp(-Energy(v));
    denom = PartitionFunction();
    return numer/denom;
  }
  
  virtual void Learn()
  {
    size_t i,j,k,m,n,t;
    double rnd;
    Layer     L_data(D_,D_), L_model(D_,D_);
    ZeroLayer L_zero(D_,D_);
    Pattern     a_data(D_), a_model(D_);
    ZeroPattern a_zero(D_), v_zero (D_);
    Pattern vf(D_,0);		// fantasy particle

    FILE *fp;
    double u_i = 0.0;
    // 1. Randomly initialize parameters \theta^0 and M fantasy Particle
    std::stringstream ss;
    ss << "dir-" << name << "/t-kl.dat";
    fp = fopen(ss.str().c_str(),"w");
    fprintf(fp, "%s\n", ExperimentCondition().c_str());
    for (t=0; t<T_; t++){
      // (a) n=1 to N, get a data dependent expection
      L_data = L_zero;
      for(n=0; n<N_; n++){
	rnd = drand48();
	for(i=0; i<cumulative_Q_.size(); i++)
	  if (rnd < cumulative_Q_[i])
	    break;
	Pattern& vi = all_v_[i];
	L_data += ub::outer_prod(vi,vi);
      }//end n

      // (b) m=1 to M, get a model dependent expection
      L_model = L_zero;
      for (m=0; m<M_; m++){
	vf = v_zero;
	// k-step gibbs sampling
	for(k=0; k<K_; k++){
	  i = rand()%D_;
	  //for(j=0; j<D_; j++) u_i += L_(i,j)*vf(j) - a_[i]*vf(j);
	  u_i = LinearSum(i, vf);
	  vf(i) = (Sigma(u_i) > drand48()) ? ONBIT : OFFBIT;
	}
	L_model += ub::outer_prod(vf,vf);
      }//end m

      // (c) update L = L + \Delta L
      //std::cout << "L_data =" << L_data << std::endl;
      //std::cout << "L_model=" << L_model << std::endl;
      //      L_ = L_ + alpha_*(L_data/N_ - L_model/M_);
      L_ = L_ + 0.01*(L_data/N_ - L_model/M_);
      symmetrization_U<double>(L_);
      alpha_ = 1.0/(t+1);
      // print KL-divergence
      for (i=0; i<all_v_.size(); i++)
	prob_P_[i] = ProbModelAssignV(all_v_[i]);
      fprintf(fp, "%d %lf\n", t, KLd(prob_P_, prob_Q_));
    }

    for (i=0; i < all_v_.size(); i++)
      fprintf(fp, "# P(x^%d) = %lf\n", i, prop_P_);
    fclose(fp);
  }
};

class FBMNoBias : public FBM
{
public:
  FBMNoBias(size_t D, size_t K, size_t M, size_t N, size_t T, size_t alpha,
	    Probability environment, Patterns &all_v) 
    : FBM(D, K, M, N, T, alpha, environment, all_v)
  {
    name = "FBMNoBias";
  }

  double Energy(Pattern &v)
  {
    double u = 0.0;
    size_t i, k;
    for (i=0; i < D_; i++)
      for (k=0; k < D_; k++)
	u += v(i)*L_(i,k)*v(k);
    return -0.5*u;
  }
  double LinearSum(int i, Pattern& vf)
  {
    double u_i=0.0;
    for(int k=0; k<D_; k++)
      u_i += L_(i,k)*vf(k);
    return u_i;
  }
};

class FBMBias : public FBM
{
public:
  FBMBias(size_t D, size_t K, size_t M, size_t N, size_t T, size_t alpha,
	  Probability environment, Patterns &all_v) 
    : FBM(D, K, M, N, T, alpha, environment, all_v)
  {
    name = "FBMBias";
  }
  virtual double Energy(Pattern &v);
  virtual double LinearSum(int i, Pattern& vf);
};

class FBMBiasModify : public FBM
{
public:
  FBMBiasModify(size_t D, size_t K, size_t M, size_t N, size_t T, size_t alpha,
		Probability environment, Patterns &all_v) 
    : FBM(D, K, M, N, T, alpha, environment, all_v)
  {
    name = "FBMBiasModify";
  }

  virtual double Energy(Pattern &v)=0;
  virtual double LinearSum(int i, Pattern& vf)=0;
  
  void Learn()
  {
    size_t i,j,k,m,n,t;
    double rnd;
    Layer     L_data(D_,D_), L_model(D_,D_);
    ZeroLayer L_zero(D_,D_);
    Pattern     a_data(D_), a_model(D_);
    ZeroPattern a_zero(D_), v_zero (D_);
    Pattern vf(D_,0);		// fantasy particle

    FILE *fp;
    double u_i = 0.0;
    // 1. Randomly initialize parameters \theta^0 and M fantasy Particle
    std::stringstream ss;
    ss << "dir-" << name << "/t-kl.dat";
    fp = fopen(ss.str().c_str(),"w");
    fprintf(fp, "%s\n", ExperimentCondition().c_str());
    for (t=0; t<T_; t++){
      // (a) n=1 to N, get a data dependent expection
      L_data = L_zero;
      a_data = a_zero;
      for(n=0; n<N_; n++){
	rnd = drand48();
	for(i=0; i<cumulative_Q_.size(); i++)
	  if (rnd < cumulative_Q_[i])
	    break;
	Pattern& vi = all_v_[i];
	L_data += ub::outer_prod(vi,vi);
	a_data += vi;
      }//end n

      // (b) m=1 to M, get a model dependent expection
      L_model = L_zero;
      a_model = a_zero;
      for (m=0; m<M_; m++){
	vf = v_zero;
	// k-step gibbs sampling
	for(k=0; k<K_; k++){
	  i = rand()%D_;
	  //for(j=0; j<D_; j++) u_i += L_(i,j)*vf(j) - a_[i]*vf(j);
	  u_i = LinearSum(i, vf);
	  vf(i) = (Sigma(u_i) > drand48()) ? ONBIT : OFFBIT;
	}
	L_model += ub::outer_prod(vf,vf);
	a_model += vf;
      }//end m

      // (c) update L = L + \Delta L
      L_ = L_ + alpha_*(L_data/N_ - L_model/M_);
      a_ = a_ + alpha_*(a_data/N_ - a_model/M_);
      symmetrization_U<double>(L_);
      alpha_ = 1.0/(t+1);

      //      std::cout << L_ << std::endl;
      // print KL-divergence
      for (i=0; i<all_v_.size(); i++)
	prob_P_[i] = ProbModelAssignV(all_v_[i]);
      fprintf(fp, "%d %lf\n", t, KLd(prob_P_, prob_Q_));
    }
    fclose(fp);
  }
};

class FBMBiasModify1 : public FBMBiasModify
{
public:
  FBMBiasModify1(size_t D, size_t K, size_t M, size_t N, size_t T, size_t alpha,
		 Probability environment, Patterns &all_v) 
    : FBMBiasModify(D, K, M, N, T, alpha, environment, all_v)
  {
    name = "FBMBiasModify1";
  }
  double Energy(Pattern &v);
  double LinearSum(int i, Pattern& vf);
};

class FBMBiasModify2 : public FBMBiasModify
{
public:
  FBMBiasModify2(size_t D, size_t K, size_t M, size_t N, size_t T, size_t alpha,
		 Probability environment, Patterns &all_v) 
    : FBMBiasModify(D, K, M, N, T, alpha, environment, all_v)
  {
    name = "FBMBiasModify2";
  }
  double Energy(Pattern &v);
  double LinearSum(int i, Pattern& vf);
};
class FBMBiasModify3 : public FBMBiasModify
{
public:
  FBMBiasModify3(size_t D, size_t K, size_t M, size_t N, size_t T, size_t alpha,
		 Probability environment, Patterns &all_v) 
    : FBMBiasModify(D, K, M, N, T, alpha, environment, all_v)
  {
    name = "FBMBiasModify3";
  }
  double Energy(Pattern &v);
  double LinearSum(int i, Pattern& vf);
};

class FBMBiasModify4 : public FBMBiasModify
{
public:
  FBMBiasModify4(size_t D, size_t K, size_t M, size_t N, size_t T, size_t alpha,
		 Probability environment, Patterns &all_v) 
    : FBMBiasModify(D, K, M, N, T, alpha, environment, all_v)
  {
    name = "FBMBiasModify4";
  }
  double Energy(Pattern &v);
  double LinearSum(int i, Pattern& vf);
};

class FBMBiasModify5 : public FBMBiasModify
{
public:
  FBMBiasModify5(size_t D, size_t K, size_t M, size_t N, size_t T, size_t alpha,
		 Probability environment, Patterns &all_v) 
    : FBMBiasModify(D, K, M, N, T, alpha, environment, all_v)
  {
    name = "FBMBiasModify5";
  }
  double Energy(Pattern &v);
  double LinearSum(int i, Pattern& vf);
};

double FBMBias::Energy(Pattern &v)
{
  double u = 0.0, bias = 0.0;
  size_t i, j;
  for (i=0; i < D_; i++){
    bias += a_[i]*v(i);
    for (j=0; j < D_; j++){
      u += L_(i,j)*v(i)*v(j);
    }
  }
  return -0.5*u + bias;
}
double FBMBias::LinearSum(int i, Pattern& vf)
{
  double u_i=0.0;
  for(int j=0; j<D_; j++)
    u_i += L_(i,j)*vf(j) - a_[i]*vf(j);    
  return u_i;
}

double FBMBiasModify1::Energy(Pattern &v)
{
  double u = 0.0, bias = 0.0;
  size_t i, j;
  for (i=0; i < D_; i++)
    bias += a_(i)*v(i);
  for (i=0; i < D_; i++)
    for (j=0; j < D_; j++)
      u += L_(i,j)*v(i)*v(j);
  return -0.5*u + bias;
}

double FBMBiasModify1::LinearSum(int i, Pattern& vf)
{
  double u_i=0.0;
  for(int j=0; j<D_; j++)
    u_i += L_(i,j)*vf(j) - a_[i]*vf(j);    
  return u_i;
}

double FBMBiasModify2::Energy(Pattern &v)
{
  double u = 0.0, bias = 0.0;
  size_t i, j;
  for (i=0; i < D_; i++)
    bias += a_(i)*v(i);
  for (i=0; i < D_; i++)
    for (j=0; j < D_; j++)
      u += L_(i,j)*v(i)*v(j);
  return -0.5*u + bias;
}

double FBMBiasModify2::LinearSum(int i, Pattern& vf)
{
  double u_i=0.0;
  for(int j=0; j<D_; j++)
    u_i += L_(i,j)*vf(j) - a_[j]*vf(j);
  return u_i;
}

double FBMBiasModify3::Energy(Pattern &v)
{
  double u = 0.0, bias = 0.0;
  size_t i, j;
  for (i=0; i < D_; i++)
    bias += a_(i)*v(i);
  for (i=0; i < D_; i++)
    for (j=0; j < D_; j++)
      u += L_(i,j)*v(i)*v(j);
  return -0.5*u + bias;
}

double FBMBiasModify3::LinearSum(int i, Pattern& vf)
{
  double u_i=0.0;
  for(int j=0; j<D_; j++)
    u_i += L_(i,j)*vf(j) - a_[i]*vf(i);
  return u_i;
}

double FBMBiasModify4::Energy(Pattern &v)
{
  double u = 0.0, bias = 0.0;
  size_t i, j;
  for (i=0; i < D_; i++)
    bias += a_(i)*v(i);
  for (i=0; i < D_; i++)
    for (j=0; j < D_; j++)
      u += L_(i,j)*v(i)*v(j);
  return -0.5*u + bias;
}

double FBMBiasModify4::LinearSum(int i, Pattern& vf)
{
  double u_i=0.0;
  for(int j=0; j<D_; j++)
    u_i += L_(i,j)*vf(j);
  return a_[i] + u_i;
}

double FBMBiasModify5::Energy(Pattern &v)
{
  double u = 0.0, bias = 0.0;
  size_t i, j;
  for (i=0; i < D_; i++)
    bias += a_(i)*v(i);
  for (i=0; i < D_; i++)
    for (j=0; j < D_; j++)
      u += L_(i,j)*v(i)*v(j);
  return -0.5*u + bias;
}

double FBMBiasModify5::LinearSum(int i, Pattern& vf)
{
  double u_i=0.0;
  for(int j=0; j<D_; j++)
    u_i += L_(i,j)*vf(j);
  return - a_[i] + u_i;
}

//可視層、隱れ層を持つ全結合のボルツマンマシン
//エネルギー関数のバイアス項 \sum_i^D a_i v_i の符号を - に修正
class FBMBiasModify1b : public FBMBiasModify1
{
public:
  FBMBiasModify1b(size_t D, size_t K, size_t M, size_t N, size_t T, size_t alpha,
		  Probability environment, Patterns &all_v)
    : FBMBiasModify1(D, K, M, N, T, alpha, environment, all_v)
  {
    name = "FBMBiasModify1b";
  }
  double Energy(Pattern &v)
  {
    double u = 0.0, bias = 0.0;
    size_t i, j;
    for (i=0; i < D_; i++)
      bias += a_(i)*v(i);
    for (i=0; i < D_; i++)
      for (j=0; j < D_; j++)
	u += L_(i,j)*v(i)*v(j);
    return -0.5*u - bias;
  }
};

class FBMBiasModify2b : public FBMBiasModify2
{
public:
  FBMBiasModify2b(size_t D, size_t K, size_t M, size_t N, size_t T, size_t alpha,
		  Probability environment, Patterns &all_v)
    : FBMBiasModify2(D, K, M, N, T, alpha, environment, all_v)
  {
    name = "FBMBiasModify2b";
  }

  double Energy(Pattern &v)
  {
    double u = 0.0, bias = 0.0;
    size_t i, j;
    for (i=0; i < D_; i++)
      bias += a_(i)*v(i);
    for (i=0; i < D_; i++)
      for (j=0; j < D_; j++)
	u += L_(i,j)*v(i)*v(j);
    return -0.5*u - bias;
  }
};

class FBMBiasModify3b : public FBMBiasModify3
{
public:
  FBMBiasModify3b(size_t D, size_t K, size_t M, size_t N, size_t T, size_t alpha,
		  Probability environment, Patterns &all_v)
    : FBMBiasModify3(D, K, M, N, T, alpha, environment, all_v)
  {
    name = "FBMBiasModify3b";
  }

  double Energy(Pattern &v)
  {
    double u = 0.0, bias = 0.0;
    size_t i, j;
    for (i=0; i < D_; i++)
      bias += a_(i)*v(i);
    for (i=0; i < D_; i++)
      for (j=0; j < D_; j++)
	u += L_(i,j)*v(i)*v(j);
    return -0.5*u - bias;
  }
};

class FBMBiasModify4b : public FBMBiasModify4
{
public:
  FBMBiasModify4b(size_t D, size_t K, size_t M, size_t N, size_t T, size_t alpha,
		  Probability environment, Patterns &all_v)
    : FBMBiasModify4(D, K, M, N, T, alpha, environment, all_v)
  {
    name = "FBMBiasModify4b";
  }

  double Energy(Pattern &v)
  {
    double u = 0.0, bias = 0.0;
    size_t i, j;
    for (i=0; i < D_; i++)
      bias += a_(i)*v(i);
    for (i=0; i < D_; i++)
      for (j=0; j < D_; j++)
	u += L_(i,j)*v(i)*v(j);
    return -0.5*u - bias;
  }
};


int main(int argc, char **argv)
{
  size_t i,j,k,m,n;
  int result;
  size_t D =  3, P = 1;
  size_t K = 5, M = 1000, N=1000, T=1000;
  double alpha = 1.0;
  std::vector<Pattern> all_v(0x01<<D, Pattern(D,0));
  std::vector<Pattern> all_h(0x01<<P, Pattern(P,0));
  for(i=0; i < (0x01<<D); i++)
    all_v[i] = dtob(i,D);
  for(i=0; i < (0x01<<P); i++)
    all_h[i] = dtob(i,P);

  // circit's all_state
  Probability prob_Q(0x01<<D);
  prob_Q[0] = 0.10;
  prob_Q[1] = 0.10;
  prob_Q[2] = 0.05;
  prob_Q[3] = 0.05;
  prob_Q[4] = 0.10;
  prob_Q[5] = 0.10;
  prob_Q[6] = 0.40;
  prob_Q[7] = 0.10;

  while ((result = getopt(argc, argv, "a:d:p:k:m:n:t:f:")) != -1){
    switch(result){
    case 'a':
      alpha = atof( optarg );
      std::cout << "learning rate alpha =" << alpha << std::endl;
      break;
    case 'd':
      //      D = atoi( optarg );
      std::cout << "visible layer degree D is const." << std::endl;
      break;
    case 'p':
      P = atoi( optarg );
      std::cout << "P =" << P << std::endl;
      break;
    case 'k':
      K = atoi( optarg );
      std::cout << "K =" << K << std::endl;
      break;
    case 'm':
      M = atoi( optarg );
      std::cout << "M =" << M << std::endl;
      break;
    case 'n':
      N = atoi( optarg );
      std::cout << "N =" << N << std::endl;
      break;
    case 't':
      T = atoi( optarg );
      std::cout << "T =" << T << std::endl;
      break;
    case 'f':
      //F = atoi( optarg );
      break;
    case '?':
      printf("unknown\n");
    }
  }

  // (x0, x1, x2)の自己相関行列
  BM *fbm_bias   = new FBMBias  (D,K, M, N, T, alpha, prob_Q, all_v);
  BM *fbm_nobias = new FBMNoBias(D,K, M, N, T, alpha, prob_Q, all_v);
  //BM	*fbm1  = new FBMBiasModify1 (D,K, M, N, T, alpha, prob_Q, all_v);
  //BM	*fbm2  = new FBMBiasModify2 (D,K, M, N, T, alpha, prob_Q, all_v);
  //BM	*fbm3  = new FBMBiasModify3 (D,K, M, N, T, alpha, prob_Q, all_v);

  // (1, x0, x1, x2)の自己相関行列 => バイアス項を修正
  BM	*fbm4  = new FBMBiasModify4 (D,K, M, N, T, alpha, prob_Q, all_v);
  //BM	*fbm5  = new FBMBiasModify5 (D,K, M, N, T, alpha, prob_Q, all_v);
  //BM	*fbm1b = new FBMBiasModify1b(D,K, M, N, T, alpha, prob_Q, all_v);
  //BM	*fbm2b = new FBMBiasModify2b(D,K, M, N, T, alpha, prob_Q, all_v);
  //BM	*fbm3b = new FBMBiasModify3b(D,K, M, N, T, alpha, prob_Q, all_v);
  //BM	*fbm4b = new FBMBiasModify4b(D,K, M, N, T, alpha, prob_Q, all_v);

  //  BM *rbm = new RBM(D, P, K, M, N, T, alpha, prob_Q, all_v, all_h);
  //  BM *rbm = new RBM(D, P, K, M, N, T, alpha, prob_Q, all_v, all_h);
  //  BM *gbm = new GBM(D, P, K, M, N, T, alpha, prob_Q, all_v, all_h);

  fbm_bias->Learn();
  fbm_nobias->Learn();
  
  //fbm1->Learn();
  //fbm2->Learn();
  //fbm3->Learn();
  fbm4->Learn();
  //fbm5->Learn();

  // fbm1b->Learn();
  // fbm2b->Learn();
  // fbm3b->Learn();
  // fbm4b->Learn();
  // rbm->Learn();
  // gbm->Learn();

  fbm_bias  ->PritProbability();
  fbm_nobias->PritProbability();
  fbm4      ->PritProbability();
  

  return 0;
}

