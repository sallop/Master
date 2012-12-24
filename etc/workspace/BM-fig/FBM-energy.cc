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

void printLayer(FILE *fp, Layer& W)
{
	int i,j;
	for(i=0; i<W.size1(); i++){
		for(j=0; j<W.size1(); j++){
			fprintf(fp, "%8.4lf\t", W(i,j));
		}
		fprintf(fp, "\n");
	}
}

template<typename T> inline
void symmetrization_U(ub::matrix<T>& W)
{
	size_t i, j, size = W.size1();
	assert(W.size1() == W.size2());
	for(i=0; i < size; i++){
		for(j=i; j < size; j++){
			W(j,i)=W(i,j);
		}
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
      a_[i] = drand48()-0.5;
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
    //std::cout << numer << "/" << denom << std::endl;
    return numer/denom;
  }
  
  void Learn()
  {
    size_t i,j,k,m,n,t;
    double rnd;
    Layer L_data(D_,D_), L_model(D_,D_);
    Pattern vf(D_,0);		// fantasy particle
    ZeroLayer L_zero(D_,D_);
    ZeroPattern v_zero(D_);
    FILE *fp, *fp_data, *fp_model, *fp_L, *fp_prob;
    double u_i = 0.0;
    // 1. Randomly initialize parameters \theta^0 and M fantasy Particle
    std::stringstream ss;
    {
	    ss.str("");
	    ss << "dir-" << name << "/t-kl.dat";
	    fp = fopen(ss.str().c_str(),"w");
	    ss.str("");
	    ss << "dir-" << name << "/DATA.dat";
	    fp_data = fopen(ss.str().c_str(), "w");
	    ss.str("");
	    ss << "dir-" << name << "/MODEL.dat";
	    fp_model = fopen(ss.str().c_str(), "w");
	    ss.str("");
	    ss << "dir-" << name << "/L.dat";
	    fp_L = fopen(ss.str().c_str(), "w");
	    ss.str("");
	    ss << "dir-" << name << "/prob.dat";
	    fp_prob = fopen(ss.str().c_str(), "w");
    }

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
	    L_ = L_ + alpha_*(L_data/N_ - L_model/M_);
	    symmetrization_U<double>(L_);
	    alpha_ = 1.0/(t+1);

	    //      std::cout << L_ << std::endl;
	    // print KL-divergence
	    for (i=0; i<all_v_.size(); i++){
		    prob_P_[i] = ProbModelAssignV(all_v_[i]);
	    }

	    fprintf(fp		, "%d %lf\n", t, KLd(prob_P_, prob_Q_));
	    // print L
	    printLayer(	fp_data	, L_data ); fprintf(fp_data , "\n");
	    printLayer(	fp_model, L_model); fprintf(fp_model, "\n");
	    printLayer(	fp_L	, L_  )   ; fprintf(fp_L    , "\n");

	    fprintf(fp_prob, "#t\tQ(v)\tP(v)\n");
	    for(i=0; i<all_v_.size(); i++){
		    fprintf(fp_prob, "%d\t%8.4lf\t%8.4lf\n", t, prob_Q_[i], prob_P_[i]);
	    }
    }
    fclose(fp      );
    fclose(fp_data );
    fclose(fp_model);
    fclose(fp_L    );
    fclose(fp_prob );
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

  double Energy(Pattern &v)
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

  double LinearSum(int i, Pattern& vf)
  {
    double u_i=0.0;
    for(int j=0; j<D_; j++)
      u_i += L_(i,j)*vf(j) - a_[i]*vf(j);    
    return u_i;
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
  	u += L_(i,k)*v(i)*v(k);
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


//可視層、隱れ層を持つ全結合のボルツマンマシン
//*バイアス項については考えないでおく*

int main(int argc, char **argv)
{
  size_t i,j,k,m,n;
  int result;
  size_t D =  3, P = 1;
  //size_t K = 5, M = 1000, N=1000, T=1000;
  size_t K = 5, M = 1000, N=1000, T=10;
  double alpha = 1.0;
  std::vector<Pattern> all_v(0x01<<D, Pattern(D,0));
  std::vector<Pattern> all_h(0x01<<P, Pattern(P,0));
  for(i=0; i < (0x01<<D); i++){
    all_v[i] = dtob(i,D);
//    all_v[i](0) = 1;
  }
  for(i=0; i < (0x01<<P); i++){
    all_h[i] = dtob(i,P);
  }
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

  double rnd;
  Layer W(D,D);
  Probability prob_P2(0x01<<D,0);
  Probability cumulative_Q_(0x01<<D, 0);
  for(i=1, cumulative_Q_[0] = prob_Q[0]; i<0x01<<D; i++)
	  cumulative_Q_[i] = cumulative_Q_[i-1] + prob_Q[i];
  for(n=0; n<N; n++){
	  //	  rnd = drand48();
	  //	  for(i=0; i<cumulative_Q_.size(); i++)
	  //		  if (rnd < cumulative_Q_[i])
	  //			  break;
	  for(i=0; i < 0x01<<D; i++)
		  if( n < (int)(N*cumulative_Q_[i]))
			  break;
	  Pattern& vi = all_v[i];
	  prob_P2[i]++;
	  W += ub::outer_prod(vi,vi);
	  std::cout << vi << std::endl;
	  printLayer(stdout,W);
	  std::cout << std::endl;
	  symmetrization_U<double>(W);
	  printLayer(stdout,W);
	  std::cout << "----------------------" << std::endl;
  }
  W /= N;
  printLayer(stdout,W);

  std::cout << std::endl;
  for(i=0; i<0x01<<D; i++) std::cout << prob_P2[i]/N << "\t";
  std::cout << std::endl;
  for(i=0; i<0x01<<D; i++) std::cout << prob_Q [i] << "\t";
  std::cout << std::endl;
  //while ((result = getopt(argc, argv, "a:d:p:k:m:n:t:f:")) != -1){
  //  switch(result){
  //  case 'a':
  //    alpha = atof( optarg );
  //  case 'd':
  //    //      D = atoi( optarg );
  //    std::cout << "visible layer degree D is const." << std::endl;
  //    break;
  //  case 'p':
  //    P = atoi( optarg );
  //    break;
  //  case 'k':
  //    K = atoi( optarg );
  //    break;
  //  case 'm':
  //    M = atoi( optarg );
  //    break;
  //  case 'n':
  //    N = atoi( optarg );
  //    break;
  //  case 't':
  //    T = atoi( optarg );
  //    break;
  //  case 'f':
  //    //F = atoi( optarg );
  //    break;
  //  case '?':
  //    printf("unknown\n");
  //  }
  //}

  //BM *fbm1 = new FBMBias  (D,K, M, N, T, alpha, prob_Q, all_v);
  BM *fbm2 = new FBMNoBias(D,K, M, N, T, alpha, prob_Q, all_v);
  //  BM *rbm = new RBM(D, P, K, M, N, T, alpha, prob_Q, all_v, all_h);
  //  BM *gbm = new GBM(D, P, K, M, N, T, alpha, prob_Q, all_v, all_h);

  //fbm1->Learn();
  fbm2->Learn();
  //  rbm->Learn();
  //  gbm->Learn();
  return 0;
}

