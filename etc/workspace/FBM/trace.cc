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

template<typename T>
void pmat(T& m)
{
  for (size_t i=0; i < m.size1(); i++){
    ub::matrix_row<T> mr(m,i);
    std::cout << mr << std::endl;
  }
  std::cout << std::endl;
}

double trace(Layer& W)
{
  using namespace boost::numeric::ublas;
  size_t D = W.size1();
  matrix_vector_range<Layer> diag(W,range(0,D),range(0,D));
  return sum(diag);
}

int main(int argc, char **argv)
{
  size_t D = 4;
  Layer W(D,D);

  for (size_t i=0; i<D; i++){
    for (size_t j=0; j<D; j++){
      W(i,j) = 10*i+j;
    }
  }

  ub::matrix_vector_range<Layer> diag(W, ub::range(0,D),ub::range(0,D));
  std::cout << "trace=" << ub::sum(diag) << std::endl;
  std::cout << "trace=" << trace(W) << std::endl;
  return 0;
}
