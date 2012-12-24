#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
//#include <boost/numeric/ublas/diagonal.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric::ublas;  
typedef vector<double> dvector;
typedef matrix<double> dmatrix;

template<typename T>
void print(T& m)
{
  for (size_t i=0; i < m.size1(); i++){
    matrix_row<T> mr(m,i);
    std::cout << "i=" << i << mr << std::endl;
  }
}

double prop(int i, dvector& v, dmatrix& L)
{
  matrix_row<dmatrix> Li(L,i);
  double term = inner_prod(v,Li);
  return term;
}

int main()
{
  size_t D = 8, P = 4;
  dmatrix W(D,P);
  dmatrix L(D,D), J(P,P);
  dvector v(D), h(P);
  for(size_t i=0; i < D; i++){
    for(size_t j=0; j < P; j++){
      W(i,j) = (i+1)*10 + (j+1);
    }
  }

  for (size_t i=0; i<D; i++){
    for (size_t j=0; j<D; j++){
      L(i,j) = (i+1)*10 + (j+1);
    }
  }
  for (size_t i=0; i<P; i++){
    for (size_t j=0; j<P; j++){
      J(i,j) = (i+1)*10 + (j+1);
    }
  }

  //std::cout << inner_prod(v,m) << std::endl;
  //unit_vector<double> v(P,i);
  std::cout << "W=" << W << std::endl;
  for (size_t i=0; i < W.size1(); ++i){
    matrix_row<dmatrix> Wr(W,i);
    std::cout << Wr << std::endl;
  }


  for (size_t i=0; i < W.size1(); ++i){
    v = unit_vector<double>(D,i);
    for (size_t j=0; j < W.size2(); ++j){
      h = unit_vector<double>(P,j);
      std::cout << "v W = " << prod(v,W) << std::endl;
      //error!      std::cout << "W v = " << prod(W,v) << std::endl;
      std::cout << "W h = " << prod(W,h) << std::endl;
      //error! std::cout << "h W = " << prod(h,W) << std::endl;
      std::cout << "(v W) h = " << inner_prod(prod(v,W),h) << std::endl;
    }
  }

  symmetric_adaptor<dmatrix,lower> sL(L);
  std::cout << sL << std::endl;
  for (size_t i=0; i < sL.size1(); i++){
    matrix_row<symmetric_adaptor<dmatrix,lower> > Lr(sL,i);
    std::cout << "i=" << i << "\t" << Lr << std::endl;
  }

  for (dmatrix::iterator1 row = L.begin1(); row != L.end1(); ++row)
    {
      std::cout << *row << std::endl;
    }
  for (dmatrix::iterator2 col = L.begin2(); col != L.end2(); ++col)
    {
      std::cout << *col << std::endl;
    }
  // for (dmatrix::iterator1 r=sL.begin1();
//        r != r.end(); ++r){
//     std::cout << *r << std::endl;
//   }

  for(dvector::iterator it=v.begin(); it != v.end(); ++it){
    *it = 9999;
  }

  std::cout << v << std::endl;


  banded_adaptor<dmatrix> bm(L,0,0);
  std::cout << bm << std::endl;
  for (size_t i=0; i < bm.size1(); i++){
    matrix_row<banded_adaptor<dmatrix> > mr(bm,i);
    std::cout << "i=" << i << "\t" << mr << std::endl;
  }
  //diagonal_matrix<dmatrix> dL(L);
  //std::cout << dL << std::endl;

  banded_adaptor<symmetric_adaptor<dmatrix,lower> > bm2(sL,0,0);
  std::cout << bm2 << std::endl;




  print<dmatrix>(L);
  print<banded_adaptor<dmatrix> >(bm);
  
  return 0;
}
