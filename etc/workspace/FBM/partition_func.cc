#include <iostream>
#include <cmath>
#include <algorithm>
#include <iterator>
#include <boost/fuction.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>


using namespace std;
using namespace boost::lambda;

typedef vector<double> dvector;
typedef vector<double> dvectors;

void dtob(dvector& v, int n)
{
  for (size_t i=0; i < v.size(); ++i)
    v[i] = (n>>i) & 0x01;
}

int main()
{
  int D = 3;
  dvectors vs(0x01<<D, dvector(D,0));
  for (dvectors::iterator it=dvectors.begin();
       it != dvectors.end(); it++)
    {
      
    }
}
