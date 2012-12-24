#include <iostream>
#include <vector>
#include <boost/algorithm/combination.hpp>

typedef std::vector<int>    ivector;
typedef std::vector<double> dvector;
using namespace std;

int main()
{
  const int r = 3;
  const int n = 10;

  ivector v(n);
  for(size_t i=0; i<n; ++i) v[i] = i;

  int N=0;
  do{
    ++N;
    if (N<10 || N>117){
      cout << "[" << v[0];
      for (size_t j=1; j<r; ++j)
	cout << "N=" << N << "," << v[j];
      cout << "]" << endl;
    } else if (N == 10) {
      cout << "..." << endl;
    }
  } while (next_combination(v.begin(), v.begin() + r, v.end()));
  return 0;
}
