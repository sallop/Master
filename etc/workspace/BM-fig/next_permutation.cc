#include <iostream>
#include <cstdio>
#include <vector>
#include <algorithm>

using namespace std;
typedef std::vector<int> ivector;
typedef std::vector<double> fvector;

int main()
{
  const int D = 3;
  size_t i,j,k;
  //  ivector data(0x01 << D, 0);
  //data.push_back(i);
// 168 pattern
//   data[0] = 1;
//   data[1] = 1;
//   data[2] = 1;
//   data[3] = 1;
//   data[4] = 1;
//   data[5] = 2;
//   data[6] = 2;
//   data[7] = 3;

  fvector data(0x01<<D, 0);
  data[0] = 0.05;
  data[1] = 0.05;
  data[2] = 0.10;
  data[3] = 0.10;
  data[4] = 0.10;
  data[5] = 0.10;
  data[6] = 0.10;
  data[7] = 0.40;





  //  j = 0;
  do{
    //    printf("%d\t",j);
//     printf("%d %d %d %d ",
// 	   data[0],data[1],data[2],data[3]);
//     printf("%d %d %d %d\n",
// 	   data[4],data[5],data[6],data[7]);
    printf("%4.4lf\t%4.4lf\t%4.4lf\t%4.4lf\t",
	   data[0],data[1],data[2],data[3]);
    printf("%4.4lf\t%4.4lf\t%4.4lf\t%4.4lf\n",
	   data[4],data[5],data[6],data[7]);
    j++;
  } while (next_permutation(data.begin(), data.end()));

  return 0;
}
