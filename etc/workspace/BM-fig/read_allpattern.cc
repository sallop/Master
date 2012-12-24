#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cstdio>
#include <cstdlib>

using namespace std;
#define DIMENSION  3
#define ALLPATTERN (0x01<<DIMENSION)


bool read_prob(double prob[ALLPATTERN])
{
  static FILE *fp = fopen("all_pattern.dat","r");
  static char str[128];
  static const char * delim = "\t";
  size_t i;
  while ((fgets(str, 128-1, fp)) != NULL){
    const char *token = strtok(str, delim);
    for (i=0;
	 token;
	 token = strtok(0, delim), i++)
      {
	//cout << token << endl;
	prob[i] = atof(token);
      }
    return true;
  }
  return false;
}


int main()
{
  double prob_Q[ALLPATTERN];
  
  
  int j=0;
  while(read_prob(prob_Q)){
    printf("j = %d\n", j);
    for (int i=0; i<ALLPATTERN; i++){
      printf("p(x^%d) = %lf\n", i, prob_Q[i]);
    }
    printf("\n");
    j++;
  }
  return 0;
}
