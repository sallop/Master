#include <stdio.h>

#define ENE(n)	energy##n
#define TE2(n)  TE1(n)
#define xstr(s)	str(s)
#define str(s)	#s
#define E(n)    xstr(ENE(n))
#define E2(n)   xstr(energy##n)
#define E3(n)   str(energy##n)


xstr(ENE(1))
str(ENE(2))

E(3)
E(4)

E2(9)
E3(8)
double _f0(double x){ return x*0;}
double _f1(double x){ return x*1;}
double _f2(double x){ return x*2;}
double _f3(double x){ return x*3;}



int main()
{
  typedef double (*func)(double);
  int i;
  func f;
  func fs[] = {_f0, _f1, _f2, _f3};//, NULL};
  //for (f = fs[0]; f != NULL; f++){
  for (i=0; i<sizeof(fs)/sizeof(fs[0]); i++){
    f = fs[i];
    printf("%lf\n", f(1.0));
  }
  return 0;
}
