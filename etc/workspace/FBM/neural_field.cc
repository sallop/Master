#include "FBM.hh"

// void set_boltzmann_neural_field(Layer &W,int n, int k, double r)
// {
//   for (int i=0; i<n; i++){
//     for (int j=i; j<n; j++){
//       W(i,j) = -r;
//       W(j,i) = -r;
//     }
//     W(i,i) = r*(k-0.5);
//   }
// }

void set_boltzmann_neural_field(Layer &W,int n, int k, double r)
{
  double b = 5;
  for (int i=0; i<n; i++){
    for (int j=i+1; j<n; j++)
      W(i,j) = W(j,i) = -r;
    for (int l=i-k; l<i+k; l++){
      int j = (l<0) ? l + n : l % n;
      W(i,j) = W(j,i) = -r*(1.0 - b);
    }
    W(i,i) = r*(k-0.5);
  }
}

int main(int argc, char *argv[])
{
  int n=16, k=atoi(argv[1]);
  double r = 100;
  Layer W(n,n);

  set_boltzmann_neural_field(W,n,k,r);

  FILE *gp = popen("gnuplot","w");
  //  FILE *gp = stdout;
  fprintf(gp,"set grid\n");
  fprintf(gp,"set xrange[%d:-1.5]\n",n);
  //fprintf(gp,"set cbrange[0.0:1.0]\n");
  fprintf(gp,"plot '-' matrix with image\n");
  for (int i=0; i<n; i++){
    for (int j=0; j<n; j++)
      fprintf(gp,"% 8.4lf ", W(i,j));
    fprintf(gp,"\n");
  }
  fprintf(gp,"e\ne\n");
  fprintf(gp,"pause %s\n", argv[2]);
  fclose(gp);


  return 0;
}
