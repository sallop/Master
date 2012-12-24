#include "Layer.hh"
#include "GBM.hh"

double energy(Layer& W, Pattern &v)
{
  double ret=0.0;
  for (size_t i=0; i<W.size1(); i++)
    for (size_t j=i; j<W.size2(); j++)
      ret += W(i,j)*v[i]*v[j];
  return -ret;
}

int main()
{
  size_t D = 16;
  size_t K = 5;
  Layer W(D,D);
  W = zerom(D,D);
  Patterns all_v(D,Pattern(D,0));
  Patterns local_minimum_patterns(D,Pattern(D,0));
  std::vector<int> local_minimal_values(D,0);
  for (size_t i=0; i<D; i++){
    Pattern& v = all_v[i];
    for (size_t k=0; k<K; k++){
      v[(i+k)%D]=1;
    }
    local_minimal_values[i] = btod(v);
    local_minimum_patterns[i] = v;
  }

  for (size_t i=0; i<D; i++){
    std::cout << local_minimal_values[i] << std::endl;
  }

//   for (size_t i=0; i<D; i++){
//     Pattern& v = all_v[i];
//     W += ub::outer_prod(v,v);
//     std::cout << "W =" << W << std::endl;
//   }
  for (size_t i=0; i<W.size1(); i++){
    Pattern v = all_v[i];
    for (size_t j=i+1; j<W.size2(); j++){
      //W(i,j) = W(j,i) = drand48();
      //W(i,j) = 10*drand48();
      W(i,j)=0.0;
    }
    W += ub::outer_prod(v,v);
    W(i,i)=0.0;
  }
  

  FILE *gp=popen("gnuplot","w");
  // fprintf(gp,"set yrange [:] reverse\n");
//   fprintf(gp,"plot '-' matrix with image\n");
//   for (size_t i=0; i<D; i++){
//     for (size_t j=0; j<D; j++){
//       fprintf(gp,"%lf ", W(i,j));
//     }
//     fprintf(gp,"\n");
//   }
//   fprintf(gp,"e\n");
//   fprintf(gp,"e\n");
//   fprintf(gp,"pause 0.5\n");

  double max_pattern = 0x01<<D;
  Pattern v(D,0);
  fprintf(gp,"set xrange [0:%lf]\n",max_pattern);
  fprintf(gp,"set yrange [%d:0]\n",-100);
  fprintf(gp,"set terminal postscript enhanced color eps\n");
  fprintf(gp,"set output 'energy_gray.eps'\n");
  //  fprintf(gp,"set multiplot layout %d,%d\n",2,1);
  fprintf(gp,"set multiplot\n");
  // fprintf(gp,"plot '-' with lines title 'i++'\n");
//   for (size_t i=0; i<max_pattern; i++){
//     size_t j=i;
//     //size_t j=i^(i>>1);
//     v = dtob(j,D);
//     //std::cout << "v(" << i << ")=" << v << std::endl;
//     fprintf(gp,"%d %lf\n",i,energy(W,v));
//   }
//   fprintf(gp,"e\n");

  fprintf(gp,"plot '-' with lines title 'gray code i^(i>>1)'\n");
  for (size_t i=0; i<max_pattern; i++){
    size_t j=i^(i>>1);
    //size_t j=i;
    v = dtob(j,D);
    fprintf(gp,"%d %lf\n",i,energy(W,v));
  }
  fprintf(gp,"e\n");
  fprintf(gp,"plot '-' with points ls 3\n");
  for (size_t i=0; i<D; i++){
    //size_t j = local_minimal_values[i];
    v = local_minimum_patterns[i];
    size_t j = btod(v);
    fprintf(gp,"%d %lf\n",j,energy(W,v));
  }
  fprintf(gp,"e\n");
  fprintf(gp,"pause 3.0\n");

  pclose(gp);
}
