#include "FBM.hh"

Pattern dtob_date(int n, size_t size)
{
  Pattern ret(size, 0);
  for (size_t i=0; i < size-1; i++)
    ret[size-1-i] = (n>>i) & 0x01;
  return ret;
}

int btod_date(const Pattern& x)
{
  size_t i, ret=0, size=x.size();
  for(i=0; i < size-1; i++)
    ret += x[size-1-i]*(0x01<<i);
  return ret;
}

#define DENAME(id,i,v,L) \
  std::cout << id << "_"<< i << "(" << v << "," << L << ")" << std::endl;

#define verb(w,s,i,j){						\
    printf("%s(%d,%d)=%4.1lf\t",#w,i,j,w(i,j));			\
    printf("%s(%d)=%4.1lf\t"   ,#s,j,s[j]);			\
    printf("%4.1lf*%4.1lf=%4.1lf\t",w(i,j),s[j],w(i,j)*s[j]);	\
    printf("u_%d=%4.1lf\n",i,u_i);				\
  }

#define restp(L,i,j,u_i)					\
  printf("+=%s(%d,%d)=%4.1lf\tu_%d=%lf\n\n",#L,i,j,L(i,j),i,u_i);

double _denergy0_verb(size_t i, Pattern& v, Layer& L)
{
  double u_i=0.0;
  DENAME("dEn0", i, v, L);
  for (size_t j=1; j<v.size(); j++){
    u_i += L(i,j)*v[j];
    verb(L,v,i,j);
  }
  u_i += L(0,i);
  restp(L,0,i,u_i);
  return u_i;
}

double _denergy02(size_t i, Pattern& v, Layer& L)
{
  double u_i=0.0;
  for (size_t j=0; j<v.size(); j++){
    u_i += L(i,j)*v[j];
  }
  return u_i;
}

double _denergy02_verb(size_t i, Pattern& v, Layer& L)
{
  double u_i=0.0;
  for (size_t j=0; j<v.size(); j++){
    u_i += L(i,j)*v[j];
    verb(L,v,i,j);
  }
  printf("dE2_%d=%4.1lf\n\n",i,u_i);
  return u_i;
}


double _denergy4_verb(size_t i, Pattern& v, Layer& L)
{
  DENAME("dEn4",i,v,L);
  double u_i=0.0;
  for (size_t j=0; j<i; j++){
    u_i += L(i,j)*v[j];
    verb(L,v,i,j);
  }
  for (size_t j=i+1; j<v.size(); j++){
    u_i += L(i,j)*v[j];
    verb(L,v,i,j);
  }
  u_i += L(i,i);// bias term
  restp(L,i,i,u_i);
  return u_i;
}

int main(int argc, char* argv[])
{
  size_t D=3, P=3;
  size_t K=5, M=1000, N=1000, T=10000;
  char prefix[32]="energy", odir[32]="dir-dat";
  g_energy_idx  = 10;
  g_denergy_idx =  4;
  g_prop_idx    = -1;
  g_annealing_schedule_idx = 1;

  std::ofstream ofs_L, ofs_t_gibbs, ofs_cfg;
  cook_args(argc, argv, D, P, K, M, N, T, prefix, odir);

  Layer L_k(D  ,D  );
  symmetric_view sL_k(L_k);
  banded_view bL_k(L_k);
  Patterns all_v_k(0x01<<D,Pattern(D,0));

  Layer L_d(D+1,D+1);
  symmetric_view sL_d(L_d);
  banded_view bL_d(L_d);
  Patterns all_v_d(0x01<<(D+1), Pattern(D+1,0));

  for (int i=0; i<0x01<<D; i++){
    all_v_k[i] = dtob(i,D);
    all_v_d[i] = dtob_date(i,D+1);
    all_v_d[i][0] = 1;
  }

  for (size_t i=1; i<L_d.size1(); i++){
    for (size_t j=i; j<L_d.size2(); j++){
//       L_k(i-1,j-1) = L_k(j-1,i-1) = i*10 + j;
//       L_d(i  ,j  ) = L_d(j  ,i  ) = i*10 + j;
      double r = drand48();
      L_k(i-1,j-1) = L_k(j-1,i-1) = r;
      L_d(i  ,j  ) = L_d(j  ,i  ) = r;
    }
    L_d(0,i) = L_d(i,0) = L_d(i,i);
  }  
  L_d = L_d - bL_d;

  printf("kurata matrix\n");
  for (size_t i=0; i<L_k.size1(); i++){
    for (size_t j=0; j<L_k.size2(); j++)
      printf("%4.1lf\t", L_k(i,j));
    printf("\n");
  }

  printf("\n");
  printf("date matrix\n");
  for (size_t i=0; i<L_d.size1(); i++){
    for (size_t j=0; j<L_d.size2(); j++)
      printf("%4.1lf\t", L_d(i,j));
    printf("\n");
  }


  std::cout << "En_d =\t";
  for (int i=0; i < 0x01<<D; i++){
    Pattern &v_d = all_v_d[i];
    std::cout << _energy0(v_d, L_d) << "\t";
  }
  std::cout << std::endl;

  std::cout << "En_k3 =\t";
  for (int i=0; i < 0x01<<D; i++){
    Pattern &v_k = all_v_k[i];
    std::cout << _energy3 (v_k, L_k) << "\t";
  }
  std::cout << std::endl;

//   std::cout << "En_k4 =\t";
//   for (int i=0; i < 0x01<<D; i++){
//     Pattern &v_k = all_v_k[i];
//     std::cout << _energy4 (v_k, L_k) << "\t";
//   }
//   std::cout << std::endl;

  std::cout << "En_k8 =\t";
  for (int i=0; i < 0x01<<D; i++){
    Pattern &v_k = all_v_k[i];
    std::cout << _energy8 (v_k, L_k) << "\t";
  }
  std::cout << std::endl;

//   std::cout << "En_k9 =\t";
//   for (int i=0; i < 0x01<<D; i++){
//     Pattern &v_k = all_v_k[i];
//     std::cout << _energy9 (v_k, L_k) << "\t";
//   }
//   std::cout << std::endl;

  std::cout << "En_k10=\t";
  for (int i=0; i < 0x01<<D; i++){
    Pattern &v_k = all_v_k[i];
    std::cout << _energy10(v_k, L_k) << "\t";
  }
  std::cout << std::endl;




  for (int i=0; i < D; i++){
    std::cout << "i=\t" << i << std::endl;    
    std::cout << "dEn_d =\t";
    for (size_t n=0; n<0x01<<D; n++){
      Pattern &v_d = all_v_d[n];
      std::cout << _denergy0(i+1, v_d, L_d) << "\t";
    }
    std::cout << std::endl;

    std::cout << "dEn_k1=\t";
    for (size_t n=0; n<0x01<<D; n++){
      Pattern &v_k = all_v_k[n];
      std::cout << _denergy1(i, v_k, L_k) << "\t";
    }
    std::cout << std::endl;

//     std::cout << "dEn_k2=\t";
//     for (size_t n=0; n<0x01<<D; n++){
//       Pattern &v_k = all_v_k[n];
//       std::cout << _denergy2(i, v_k, L_k) << "\t";
//     }
//     std::cout << std::endl;

//     std::cout << "dEn_k3=\t";
//     for (size_t n=0; n<0x01<<D; n++){
//       Pattern &v_k = all_v_k[n];
//       std::cout << _denergy3(i, v_k, L_k) << "\t";
//     }
//     std::cout << std::endl;

    std::cout << "dEn_k4=\t";
    for (size_t n=0; n<0x01<<D; n++){
      Pattern &v_k = all_v_k[n];
      std::cout << _denergy4(i, v_k, L_k) << "\t";
    }
    std::cout << std::endl;
  }

  for (int i=0; i < D; i++){
    //    std::cout << "dEn_" << i << "=\n";
    std::cout << "i =" << i << std::endl;
    for (size_t n=0; n<0x01<<D; n++){
      Pattern &v_d = all_v_d[n];
      Pattern &v_k = all_v_k[n];
      std::cout << "v^" << n << "=" << v_d << std::endl;
      std::cout << "v^" << n << "=" << v_k << std::endl;
      std::cout << "En0=" <<_denergy0(i+1, v_d, L_d) << std::endl;
      std::cout << "En4=" <<_denergy4(i  , v_k, L_k) << std::endl;
      //      _denergy0_verb (i+1, v_d, L_d);
      //      _denergy4_verb (i  , v_k, L_k);
//       printf("dEn02_%d=%4.1lf\n"  ,i,_denergy02(i,v_d,L_d));
//       printf("dEn4k_%d=%4.1lf\n\n",i,_denergy4 (i,v_k,L_k) );
    }
    std::cout << std::endl;
  }
  

  return 0;
}
