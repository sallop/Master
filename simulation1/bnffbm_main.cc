#include "BM.hh"
#include "GBM.hh"
#include "config.hh"

#define TEMPLATE 1.0
//#define DEBUG
// base program is GBM

template<typename MATRIX> void
set_boltzmann_neural_field(MATRIX &W,int n, int k, double r)
{
  // n: number of units
  // k: k-near bit is on
  // r: if r -> \infty then get a boltzmann distribution.
  double b = 10;
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

template<typename MATRIX> void
plot_matrix(FILE *gp, MATRIX &W, const char* title)
{
  fprintf(gp, "plot '-' matrix with image title '%s'\n", title);
  for (size_t i=0; i < W.size1(); i++){
    for (size_t j=0; j < W.size2(); j++)
      fprintf(gp, "%lf ", W(i,j));
    fprintf(gp, "\n");
  }
  fprintf(gp, "e\n");
  fprintf(gp, "e\n");
}

#define _s(x) x,#x

FILE* init_gnuplot()
{
  FILE *gp = popen("gnuplot -persist","w");
#ifdef DEBUG
  fprintf(gp, "set terminal x11\n");
#else
  fprintf(gp, "set terminal gif\n");
#endif
  fprintf(gp, "set palette defined (0 'white', 1 'blue')\n");
  fprintf(gp, "unset key\n");
  fprintf(gp, "unset xtics\n");
  fprintf(gp, "unset ytics\n");
  //fprintf(gp, "unset cbtics\n");
  fprintf(gp, "set yrange [:] reverse\n");
  return gp;
}

FILE* plot_status(FILE *gp, char *prefix, int time,
		  SubLayer& L, Layer& La, Layer& Ld, Layer& Lm,
		  SubLayer& W, Layer& Wa, Layer& Wd, Layer& Wm,
		  SubLayer& J, Layer& Ja, Layer& Jd, Layer& Jm)
{
  fprintf(gp, "set output 'dir-dat/%s-LWJ-%d.gif'\n",prefix,time);
  fprintf(gp, "set multiplot layout 3,4\n");
  // L
  fprintf(gp, "set title 'L'\n");
  plot_matrix<SubLayer>(gp,_s(L ));
  fprintf(gp, "set title 'La'\n");
  plot_matrix<Layer>   (gp,_s(La));
  fprintf(gp, "set title 'Ld'\n");
  plot_matrix<Layer>   (gp,_s(Ld));
  fprintf(gp, "set title 'Lm'\n");
  plot_matrix<Layer>   (gp,_s(Lm));
  // W
  
  //  fprintf(gp, "set cbrange[-1:1]\n");
  fprintf(gp, "set title 'W'\n");
  plot_matrix<SubLayer>(gp,_s(W ));
  //  fprintf(gp, "set cbrange[*:*]\n");
  fprintf(gp, "set title 'Wa'\n");
  plot_matrix<Layer>   (gp,_s(Wa));
  fprintf(gp, "set title 'Wd'\n");
  plot_matrix<Layer>   (gp,_s(Wd));
  fprintf(gp, "set title 'Wm'\n");
  plot_matrix<Layer>   (gp,_s(Wm));


  // J
  fprintf(gp, "set title 'J'\n");
  plot_matrix<SubLayer>(gp,_s(J ));
  fprintf(gp, "set title 'Ja'\n");
  plot_matrix<Layer>   (gp,_s(Ja));
  fprintf(gp, "set title 'Jd'\n");
  plot_matrix<Layer>   (gp,_s(Jd));
  fprintf(gp, "set title 'Jm'\n");
  plot_matrix<Layer>   (gp,_s(Jm));
  fprintf(gp, "unset multiplot\n");

  fprintf(gp, "set title 'time = %d'\n");
  return gp;
}

FILE* plot_status(FILE *gp, char *prefix, int time, Layer& F)
{
  fprintf(gp, "set output 'dir-dat/%s-F-%d.gif'\n",prefix,time);
  fprintf(gp, "set title 'F'\n");
  plot_matrix<Layer>(gp,_s(F));
  return gp;
}


typedef double (*transition_func)(double alpha, double prob);

double transition_on (double alpha, double q){ return alpha*q;}
double transition_off(double alpha, double q){ return     0.0;}
transition_func g_transition[2][2] = {
  {transition_off, transition_off},
  {transition_off, transition_on },
};


int
main(int argc, char **argv)
{
  int D=16, P=16;
  int K=256, M=10, N=1000, T=10000;

  Layer F(D+P,D+P);

  char prefix[32]="bnffbm", odir[32] = "dir-dat";
  double connection_scale = 1.0;
  g_annealing_schedule_idx = 1;
  g_boltzmann_neural_field_idx = 0;
  g_denergy_idx    = 0;
  g_energy_idx     = 0;
  g_init_layer_idx = 1;
  g_prop_idx       = 0;

  srand(g_seed);
  srand48(g_seed);

  std::fstream ofs_F, ofs_t_gibbs, ofs_cfg;

  cook_args(argc, argv, D, P, K, M, N, T,
	    connection_scale, prefix, odir);

  init_plot_table(odir, prefix,
		  ofs_F,
		  ofs_t_gibbs,
		  ofs_cfg,
		  std::ios_base::out);

  banded_view Fb(F,0);
  SubLayer L(F ,ub::range(0,D  ),ub::range(0,D  ));
  SubLayer W(F ,ub::range(0,D  ),ub::range(D,D+P));
  SubLayer J(F ,ub::range(D,D+P),ub::range(D,D+P));
  Layer La(D,D), Ld(D,D), Lm(D,D);
  Layer Wa(D,P), Wd(D,P), Wm(D,P);
  Layer Ja(P,P), Jd(P,P), Jm(P,P);
  const Layer L0 = zerom(D,D);
  const Layer W0 = zerom(D,P);
  const Layer J0 = zerom(P,P);
  symmetric_view Fs(F);
  symmetricsub_view Ls(L), Js(J);

  Probability prob_Q (D,0); // environment
  Probability accm_Q (prob_Q.size(),0);
  Probability boltz_P(prob_Q.size(),0);
  Probability gibbs_P(D,0); // can't print KLd. so alternate vf mean
  Probability gibbs_P0 = zerov(gibbs_P.size());
  //double alpha = 0.75;	     // learning rate
  double alpha       = g_init_alpha; // 
  double temperature = 1.0; // temperature
  double kl;// learning rate, store the kl divergence
  
  Pattern v(D,0), h(P,0);
  const Pattern v0(D,0), h0(P,0);
  Pattern tmp_v(D,0), tmp_h(P,0);
  Pattern pv(D,0), aqv(D,0);
  Pattern ph(P,0), aqh(P,0);
  // all state
  Patterns all_v (D,Pattern(D,0));
  Patterns all_h (P,Pattern(P,0));

  Patterns all_mu(all_v.size(),Pattern(P,0));
  DtoB_table dtob_tbl;
  // fantasy particle
  Patterns all_vf(M,Pattern(D,0));// fantasy particle
  Patterns all_hf(M,Pattern(P,0));// fantasy particle

  std::fill(prob_Q.begin(), prob_Q.end(), 1./prob_Q.size());
  std::partial_sum(prob_Q.begin(),prob_Q.end(),accm_Q.begin());

  std::cout << _eq(prob_Q) << std::endl;
  std::cout << _eq(accm_Q) << std::endl;

  assert(accm_Q[accm_Q.size()-1] > 0.99 &&
	 accm_Q[accm_Q.size()-1] < 1.01);
  accm_Q[accm_Q.size()-1];
  for (size_t i=0; i<all_v.size(); i++){
    Pattern &v = all_v[i];
    for (size_t j=0; j<g_neighbor_k; j++){
      v[(i+j)%D]=ONBIT;
    }
    dtob_tbl[btod(v)] = make_pair(i,v);
  }
  for (size_t i=0; i<all_h.size(); i++){
    Pattern &h = all_h[i];
    for (size_t j=0; j<g_neighbor_k; j++){
      h[(i+j)%P]=ONBIT;
    }
  }

  for (size_t i=0; i<all_v.size(); i++){
    std::cout << _eq(i) << " " << all_v[i] << std::endl;
    std::cout << _eq(i) << " " << all_h[i] << std::endl;
  }

 
  //1. Randomly initialize mu and
  F = init_layer(D+P,connection_scale);
  L = L0;
  //  W = W0;
  if (g_boltzmann_neural_field_idx == 0){
    J = J0;
  } else if(g_boltzmann_neural_field_idx == 1){
    set_boltzmann_neural_field<SubLayer>(J, J.size1(), g_neighbor_k-1, 10);
  } else {
    std::cout << "Not support flag" << std::endl;
    std::cout << _eq(g_boltzmann_neural_field_idx) << std::endl;
    assert(false);
  }
  // n: number of units
  // k: k-bit is on
  // r: if r -> \infty then get a boltzmann distribution.
  //set_boltzmann_neural_field(F, F.size1(), 2, 16); -> wrong
  //  set_boltzmann_neural_field<SubLayer>(J, J.size1(), 2, 16);
  
  std::cout << _eq(F) << std::endl;
  std::cout << _eq(L) << std::endl;
  std::cout << _eq(W) << std::endl;
  std::cout << _eq(J) << std::endl;

  std::cout << _eq(alpha) << std::endl;
  std::cout << _eq(g_boltzmann_neural_field_idx) << std::endl;
  std::cout << _eq(g_neighbor_k) << std::endl;

  std::cout << _eq(argc) << std::endl;
  std::cout << _eq(argv[1]) << std::endl;
  //  assert(false);
  for(size_t m=0; m<M; m++){
    for(size_t i=0; i<D; i++)
      all_vf[m][i] = drand48()>0.5?ONBIT:OFFBIT;
    for(size_t i=0; i<P; i++)
      all_hf[m][i] = drand48()>0.5?ONBIT:OFFBIT;
  }
  
  {// plot configuration
    symcfg( D, P, K, M, N, T,
	    prefix,
	    F,L,W,J,
	    prob_Q,accm_Q,
	    alpha,
	    all_v,all_h,all_mu,
	    connection_scale,
	    ofs_cfg);
  }

  symcfg( D, P, K, M, N, T,
	  prefix,F,L,W,J,
	  prob_Q,accm_Q,alpha,
	  all_v,all_h,all_mu,
	  connection_scale,std::cout);

  //FILE *gp = popen("gnuplot -persist","w");
  FILE *gp = init_gnuplot();
  // Full connected boltzmann machine kurata9_main.cc


  for(size_t t=0; t<T; t++){
    double Z;
    Ld = Lm = L0;
    Wd = Wm = W0;

    {// hebb learning
      v = all_v[rand()%D];
      // upper side probability
      Z = 0.0;
      for (size_t j=0; j<P; j++){
	ph[j] = exp(-denergy_up(j,v,all_h[j],L,J,W));
	Z += ph[j];
      }
      ph /= Z;
      //std::cout << __LINE__ << "-" << _eq(ph) << std::endl;
      //std::cout << __LINE__ << "-"
      //	      << std::accumulate(ph.begin(), ph.end(), 0.0) << std::endl;
      //assert(false);
    
      // down side probability
      Z = 0.0;
      for (size_t i=0; i<D; i++){
	pv[i] = exp(-denergy_dw(i,all_v[i],h,L,J,W));
	Z += pv[i];
      }
      pv /= Z;

      // down to up
      for (size_t j=0; j<P; j++){
	Pattern &hj = all_h[j];
	for (size_t k=0; k < D; k++){
	  for (size_t l=0; l < P; l++){
	    Wd(k,l) += g_transition
	      [static_cast<int>(hj[k])]
	      [static_cast<int>(hj[l])]
	      (alpha,ph[j]);
	  } // l
	}	// k
      } // j

      //    std::cout << "Wd" << std::endl;
      //assert(false);
    }
    std::cout << __LINE__ << "\n" << _eq(Wd) << std::endl;

    {// chose a cell
      size_t k = 0;
      double r = drand48();
      std::partial_sum(pv.begin(),pv.end(),aqv.begin());
      aqv[aqv.size()-1]=1.0;
      
      for(k = 0; k < aqv.size(); k++)
	if (r < aqv[k])
	  break;

      h = all_h[k];
      for (size_t i=0; i<D; i++){
	double tmp = 0.0;
	v[i] = ((tmp=denergy_dw(i,v,h,L,J,W)) > drand48()) ? ONBIT : OFFBIT;
      }
    }// chose a cell
    //assert (false);

    // upper side probability
    Z = 0.0;
    for (size_t j=0; j<P; j++){
      ph[j] = exp(-denergy_up(j,v,all_h[j],L,J,W));
      Z += ph[j];
    }
    ph /= Z;

    // down side probability
    Z = 0.0;
    for (size_t i=0; i<D; i++){
      pv[i] = exp(-denergy_dw(i,all_v[i],h,L,J,W));
      Z += pv[i];
    }
    pv /= Z;

    // up to down
    for (size_t j=0; j<P; j++){
      Pattern &hj = all_h[j];
      for (size_t k=0; k < D; k++){
	for (size_t l=0; l < P; l++){
	  Wm(k,l) += g_transition
	    [static_cast<int>(hj[k])]
	    [static_cast<int>(hj[l])]
	    (alpha,ph[j]);
#ifdef DEBUG
// 	  std::cout << _eq(j) << _eq(k) << _eq(l) << std::endl;
// 	  std::cout << _eq(hj[k]) << _eq(hj[l]) << std::endl;	  
// 	  std::cout << _eq(g_transition
// 			   [static_cast<int>(hj[k])]
// 			   [static_cast<int>(hj[l])]
// 			   (alpha,ph[j]))
// 		    << std::endl;
// 	  plot_matrix<Layer>(gp, Wm, "Wm");
#endif
	} // l
      }	// k
      //      std::cout << _eq(hj) << std::endl;
    } // j
    //    std::cout << "Wm" << std::endl;
    //    assert(false);

    // std::cout << __LINE__ << _eq(Wm) << std::endl;
    //     std::cout << __LINE__ << _eq(v) << std::endl;
    //     std::cout << __LINE__ << _eq(h) << std::endl;
    //     assert(false);
    {// kurata's impliment is symmetry
      Wa = alpha*(Wd - Wm);
      La = alpha*(Ld - Lm);
      Ja = alpha*(Jd - Jm);

      W = W + Wa;
      L = L + La;
      J = J + Ja;
    }

    if (g_annealing_schedule_idx == 0){
      alpha = annealing_schedule(alpha);
    } else {
      alpha = annealing_schedule(t);
    }

    // output data
    if (t%10 == 0){
#ifdef DEBUG
      plot_status(gp,prefix,t,Wa);
#else
      plot_status(gp, prefix, t,
		  L, La, Ld, Lm,
		  W, Wa, Wd, Wm,
		  J, Ja, Jd, Jm);
#endif
      // plot_status(gp, prefix, t, F);
    }

    ofs_F       << t << "\t" << F       << std::endl;
    ofs_t_gibbs << t << "\t" << gibbs_P << std::endl;

    //    assert(t < 3);
  }
  pclose(gp);
  finalize_plot_table();
    
  return 0;
}
