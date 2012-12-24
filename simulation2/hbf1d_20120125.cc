/*

  1D Boltzman Neural Fields: Simplified Kurata (1988) model
  2 Jan 2012
  A.Date

  Input :  X[1] -- X[N_INPUTS]
  Output:  X[N_INPUTS+1] -- X[N_NEURONS]
  
  neuron 0: X[0] is always active. 
  i.e., -w[i][0] is threshoold of i-th unit.

*/

#include <cstdio>
#include <cstring>
#include <cmath>
#include <cstdlib>  /* for drand48 */
#include <cassert>
#include <unistd.h>
#include "Layer.hh"

#define N_INPUTS   16 
#define N_OUTPUTS  16 
#define N_NEURONS  32   /* should be  N_INPUTS + N_OUTPUTS */
#define IMAGE_SIZE 480  /* window size. arbitrary number */
int RAND_SEED = 4423;
int N_LEARNING = 10000;



double g_alpha = 0.01;
double g_sigma = 0.1;
char g_prefix[32]="d1f";
int g_PARAMETER_k = 5;   /* 第2層目において隣接する素子を何個興奮させるか．*/
int g_cmds_idx = 1;
int g_terminal_idx = 1;	     // 0:eps 1:gif 2:x11
/* should be symmetry, W[i][j]=W[j][i] */


//double W[N_NEURONS+1][N_NEURONS+1];
//double X[N_NEURONS+1];
//double F[N_NEURONS+1][N_NEURONS+1];
//double G[N_NEURONS+1][N_NEURONS+1];

//double q[N_OUTPUTS]; wrong!
//double q[N_NEURONS+1]; right

Layer W ;
Layer Wa;
Layer Wd;
Layer Wm;
Layer W0;
Pattern X;
Pattern X0;
Layer F;
Layer G;
Probability q;
Patterns all_v;
Patterns all_h;

typedef double (*func_ptr)(double);
double transition_off(double prob){ return 0.0; }
double transition_on (double prob){ return g_alpha*prob; }
func_ptr transition_tbl[2][2]={
  {transition_off, transition_off},
  {transition_off, transition_on },
};

enum TerminalType{EPS=0,GIF=1,X11=2};

void init_global_object()
{
  W0= zerom(N_NEURONS+1,N_NEURONS+1);
  W = Layer(N_NEURONS+1,N_NEURONS+1);
  Wa= Layer(N_NEURONS+1,N_NEURONS+1);
  Wd= Layer(N_NEURONS+1,N_NEURONS+1);
  Wm= Layer(N_NEURONS+1,N_NEURONS+1);

  F = Layer(N_NEURONS+1,N_NEURONS+1);
  G = Layer(N_NEURONS+1,N_NEURONS+1);

  X = Pattern(N_NEURONS+1,0);
  X0= Pattern(N_NEURONS+1,0);
  q = Probability(N_NEURONS+1);

  all_v = Patterns(N_INPUTS ,Pattern(N_INPUTS,0));
  for (size_t i=0; i<N_INPUTS; i++){
    Pattern &v = all_v[i];
    for (size_t k=0; k<g_PARAMETER_k; k++){
      v[(i+k)%N_INPUTS] = ONBIT;
    }
    std::cout << v << std::endl;
  }
  std::cout << "------------------------------------" << std::endl;
  all_h = Patterns(N_OUTPUTS,Pattern(N_OUTPUTS,0));
  for (size_t i=0; i<N_OUTPUTS; i++){
    Pattern &h = all_h[i];
    for (size_t k=0; k<g_PARAMETER_k; k++){
      h[(i+k)%N_OUTPUTS] = ONBIT;
    }
    std::cout << h << std::endl;
  }

  char buf[126]="";
  sprintf(buf,
	  "%s-"		// prefix
	  "a%4.2lf"		// g_alpha
	  "k%d"		// g_PARAMETER_k
	  "m%4.2lf"		// g_sigma
	  "s%d"		// RAND_SEED
	  "r%d",		// g_cmds_idx
	  g_prefix,
	  g_alpha ,
	  g_PARAMETER_k,
	  g_sigma,
	  RAND_SEED,
	  g_cmds_idx
	  );
  for (char *cp=buf; *cp ; cp++)
    if (*cp == '.') *cp = '_';
  strcpy(g_prefix,buf);
  //    printf("g_prefix=%s\n",g_prefix);
  //assert(false);
  std::cout << "prefix =" << g_prefix << std::endl;
  //    assert(false);
  sprintf(buf,"dir-dat/%s-state.dat",g_prefix);
}

FILE *gp;
char buf[100];

//int RAND_SEED = 5666;

void init_1d_dk2011bnf();
void generate_1d_input_pattern(SubPattern& v, SubPattern& h);
void one_cycle(SubPattern& v, SubPattern& h);
void test_symmetric_connections(int k);
void clear_counts();
int choose_unit_randomly(int start, int end);
double update_a_unit(int k);
void count_states_F();
void count_states_G();
void learning(int n_iter);
void set_1d_ordered_connection();
double compute_energy();

void learning_forward(int hebb);
//void learning_forward_hebb();
//void learning_forward_antihebb();
void learning_forward_hebb(SubPattern& v, SubPattern& h);
void learning_forward_antihebb(SubPattern& v, SubPattern& h);
void generate_input_from_upper_layer(SubPattern& v, SubPattern& h);
void init_1d_dk2011_bnf();

void write_2d_data(Layer &M);
void init_gnuplot ();
void plot3d_gnuplot (char *title);
void plot3d_gnuplot_file(char *title, int t);
void write_2d_data();
void print_activity_pattern();
double nrand();

int
main (int argc, char *argv[])
{

  int i;
  int result;
  while ((result = getopt(argc,argv,"a:k:m:s:r:t:")) != -1){
    switch (result){
    case 'a':
      // alpha: learning rate
      g_alpha = atof(optarg);
      break;
    case 'm':
      // sigma: 
      g_sigma = atof(optarg);
      break;
    case 'k':
      g_PARAMETER_k = atoi(optarg);
      break;
    case 's':
      // seed: random seed
      RAND_SEED = atof(optarg);
      break;
    case 'r':
      g_cmds_idx = atoi(optarg);
      break;
    case 't':
      g_terminal_idx = atoi(optarg);
      break;
    default:
      std::cout << optarg << " is not implimented." << std::endl;
      assert (false);
    }
  }

  srand48(RAND_SEED);

  init_global_object ();

  init_1d_dk2011_bnf();

  // set_1d_ordered_connection();
  // test_symmetric_connections(1);
  init_gnuplot();

  SubPattern v(X, ub::range(1         , N_INPUTS +1));
  SubPattern h(X, ub::range(N_INPUTS+1, N_NEURONS+1));
  for (i = 0; i < N_LEARNING; i++) {
    if (i % 10 == 0) {
      sprintf (buf, "t=%d", i);
      plot3d_gnuplot_file(buf,i);
    }
    one_cycle(v,h);
  }
}

void init_1d_dk2011_bnf()
{
  int i, j;
  int start, end;
  double k = (double)g_PARAMETER_k;

  X[0] = 0.0; /* 0th unit is always active */

  /* 結合係数の初期化 */

  for (i = 0; i <= N_NEURONS; i++)
    W(i,i) = 0.0; /* すべての素子は自己結合なし */

  /* 結合係数の初期化：入力層の素子間には結合がない */
  for (i = 1; i <= N_INPUTS; i++){
    for (j =i+1; j <= N_INPUTS; j++) {
      W(i,j) = 0.0;   
      W(j,i) = 0.0;
    }
    W(i,0) = drand48() -0.5;
    //    W(i,0) = nrand();
    W(0,i) = W(i,0);
  }
  
  /* 結合係数の初期化：第一層と第二層の素子間の結合の初期値は乱数 */
  for (i = 1; i <= N_INPUTS; i++) {
    for (j = N_INPUTS+1; j <= N_NEURONS; j++) {
      // W(i,j) = nrand();
      W(i,j) = g_sigma*(drand48() - 0.5);
      W(j,i) = W(i,j);
    }
  }

  /* 結合係数の初期化：第二層の素子間にも結合を考えない．
　　ボルツマンマシンとしては動かさない */
  for (i = N_INPUTS+1; i <= N_NEURONS; i++) {
    for (j = i+1; j <= N_NEURONS; j++) {
      W(i,j) = 0.0;   
      W(j,i) = 0.0;
    }
    W(i,0) = g_sigma*(drand48() - 0.5);
    // W(i,0) = nrand();
    W(0,i) = W(i,0);
  }
}

void set_1d_ordered_connection()
{
  int i,j;
  for (i = 1; i <= N_INPUTS; i++) {
    for (j = N_INPUTS+1; j <= N_NEURONS; j++) {
      if ( j == i + N_INPUTS ){ 
	W(i,j) = 1.0;
      }
      else{
	W(i,j) = 0.0;
      }
      W(j,i) = W(i,j);
    }
  }
  
  /* add noise */
  for (i = 1; i <= N_INPUTS; i++) {
    for (j = N_INPUTS+1; j <= N_NEURONS; j++) {
      W(i,j) += g_sigma*nrand();
      W(j,i) = W(i,j);
    }
  }
}

double
compute_energy()
{
  int i,j;
  double e = 0.0;
  for (i = 0; i <= N_NEURONS; i++) {
    for (j = i+1; j <= N_NEURONS; j++) {
      e += W(i,j)*X[i]*X[j];
    }
  }
  return -e;
}

void one_cycle(SubPattern& v, SubPattern& h)
{
  int i, j;

  static const int n_iterate = 100;
  int hebb;

  X=X0;
  Wd = Wm = W0;
  for (i = 0; i < n_iterate; i++){
    generate_1d_input_pattern(v,h);
    learning_forward_hebb(v,h);
    generate_input_from_upper_layer(v,h); /* check ! */
    learning_forward_antihebb(v,h);
  }
  Wd /= static_cast<double>(n_iterate);
  Wm /= static_cast<double>(n_iterate);
  Wa = Wd - Wm;
}

void  learning_forward(int hebb)
{
  int i, j, k, l;
  int start,end;
  int argmin;
  double e;
  double min;
  double sum;
  /*
   *　出力層に出現するのは, 隣あう PARAMETER_k 個が発火する 
   * N_OUTPUTS 個パターンのみと仮定．
   */
    
  /* 出力パターン：その1 */
  for (j = N_INPUTS+1; j <= N_INPUTS+g_PARAMETER_k ; j++){
    X[j]=1.0;
  }
  for (j = N_INPUTS+g_PARAMETER_k+1; j <= N_NEURONS; j++){
    X[j]=0.0;
  }
  argmin = N_INPUTS+1;
  min = compute_energy();
  q[N_INPUTS+1] = min; 
  
  /* 出力パターン：その2以降 */
  for (j = N_INPUTS+2; j <= N_NEURONS; j++) {
    for (k = N_INPUTS+1; k <= N_NEURONS; k++){
      X[k]=0.0;
    }
    for (k = j; k < j+g_PARAMETER_k; k++){
      if ( k <= N_NEURONS ){
	X[k]=1.0;
      }
      else{
	X[N_INPUTS + (k % N_NEURONS) ]=1.0;
      }
    }
    e = compute_energy();
    q[j] = e;
    if ( e < min ){
      min = e;
      argmin = j;
    }
  }
  
  for (j = N_INPUTS+1; j <= N_NEURONS; j++){
    q[j] = exp(-q[j]);
  }
  
  sum = 0.0;
  for (j = N_INPUTS+1; j <= N_NEURONS; j++){
    sum += q[j];
  }
  //    printf("sum=%.2lf\n",sum);
  for (j = N_INPUTS+1; j <= N_NEURONS; j++){
    q[j] = q[j]/sum;
    //     printf("%.2lf ",q[j]);
  }
  //  printf("\n");
  
  /* 学習 */
  for (j = N_INPUTS+1; j <= N_NEURONS; j++) {
    /* 出力層にパターンをひとつづつあてはめてみる */
    for (k = N_INPUTS+1; k <= N_NEURONS; k++){
      X[k]=0.0;
    }
    for (k = j; k < j+g_PARAMETER_k; k++){
      if ( k <= N_NEURONS ){
	X[k]=1.0;
      }
      else{
	X[N_INPUTS + (k % N_NEURONS) ]=1.0;
      }
    }
    /* k番目の素子（入力層）と l番目の素子の学習 */
    for (k=1; k<=N_INPUTS; k++) {
      for (l=N_INPUTS+1; l<=N_NEURONS; l++){
	double weight = transition_tbl
	  [static_cast<int>(X[k])]
	  [static_cast<int>(X[l])]
	  (q[j]);
	if ( hebb == 1 ){
	  Wd(k,l) += weight;
 	  W (k,l) += weight;
	  Wd(l,k) = Wd(k,l);
	}
	else{
	  Wm(k,l) += weight;
 	  W (k,l) -= weight;
	  Wm(l,k) = Wm(k,l);	  
	}
	W(l,k) = W(k,l);
      }
    }    
  } /* end of for j */

  for (j = N_INPUTS+2; j <= N_NEURONS; j++){
    q[j] = q[j] + q[j-1];
    printf("%.2lf ",q[j]);
  }
  printf("\n");
}


void learning_forward_hebb(SubPattern &v, SubPattern &h)
{
  int i, j, k, l;
  int start,end;
  int argmin;
  double e;
  double min;
  double sum;
  /*
   *　出力層に出現するのは, 隣あう PARAMETER_k 個が発火する 
   * N_OUTPUTS 個パターンのみと仮定．
   */
    
  /* 出力パターン：その1 */
  //or (j = N_INPUTS+1; j <= N_INPUTS+g_PARAMETER_k ; j++){
  // X[j]=1.0;
  //
  //or (j = N_INPUTS+g_PARAMETER_k+1; j <= N_NEURONS; j++){
  // X[j]=0.0;
  //
  h = all_h[0];

  argmin = N_INPUTS+1;
  min = compute_energy();
  q[N_INPUTS+1] = min; 
  
  /* 出力パターン：その2以降 */
  for (j = N_INPUTS+2; j <= N_NEURONS; j++) {
    // for (k = N_INPUTS+1; k <= N_NEURONS; k++){
    //       X[k]=0.0;
    //     }
    //     for (k = j; k < j+g_PARAMETER_k; k++){
    //       if ( k <= N_NEURONS ){
    // 	X[k]=1.0;
    //       }
    //       else{
    // 	X[N_INPUTS + (k % N_NEURONS) ]=1.0;
    //       }
    //     }
    h = all_h[j-(N_INPUTS+2)];

    e = compute_energy();
    q[j] = e;
    if ( e < min ){
      min = e;
      argmin = j;
    }
  }
  
  for (j = N_INPUTS+1; j <= N_NEURONS; j++){
    q[j] = exp(-q[j]);
  }
  
  sum = 0.0;
  for (j = N_INPUTS+1; j <= N_NEURONS; j++){
    sum += q[j];
  }
  //    printf("sum=%.2lf\n",sum);
  for (j = N_INPUTS+1; j <= N_NEURONS; j++){
    q[j] = q[j]/sum;
    //     printf("%.2lf ",q[j]);
  }
  //  printf("\n");
  
  /* 学習 */
  for (j = N_INPUTS+1; j <= N_NEURONS; j++){
    /* 出力層にパターンをひとつづつあてはめてみる */
    // for (k = N_INPUTS+1; k <= N_NEURONS; k++){
    //       X[k]=0.0;
    //     }
    //     for (k = j; k < j+g_PARAMETER_k; k++){
    //       if ( k <= N_NEURONS ){
    // 	X[k]=1.0;
    //       }
    //       else{
    // 	X[N_INPUTS + (k % N_NEURONS) ]=1.0;
    //       }
    //     }
    h = all_h[j-(N_INPUTS+1)];

    /* k番目の素子（入力層）と l番目の素子の学習 */
    for (k=1; k<=N_INPUTS; k++) {
      for (l=N_INPUTS+1; l<=N_NEURONS; l++){
	double weight = transition_tbl
	  [static_cast<int>(X[k])]
	  [static_cast<int>(X[l])]
	  (q[j]);
	Wd(k,l) += weight;
	W (k,l) += weight;
	Wd(l,k) = Wd(k,l);
	W (l,k) = W(k,l);
      }
    }
  } /* end of for j */

  for (j = N_INPUTS+2; j <= N_NEURONS; j++){
    q[j] = q[j] + q[j-1];
    printf("%.2lf ",q[j]);
  }
  printf("\n");
}

void learning_forward_antihebb(SubPattern& v, SubPattern& h)
{
  int i, j, k, l;
  int start,end;
  int argmin;
  double e;
  double min;
  double sum;
  /*
   *　出力層に出現するのは, 隣あう PARAMETER_k 個が発火する 
   * N_OUTPUTS 個パターンのみと仮定．
   */

  /* 出力パターン：その1 */
  h = all_h[0];
  argmin = N_INPUTS+1;
  min = compute_energy();
  q[N_INPUTS+1] = min; 
  
  /* 出力パターン：その2以降 */
  for (j = N_INPUTS+2; j <= N_NEURONS; j++){
    h = all_h[j-(N_INPUTS+2)];
    e = compute_energy();
    q[j] = e;
    if ( e < min ){
      min = e;
      argmin = j;
    }
  }
  
  for (j = N_INPUTS+1; j <= N_NEURONS; j++){
    q[j] = exp(-q[j]);
  }
  
  sum = 0.0;
  for (j = N_INPUTS+1; j <= N_NEURONS; j++){
    sum += q[j];
  }
  for (j = N_INPUTS+1; j <= N_NEURONS; j++){
    q[j] = q[j]/sum;
  }
  /* 学習 */
  for (j = N_INPUTS+1; j <= N_NEURONS; j++) {
    /* 出力層にパターンをひとつづつあてはめてみる */
    h = all_h[j-(N_INPUTS+1)];
    /* k番目の素子（入力層）と l番目の素子の学習 */
    for (k=1; k<=N_INPUTS; k++) {
      for (l=N_INPUTS+1; l<=N_NEURONS; l++){
	double weight = transition_tbl
	  [static_cast<int>(X[k])]
	  [static_cast<int>(X[l])]
	  (q[j]);
	Wm(k,l) += weight;
	W (k,l) -= weight;
	Wm(l,k) = Wm(k,l);	  
	W(l,k) = W(k,l);
      }
    }
  } /* end of for j */
  for (j = N_INPUTS+2; j <= N_NEURONS; j++){
    q[j] = q[j] + q[j-1];
    printf("%.2lf ",q[j]);
  }
  printf("\n");
}

void generate_input_from_upper_layer(SubPattern &v, SubPattern &h)
{
  int i;
  int k;
  double r;
  r = drand48();
  i = N_INPUTS+1;
  while ( q[i] < r ){
    i++;
  }
  k = i;
  h = all_h[k-(N_INPUTS+1)];
  for (i = 1; i <= N_INPUTS; i++){
    update_a_unit(i);
  }
}

/* 第一層目の素子間に結合は存在しない．
   第二層内の素子間の結合は学習しない．これらは数える必要なし 
　 閾値は固定（学習しない）
*/
void
learning(int n_iter)
{
  int i,j;

  for (i=1; i<=N_INPUTS; i++) {
    for (j=N_INPUTS+1; j<=N_NEURONS; j++) {
      //	printf("%lf\n",F[i][j] - G[i][j]);
      F(i,j) = F(i,j)/(double)n_iter;
      G(i,j) = G(i,j)/(double)N_OUTPUTS;
    }
  }
    
  for (i=1; i<=N_INPUTS; i++) {
    for (j=N_INPUTS+1; j<=N_NEURONS; j++) {
      W(i,j) = W(i,j) + g_alpha*(F(i,j) - G(i,j));
      W(j,i) = W(i,j);
    }
  }
}

void clear_counts()
{
  int i, j;
  for (i=1; i<=N_INPUTS; i++) {
    for (j=N_INPUTS+1; j<=N_NEURONS; j++) {
      F(i,j)  = 0.0;
      G(i,j)  = 0.0;
    }
    F(0,i) = 0.0;
    G(0,i) = 0.0;
  }
}

void count_states_F()
{
  int i, j;
  for (i=1; i<=N_INPUTS; i++) {
    for (j=N_INPUTS+1; j<=N_NEURONS; j++) {
      if ( X[i] > 0.5 && X[j] > 0.5 ){
	F(i,j)  += 1.0;
      }
    }
    if ( X[i] > 0.5 ){
      F(0,i) += 1.0;
    }
  }
}

void count_states_G()
{
  int i, j, c;
  for (i=1; i<=N_INPUTS; i++) {
    for (j=N_INPUTS+1; j<= N_NEURONS; j++) {
      if ( X[i] > 0.5 && X[j] > 0.5 ){
	G(i,j)  += 1.0;
      }
    }
    if ( X[i] > 0.5 ){
      G(0,i) += 1.0;
    }
  }
}

void
generate_1d_input_pattern(SubPattern& v, SubPattern& h)
{
  int i, p;
  p = (int)(drand48()*(double)N_INPUTS) +1;/* choose a position */
  v = all_v[p-1];
}

int
choose_unit_randomly(int start, int end)
{
  int k;
  k = (int)( drand48()*(double)(end-start+1) );
  return start+k;
}

double update_a_unit(int k)
{
  int i;
  double u = 0.0;
  double p1;

  for (i=0; i<=N_NEURONS; i++) {
    u += W(k,i)*X[i];
  }

  p1 = 1.0/(1.0 + exp(-u) );

  if ( drand48() < p1 ){ 
    X[k] = 1.0; 
  }
  else{
    X[k] = 0.0;
  }
  return p1;
}

void
test_symmetric_connections(int k)
{
  int i, j;
  if ( k == 1 ){
    for (i = 0; i <= N_NEURONS; i++) {
      for (j = i+1; j <= N_NEURONS; j++) {
	printf("W(%d,%d) = %.5lf\n",i,j,W(i,j));
      }
    }
  }
  else{
    for (i = 0; i <= N_NEURONS; i++) {
      for (j = i+1; j <= N_NEURONS; j++) {
	if ( W(i,j) != W(j,i) ){
	  printf("Somthing wrong! (i:%d,j:%d) w_ij=%.5lf, w_ji=%.5lf\n"
		 ,i,j,W(i,j),W(j,i));
	}
      }
    }
  }
}

void
init_gnuplot ()
{
  sprintf (buf, "gnuplot -geometry %dx%d", IMAGE_SIZE,IMAGE_SIZE);
  gp = popen(buf,"w");
  //  fprintf(gp, "set terminal postscript enhanced color eps\n");
  fprintf(gp, "set terminal gif\n");
  fprintf(gp, "set cbrange [*:*]\n");
  fprintf(gp, "set yrange [] reverse\n");
  fprintf(gp, "set palette defined (0 \"white\", 1 \"blue\")\n"); 
}

void
plot3d_gnuplot_file(char *title, int t)
{
  double n_matrix = 4.0;
  double x_offset = 0.0;
#define y_offset(n) ((n_matrix - (n))*(1.0/n_matrix))
  fprintf(gp, "set output 'dir-dat/%s-t%d.gif\n"
	  , g_prefix
	  , t);
  //, g_terminal_suffix[g_terminal_idx]);

  fprintf(gp, "set rmargin %lf\n", 2.5);
  fprintf(gp, "set lmargin %lf\n", 1.8);
  fprintf(gp, "unset key\n");
  fprintf(gp, "unset xtics\n");
  fprintf(gp, "unset ytics\n");
  
  fprintf(gp, "set size square\n");
  fprintf(gp, "set multiplot layout %d,%1.0lf\n",1,n_matrix);

  fprintf(gp, "set cbrange[-10:10]\n");
  fprintf(gp, "plot '-' matrix with image title 'W'\n");
  write_2d_data(W);
  // --------------------------------------------------
  //  fprintf(gp, "set origin %lf, %lf\n", x_offset, y_offset(2));
  fprintf(gp, "set cbrange[-0.01:0.01]\n");
  fprintf(gp, "plot '-' matrix with image title 'Wa'\n");
  write_2d_data(Wa);
  // --------------------------------------------------
  //  fprintf(gp, "set origin %lf, %lf\n", x_offset, y_offset(3));
  fprintf(gp, "set cbrange[-0.01:0.01]\n");
  fprintf(gp, "plot '-' matrix with image title 'Wd'\n");
  write_2d_data(Wd);
  // --------------------------------------------------
  //  fprintf(gp, "set origin %lf, %lf\n", x_offset, y_offset(4));
  fprintf(gp, "set cbrange[-0.01:0.01]\n");
  fprintf(gp, "plot '-' matrix with image title 'Wm'\n");
  write_2d_data(Wm);
  // --------------------------------------------------
  fprintf(gp, "unset multiplot\n");
  fflush (gp);
#undef y_offset
}


/* write data for gnuplot */
void write_2d_data(Layer &M)
{
  int i,j;
  for (i=N_INPUTS+1; i<=N_NEURONS; i++){
    for (j=1; j<=N_INPUTS; j++){
      fprintf(gp, "%.8lf ", M(i,j));
    }
    fprintf(gp, "\n");
  }
  fprintf(gp,"e\n");
  fprintf(gp,"e\n");
  fflush(gp); 
}

void print_activity_pattern(){

  int i;

  for (i=1; i<=N_INPUTS; i++){
    printf("%.0lf",X[i]);
  }
  printf("\n ");
  for (i=N_INPUTS+1; i<=N_NEURONS; i++){
    printf("%.0lf", X[i]);
  }
  printf("\n");
}

/* 
  nrand() :  a random sample from standard normal distribution  
  一様分布 drand48() を使い標準正規分布に従うデータを出力する関数 
*/
double nrand()
{
  static int sw=0;
  static double r1,r2,s;
  
  if (sw==0){
    sw=1;
    do {
      r1=2.0*drand48()-1.0;
      r2=2.0*drand48()-1.0;
      s=r1*r1+r2*r2;
    } while (s>1.0 || s==0.0);
    s=sqrt(-2.0*log(s)/s);
    return(r1*s);
  }
  else {
    sw=0;
    return(r2*s);
  }
}
