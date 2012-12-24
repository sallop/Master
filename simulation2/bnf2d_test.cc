/*

  2D Boltzman Neural Fields: Simplified Kurata (1988) Model
  2 Jan 2012
  A.Date

  2次元から1次元への写像


  データ構造:
  Input :  X[1] -- X[N_INPUTS]      4x4
  Output:  X[N_INPUTS+1] -- X[N_NEURONS]

*/

#include <cstdio>
#include <cmath>
#include <cstdlib>  /* for drand48 */
#include <boost/numeric/ublas/matrix.hpp>
#include "Layer.hh"

#define N_INPUTS  64   /* should be  N_XUNITS* N_YUNITS */
#define N_OUTPUTS 32 
#define N_NEURONS 96  /* should be  N_INPUTS + N_OUTPUTS */
#define WINDOW_SIZE 480  /* window size. arbitrary number */

# define MAX_IMAGE_SIZE 10000 /* larger than (N_XUNITS + MARGIN) x (N_YUNITS + MARGIN) x N_OUTPUTS */

int N_XUNITS = 8;
int N_YUNITS = 8;
int WIDTH_INPUT = 3;

int MARGIN = 1;

int g_PARAMETER_k = 3;   /* 第2層目において隣接する素子を何個興奮させるか．*/
double g_sigma = 0.5;
double g_alpha = 0.01;


Layer		W ;
Layer		Wa;
Layer		Wd;
Layer		Wm;
Layer		W0;
Pattern		X;
Pattern		X0;
Pattern		v0;
Pattern		h0;
Layer		F;
Layer		G;
Probability	q;
Patterns	all_v;
Patterns	all_h;

void init_global_object()
{
//double W[N_NEURONS+1][N_NEURONS+1];  /* should be symmetry, W[i][j]=W[j][i] */
//double X[N_NEURONS+1];
//double F[N_NEURONS+1][N_NEURONS+1];
//double G[N_NEURONS+1][N_NEURONS+1];
//double q[N_OUTPUTS];
//double img[MAX_IMAGE_SIZE];
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
  v0 = Pattern(N_INPUTS,0);
  h0 = Pattern(N_INPUTS,0);

  all_v = Patterns(N_INPUTS ,Pattern(N_INPUTS,0));
  for (size_t p=0; p < N_INPUTS; p++){
    size_t x, y, i, j, r, c, index;
    Pattern &v = all_v[p];
    x = p%N_XUNITS;
    y = p/N_XUNITS;
    v = v0;

    for (i = 0; i < WIDTH_INPUT; i++){
      for (j = 0; j < WIDTH_INPUT; j++){
	r = x + i;
	c = y + j;
	if (r > N_XUNITS-1)
	  r = r - N_XUNITS;
	if (c > N_YUNITS-1)
	  c = c - N_YUNITS;
	//	index = c*N_XUNITS + r + 1;
	index = c*N_XUNITS + r;
	//X[index] = 1.0;
	v[index]=ONBIT;
      }	// for j=0; j < WIDTH_INPUT; j++
    } // for i=0; i< WIDTH_INPUT; i++
  } // for p=1; p <= N_INPUTS+1; p++

  all_h = Patterns(N_OUTPUTS,Pattern(N_OUTPUTS,0));
  for (size_t i=0; i<N_OUTPUTS; i++){
    Pattern &h = all_h[i];
    for (size_t k=0; k<g_PARAMETER_k; k++){
      h[(i+k)%N_OUTPUTS] = ONBIT;
    }
  }

  {
    char buf[126];
    sprintf(buf,
	    "%s-a%4.2lfk%dm%4.2lf",
	    g_prefix, g_alpha , g_PARAMETER_k, g_sigma);
    for (char *cp=buf; *cp ; cp++)
      if (*cp == '.') *cp = '_';
    strcpy(g_prefix,buf);
    //    printf("g_prefix=%s\n",g_prefix);
    //assert(false);
    sprintf(buf,"%s-state.dat",g_prefix);
    FILE *fp = fopen(buf,"w");
    for (size_t i=0; i<all_v.size(); i++){
      Pattern &v = all_v[i];
      for (size_t j=0; j<v.size(); j++){
	fprintf(fp,"%lf ", v[j]);
      }
      fprintf(fp,"\n");
    }
    fprintf(fp,"\n");
    for (size_t i=0; i<all_h.size(); i++){
      for (size_t j=0; j<h.size(); j++){
	fprintf(fp,"%lf ", h[j]);
      }
      fprintf(fp,"\n");
    }
    fclose(fp);
  }

}

FILE *gp;
char buf[100];

int RAND_SEED = 5666;
int N_LEARNING = 100000;

void init_1d_dk2011bnf();
void generate_1d_input_pattern();
void generate_2d_input_pattern();
void one_cycle();
void test_symmetric_connections(int k);
int choose_unit_randomly(int start, int end);
double update_a_unit(int k);
void set_1d_ordered_connection();
double compute_energy();
void learning_forward(int hebb);
void generate_input_from_upper_layer();
void init_1d_dk2011_bnf();

void write_multi_images(int nx, int ny);
void init_gnuplot ();
void plot3d_gnuplot (char *buf);
void plot3d_gnuplot_file(char *title, int t);
void print_activity_pattern();
double nrand();

void write_2d_data()
{
  int i,j;
  for (i=N_INPUTS+1; i<=N_NEURONS; i++){
    for (j=1; j<=N_INPUTS; j++){
      fprintf(gp, "%.8lf ", W[i][j]);
    }
    fprintf(gp, "\n");
  }
  fprintf(gp,"e\n");
  fprintf(gp,"e\n");
  fflush(gp); 
}

int
main (int argc, char *argv[])
{

  int i,j;
  long seed = RAND_SEED;

  for (i = 1; i < argc; i++) {
    switch (*(argv[i] + 1)) {
    case 'r':
      seed = atoi (argv[++i]);
      break;
    default:
      fprintf (stderr, "Usage : %s\n", argv[0]);
      fprintf (stderr, "\t-r : random-seed(%ld)\n", seed);
      exit (0);
      break;
    }
  }

  srand48 (seed);
  init_1d_dk2011_bnf();

  // set_1d_ordered_connection(); // 一つの解＋ノイズからはじめる．
  // test_symmetric_connections(1); // 対称結合になっているかのチェック．

  init_gnuplot ();

  for (i = 0; i < N_LEARNING; i++) {
    if (i % 100 == 0) {
      sprintf (buf, "t=%d", i);
      plot3d_gnuplot_file (buf, i);
    }
    one_cycle();
  }
  
}


void init_1d_dk2011_bnf(){

  int i, j;
  int start, end;
  double k = (double)g_PARAMETER_k;

  X[0] = 0.0; /* 0th unit is always active if we consider threshold */

  /* 結合係数の初期化 */

  for (i = 0; i <= N_NEURONS; i++) {
    W[i][i] = 0.0; /* すべての素子は自己結合なし */
  }

  /* 結合係数の初期化：入力層の素子間には結合がない */
  for (i = 1; i <= N_INPUTS; i++) {
    for (j =i+1; j <= N_INPUTS; j++) {
      W[i][j] = 0.0;   
      W[j][i] = 0.0;
    }
    W[i][0] = drand48() -0.5;
    //    W[i][0] = nrand();
    W[0][i] = W[i][0];
  }
  
  /* 結合係数の初期化：第一層と第二層の素子間の結合の初期値は乱数 */
  for (i = 1; i <= N_INPUTS; i++) {
    for (j = N_INPUTS+1; j <= N_NEURONS; j++) {
      // W[i][j] = nrand();
      W[i][j] = drand48() - 0.5;
      W[j][i] = W[i][j];
    }
  }

  /* 結合係数の初期化：第二層の素子間にも結合を考えない．
　　ボルツマンマシンとしては動かさない */
  for (i = N_INPUTS+1; i <= N_NEURONS; i++) {
    for (j = i+1; j <= N_NEURONS; j++) {
      W[i][j] = 0.0;   
      W[j][i] = 0.0;
    }
    W[i][0] = drand48() - 0.5;
    // W[i][0] = nrand();
    W[0][i] = W[i][0];
  }

}



void set_1d_ordered_connection(){

  int i,j;
  for (i = 1; i <= N_INPUTS; i++) {
    for (j = N_INPUTS+1; j <= N_NEURONS; j++) {
      if ( j == i + N_INPUTS){
	W[i][j] = 1.0;
      }
      else{
	W[i][j] = 0.0;
      }
      W[j][i] = W[i][j];
    }
  }
  
  /* add noise */
  for (i = 1; i <= N_INPUTS; i++) {
    for (j = N_INPUTS+1; j <= N_NEURONS; j++) {
      W[i][j] += g_sigma*nrand();
      W[j][i] = W[i][j];
    }
  }

}


double
compute_energy(){

  int i,j;
  double e = 0.0;

  for (i = 0; i <= N_NEURONS; i++) {
    for (j = i+1; j <= N_NEURONS; j++) {
      e += W[i][j]*X[i]*X[j];
    }
  }

  return -e;
}


void one_cycle(){

  int i, j;
  int hebb;
  
  generate_2d_input_pattern();
  hebb = 1;
  learning_forward(hebb);

  // print_activity_pattern();

  generate_input_from_upper_layer(); /* check ! */
  hebb = -1;
  learning_forward(hebb);

}



/*  第1層になんらかの入力パターンを入力し，この関数を呼び出す．
    引数 hebb は 1:学習， -1:反学習

    1.  与えられた入力をもとに，16通りの出力パターン，それぞれの出現確率 q[i] を計算．
    2.  q[i] をもとに，層間の結合係数を学習もしくは反学習　
*/
void  learning_forward(int hebb){

  int i, j, k, l;
  int start,end;
  int argmin;
  
  double e;
  double min;
  double sum;
  
  /*　出力層に出現するのは，隣あう g_PARAMETER_k 個が発火する N_OUTPUTS 個パターンのみと仮定．*/
  
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
      for (l=N_INPUTS+1; l<=N_NEURONS; l++) {
	if ( X[k] > 0.5 && X[l] > 0.5 ){
	  if ( hebb == 1 ){
	    W[k][l] += g_alpha*q[j];
	  }
	  else{
	    W[k][l] -= g_alpha*q[j];
	  }
	  W[l][k] = W[k][l];
	}
      }
    }

    /* 各素子のバイアス項の学習 */
    for (k=1; k<=N_NEURONS; k++) {
      if ( X[k] > 0.5 ){
	if ( hebb == 1 ){
	  W[0][k] += g_alpha*q[j];
	}
	else{
	  W[0][k] -= g_alpha*q[j];
	}
	W[k][0] = W[0][k];
      }
    }

    
  } /* end of for j */

  /* q[i] を累積密度関数に変換しておく */
  for (j = N_INPUTS+2; j <= N_NEURONS; j++){
    q[j] = q[j] + q[j-1];
    //    printf("%.2lf ",q[j]);
  }
  //  printf("\n");
  
}


/* すべての出力パターン q[i] の確率をもとに，各入力素子の活動を更新 */
void generate_input_from_upper_layer(){

  int i;
  int k;
  double r;

  r = drand48();

  i = N_INPUTS+1;
  while ( q[i] < r ){
    i++;
  }
  k = i;

  for (i = N_INPUTS+1; i <= N_NEURONS; i++){
    X[i]=0.0;
  }
  for (i = k; i < k+g_PARAMETER_k; i++){
    if ( i <= N_NEURONS ){
      X[i]=1.0;
    }
    else{
      X[N_INPUTS + (i % N_NEURONS) ]=1.0;
    }
  }
  
  for (i = 1; i <= N_INPUTS; i++){
    update_a_unit(i);
  }
  
}


void
generate_1d_input_pattern()
{

  int i, p;

  p = (int)(drand48()*(double)N_INPUTS) + 1;  /* choose a position */

  for (i = 1; i <= N_INPUTS; i++) {
    X[i]=0.0;
  }
  for (i = p; i < p + g_PARAMETER_k; i++){
    if ( i <= N_INPUTS ){
      X[i]=1.0;
    }
    else{
      X[ i % N_INPUTS ] = 1.0;
    }
  }

}



void
generate_2d_input_pattern ()
{

  int i, j, k, p;
  int x, y;
  int r,c;
  int index;

  p = (int)(drand48()*(double)N_INPUTS) +1;  /* choose a position (upper left) */
  x = (p-1) % N_XUNITS;
  y = (p-1) / N_XUNITS;
  /* x,y は左上の座標値 (0,0) .... (N_XUNITS-1,N_YUNITS-1) */

  //     index = y*N_XUNITS + x + 1;
  //     printf("i=%d,%d; (x,y) = (%d,%d) \n",p, index, x,y);

  for (i = 1; i <= N_INPUTS; i++) {
    X[i]=0.0;
  }

  for (i = 0; i < WIDTH_INPUT; i++) {
    for (j = 0; j < WIDTH_INPUT; j++) {
      r = x + i;
      c = y + j;
      if ( r > N_XUNITS-1 ){
        r = r - N_XUNITS;
      }
      if ( c > N_YUNITS-1 ){
        c = c - N_YUNITS;
      }
      index = c*N_XUNITS + r + 1;
      X[index] = 1.0;
      //      printf("%d ",index);
    }
  }
  //  printf("\n");

}



int choose_unit_randomly(int start, int end){
    
    int k;
    k = (int)( drand48()*(double)(end-start+1) );

    return start+k;

}


double update_a_unit(int k){
    
    int i;
    double u = 0.0;
    double p;

     for (i=0; i<=N_NEURONS; i++) {
         u += W[k][i]*X[i];
     }

     p = 1.0/(1.0 + exp(-u) );

     if ( drand48() < p ){ 
         X[k] = 1.0; 
     }
     else{
         X[k] = 0.0;
     }

     return p;

}



void  test_symmetric_connections(int k){
  int i, j;
  if ( k == 1 ){
    for (i = 0; i <= N_NEURONS; i++) {
      for (j = i+1; j <= N_NEURONS; j++) {
	printf("W[%d][%d] = %.5lf\n",i,j,W[i][j]);
      }
    }
  }
  else{
    for (i = 0; i <= N_NEURONS; i++) {
      for (j = i+1; j <= N_NEURONS; j++) {
	if ( W[i][j] != W[j][i] ){
	  printf("Somthing wrong! (i:%d,j:%d) w_ij=%.5lf, w_ji=%.5lf\n",i,j,W[i][j],W[j][i]);
	}
      }
    }
  }

}




void
init_gnuplot ()
{

  sprintf (buf, "gnuplot -geometry %dx%d", WINDOW_SIZE,WINDOW_SIZE);
  gp = popen(buf,"w");
  //  fprintf(gp, "set terminal postscript eps color \"Times\" 20\n");
  //  fprintf(gp, "set terminal tgif\n");
  fprintf(gp, "set term x11\n");
  fprintf(gp, "set cbrange[-10:10]\n");
  fprintf(gp, "set yrange [] reverse\n");
  fprintf(gp, "unset key\n");
  fprintf(gp, "set palette defined (0 \"white\", 1 \"blue\")\n"); 


}

void
plot3d_gnuplot (char *title)
{
  fprintf(gp, "plot '-' matrix with image\n");
  write_multi_images(4,4); 
  //  write_multi_images(2,2); 
  fprintf (gp, "set title \"%s\"\n", title);
  fflush (gp);
}

void
plot3d_gnuplot_file(char *title, int t)
{
  fprintf(gp, "set output 'dir-dat/bnf2d-%d.eps'\n", t);
  fprintf(gp, "plot '-' matrix with image\n");
  //  write_multi_images(4,4); 
  //  write_multi_images(2,2);
  write_2d_data();
  fprintf (gp, "set title \"%s\"\n", title);
  fflush (gp);
}


void write_multi_images(int nx, int ny)
{
  int i,j;
  int n_pixels = (N_XUNITS + MARGIN)*(N_YUNITS + MARGIN)*nx*ny;
  int width = (N_XUNITS + MARGIN)*nx;
  int index;
  int r,c;
  int a;
  int start_x, start_y;

  for (i=0; i<n_pixels; i++){
    img[i]=0.0;
  }

  for (i = N_INPUTS+1; i <= N_NEURONS; i++){
    a = (i - (N_INPUTS+1) ); 
    start_x = (a % nx)*(N_XUNITS + MARGIN);
    start_y = (a / nx)*(N_YUNITS + MARGIN);
    for (j=1; j <= N_INPUTS; j++){
      r = start_x + (j-1) % N_XUNITS;
      c = start_y + (j-1) / N_XUNITS;
      a = r + c*width;
      img[a] = W[i][j];
    }
  }

  
  for (i=0; i < n_pixels; i++){
    fprintf(gp, "%.8lf ", img[i]);
    if (  (i+1) % width  ==  0 ){
      fprintf(gp, "\n");
    }
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
