/*

  1D Boltzman Neural Fields: Simplified Kurata (1988) model
  2 Jan 2012
  A.Date

  Input :  X[1] -- X[N_INPUTS]
  Output:  X[N_INPUTS+1] -- X[N_NEURONS]
  
  neuron 0: X[0] is always active. 
  i.e., -w[i][0] is threshoold of i-th unit.

*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>  /* for drand48 */
#include <assert.h>

#define ALPHA 0.01
#define N_INPUTS 16 
#define N_OUTPUTS 16 
#define N_NEURONS 32  /* should be  N_INPUTS + N_OUTPUTS */
#define SIGMA 0.1

#define IMAGE_SIZE 480  /* window size. arbitrary number */

int PARAMETER_k = 5;   /* 第2層目において隣接する素子を何個興奮させるか．*/

double W[N_NEURONS+1][N_NEURONS+1];  /* should be symmetry, W[i][j]=W[j][i] */
double X[N_NEURONS+1];

double F[N_NEURONS+1][N_NEURONS+1];
double G[N_NEURONS+1][N_NEURONS+1];

double q[N_OUTPUTS];

FILE *gp;
char buf[100];

//int RAND_SEED = 5666;
int RAND_SEED = 4423;
int N_LEARNING = 10000;

void init_1d_dk2011bnf();
void generate_1d_input_pattern();
void one_cycle();
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
void generate_input_from_upper_layer();
void init_1d_dk2011_bnf();

void write_data ();
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

  // set_1d_ordered_connection();
  // test_symmetric_connections(1);
  init_gnuplot ();

  for (i = 0; i < N_LEARNING; i++) {
    if (i % 10 == 0) {
      sprintf (buf, "t=%d", i);
      //plot3d_gnuplot (buf);
      plot3d_gnuplot_file(buf,i);
    }
    one_cycle();
  }
  
}


void init_1d_dk2011_bnf(){

  int i, j;
  int start, end;
  double k = (double)PARAMETER_k;

  X[0] = 0.0; /* 0th unit is always active */

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
      W[i][j] += SIGMA*nrand();
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
  int n_iterate = 100;

  int hebb;

  
  for (i = 1; i <= N_NEURONS; i++){
    X[i]=0.0;
  }

  for (i = 0; i < n_iterate; i++) {
    
    generate_1d_input_pattern();
    hebb = 1;
    learning_forward(hebb);


    printf("+");
    print_activity_pattern();


    generate_input_from_upper_layer(); /* check ! */
    hebb = -1;
    learning_forward(hebb);
    
    
    printf("-");
    print_activity_pattern();

    printf("\n");
    
    //    assert(i < 3);
  }

}




/*  第1層になんらかの入力パターンを入力してから，
　　この関数を呼び出す．

    引数は，学習，反学習のフラグ．　1 or -1　

　　この関数の仕事：

1.  入力をもとに，16通りの出力パターン，それぞれの
　　出現確率 q[i] を計算．

2.  q[i] をもとに，層間の結合係数を学習もしくは反学習　

以上．

*/
void  learning_forward(int hebb){

  int i, j, k, l;
  int start,end;
  int argmin;
  
  double e;
  double min;
  double sum;
  
  /*　出力層に出現するのは，隣あう PARAMETER_k 個が発火する N_OUTPUTS 個パターンのみと仮定．*/
  
  
  /* 出力パターン：その1 */
  for (j = N_INPUTS+1; j <= N_INPUTS+PARAMETER_k ; j++){
    X[j]=1.0;
  }
  for (j = N_INPUTS+PARAMETER_k+1; j <= N_NEURONS; j++){
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
    for (k = j; k < j+PARAMETER_k; k++){
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
    for (k = j; k < j+PARAMETER_k; k++){
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
	    W[k][l] += ALPHA*q[j];
	  }
	  else{
	    W[k][l] -= ALPHA*q[j];
	  }
	  W[l][k] = W[k][l];
	}
      }
    }
    
  } /* end of for j */


/*
  sum = 0.0;
  for (k=1; k<=N_INPUTS; k++) {
    for (l=N_INPUTS+1; l<=N_NEURONS; l++) {
      sum += W[k][l]*W[k][l];
    }
  }
  for (k=1; k<=N_INPUTS; k++) {
    for (l=N_INPUTS+1; l<=N_NEURONS; l++) {
      W[k][l] = W[k][l]/sqrt(sum);
      W[l][k] = W[k][l];
    }
  }
*/
  
  for (j = N_INPUTS+2; j <= N_NEURONS; j++){
    q[j] = q[j] + q[j-1];
        printf("%.2lf ",q[j]);
  }
   printf("\n");

}


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
  for (i = k; i < k+PARAMETER_k; i++){
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



/* 第一層目の素子間に結合は存在しない．
   第二層内の素子間の結合は学習しない．これらは数える必要なし 
　 閾値は固定（学習しない）
*/
void
learning(int n_iter){

    int i,j;

    for (i=1; i<=N_INPUTS; i++) {
      for (j=N_INPUTS+1; j<=N_NEURONS; j++) {
//	printf("%lf\n",F[i][j] - G[i][j]);
	F[i][j] = F[i][j]/(double)n_iter;
	G[i][j] = G[i][j]/(double)N_OUTPUTS;
      }
    }
    
    for (i=1; i<=N_INPUTS; i++) {
      for (j=N_INPUTS+1; j<=N_NEURONS; j++) {
            W[i][j] = W[i][j] + ALPHA*(F[i][j] - G[i][j]);
            W[j][i] = W[i][j];
      }
    }
}



void clear_counts(){

  int i, j;
  for (i=1; i<=N_INPUTS; i++) {
    for (j=N_INPUTS+1; j<=N_NEURONS; j++) {
      F[i][j]  = 0.0;
      G[i][j]  = 0.0;
    }
    F[0][i] = 0.0;
    G[0][i] = 0.0;
  }


  
}

void count_states_F(){

  int i, j;
  for (i=1; i<=N_INPUTS; i++) {
    for (j=N_INPUTS+1; j<=N_NEURONS; j++) {
      if ( X[i] > 0.5 && X[j] > 0.5 ){
	F[i][j]  += 1.0;
      }
    }
    
    if ( X[i] > 0.5 ){
      F[0][i] += 1.0;
    }
  }

}


void count_states_G()
{

  int i, j, c;

  for (i=1; i<=N_INPUTS; i++) {
    for (j=N_INPUTS+1; j<= N_NEURONS; j++) {
      if ( X[i] > 0.5 && X[j] > 0.5 ){
	G[i][j]  += 1.0;
      }
    }
    if ( X[i] > 0.5 ){
      G[0][i] += 1.0;
    }
  }

}



void
generate_1d_input_pattern()
{

  int i, p;

  p = (int)(drand48()*(double)N_INPUTS) +1;  /* choose a position */

  for (i = 1; i <= N_INPUTS; i++) {
    X[i]=0.0;
  }
  for (i = p; i < p + PARAMETER_k; i++){
    if ( i <= N_INPUTS ){
      X[i]=1.0;
    }
    else{
      X[ i % N_INPUTS ] = 1.0;
    }
  }

}

int choose_unit_randomly(int start, int end){
    
    int k;
    k = (int)( drand48()*(double)(end-start+1) );

    return start+k;

}


double update_a_unit(int k){
    
    int i;
    double u = 0.0;
    double p1;

     for (i=0; i<=N_NEURONS; i++) {
         u += W[k][i]*X[i];
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
  sprintf (buf, "gnuplot -geometry %dx%d", IMAGE_SIZE,IMAGE_SIZE);
  gp = popen(buf,"w");
  fprintf(gp, "set terminal postscript eps color \"Times\" 20\n");
//	fprintf(gp, "set terminal tgif\n");
  fprintf(gp, "set term x11\n");
//  fprintf(gp, "set cbrange [-10:10]\n");
  fprintf(gp, "set cbrange [*:*]\n");
  fprintf(gp, "set yrange [] reverse\n");
  fprintf(gp, "unset key\n");
  fprintf(gp, "set palette defined (0 \"white\", 1 \"blue\")\n"); 
}


void
plot3d_gnuplot (char *title)
{
  fprintf(gp, "plot '-' matrix with image\n");
  write_2d_data(); 
  fprintf (gp, "set title \"%s\"\n", title);
  fflush (gp);
}

void
plot3d_gnuplot_file(char *title, int t)
{
  fprintf(gp, "set output 'dir-dat/bnf1d-%d.eps'\n", t);
  fprintf(gp, "plot '-' matrix with image\n");
  write_2d_data();
  fprintf (gp, "set title \"%s\"\n", title);
  fflush (gp);
}


/* write data for gnuplot */
void write_2d_data(){

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
