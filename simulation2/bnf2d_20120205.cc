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
#include <cstring>
#include <cmath>
#include <cstdlib>  /* for drand48 */
#include <numeric>
#include <boost/numeric/ublas/matrix.hpp>
#include "Layer.hh"

void	generate_2d_input_pattern(SubPattern& v, SubPattern &h);
void	one_cycle(SubPattern &v, SubPattern &h);
void	test_symmetric_connections(int k);
int	choose_unit_randomly(int start, int end);
double	update_a_unit(int k);
void	set_1d_ordered_connection();
double	compute_energy();
void	learning_forward(int hebb);
void	learning_forward_hebb(SubPattern &v, SubPattern &h);
void	learning_forward_antihebb(SubPattern &v, SubPattern &h);
void	generate_input_from_upper_layer(SubPattern &v, SubPattern &h);
void	init_1d_dk2011_bnf(SubPattern &v, SubPattern &h);
void	write_2d_data(Layer &M);

void	write_multi_images(int nx, int ny, Layer& W);
void	init_gnuplot ();
void	plot3d_gnuplot (char *buf);
void	plot3d_gnuplot_file(char *title, int t);
void	print_activity_pattern();
double	nrand();

//#define N_INPUTS	64	/* should be  N_XUNITS* N_YUNITS */
//#define N_OUTPUTS	32
//#define N_NEURONS	96	/* should be  N_INPUTS + N_OUTPUTS */
#define N_INPUTS	64	/* should be  N_XUNITS* N_YUNITS */
#define N_OUTPUTS	16
#define N_NEURONS	80	/* should be  N_INPUTS + N_OUTPUTS */
#define WINDOW_SIZE	480	/* window size. arbitrary number */
#define MAX_IMAGE_SIZE	10000

int	N_XUNITS      = 8;
int	N_YUNITS      = 8;

int	WIDTH_INPUT   = 3;
int	MARGIN	      = 1;
int	g_PARAMETER_k = 3;/* 第2層目において隣接する素子を何個興奮させるか．*/
double	g_sigma	      = 1.0;
double	g_alpha	      = 0.01;
int     g_cmds_idx    = 0;
int     g_terminal_idx= 0;
double	img[MAX_IMAGE_SIZE];
int	RAND_SEED  = 5666 ;
int	N_LEARNING = 10000;
//int	N_LEARNING = 100;

int g_counter = 0;

FILE *gp;
char buf[100];
//char g_prefix[32]="d2";
char g_prefix[32]="d2f";
//char g_prefix[32]="d2v";
Layer		W ;
Layer		Wa;
Layer		Wd;
Layer		Wm;
Layer		W0;
Pattern		X ;
Pattern		X0;
Pattern		v0;
Pattern		h0;
Layer		F;
Layer		G;
Probability	q;
Patterns	all_v;
Patterns	all_h;

typedef double (*func_ptr)(double);
double transition_off(double prob){ return 0.0; }
double transition_on (double prob){ return g_alpha*prob; }
func_ptr transition_tbl[2][2]={
  {transition_off, transition_off},
  {transition_off, transition_on },
};

enum TerminalType{EPS=0,GIF=1,X11=2};

void write_2d_data(Layer &W)
{
  int i,j;
  for (i=N_INPUTS+1; i<=N_NEURONS; i++){
    for (j=1; j<=N_INPUTS; j++){
      fprintf(gp, "%.8lf ", W(i,j));
    }
    fprintf(gp, "\n");
  }
  fprintf(gp,"e\n");
  fprintf(gp,"e\n");
  fflush(gp); 
}
    
void write_2d_all(Layer &W)
{
  int i,j;
  for (i=0; i<=N_NEURONS; i++){
    for (j=0; j<=N_NEURONS; j++){
      fprintf(gp, "%.8lf ", W(i,j));
    }
    fprintf(gp, "\n");
  }
  fprintf(gp,"e\n");
  fprintf(gp,"e\n");
  fflush(gp); 
}

void init_global_object()
{
  W0 = zerom(N_NEURONS+1,N_NEURONS+1);
  W  = Layer(N_NEURONS+1,N_NEURONS+1);
  Wa = Layer(N_NEURONS+1,N_NEURONS+1);
  Wd = Layer(N_NEURONS+1,N_NEURONS+1);
  Wm = Layer(N_NEURONS+1,N_NEURONS+1);
  F  = Layer(N_NEURONS+1,N_NEURONS+1);
  G  = Layer(N_NEURONS+1,N_NEURONS+1);
  X  = Pattern(N_NEURONS+1,0);
  X0 = Pattern(N_NEURONS+1,0);
  q  = Probability(N_NEURONS+1,0.0);
  v0 = Pattern(N_INPUTS,0);
  h0 = Pattern(N_INPUTS,0);

  all_v = Patterns(N_INPUTS ,Pattern(N_INPUTS,0));
  for (size_t p=0; p < N_INPUTS; p++){
    size_t x, y, i, j, r, c, index;
    Pattern &v = all_v[p];
    x = p%N_XUNITS;
    y = p/N_XUNITS;
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

  {//
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
    FILE *fp = fopen(buf,"w");
    Pattern& h=all_h[0];
    for (size_t i=0; i<h.size(); i++){
      fprintf(fp,"%lf ", h[i]);
      if((i+1)%N_XUNITS==0)
	fprintf(fp,"\n");
    }
    fclose(fp);
  }
  //  assert (false);
}

int
main (int argc, char *argv[])
{
  int i,j;
  int result;
  while ((result = getopt(argc,argv,"a:k:m:s:p:r:t:c:")) != -1){
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
    case 't':
      // set terminal
      g_terminal_idx = atoi(optarg);
      break;
    case 'p':
      strcpy(g_prefix,optarg);
      break;
    case 'r':
      // rayout - plot commands
      g_cmds_idx = atoi(optarg);
      break;
    default:
      std::cout << optarg << " is not implimented." << std::endl;
      assert (false);
    }
  }

  init_global_object();
  
  SubPattern v(X, ub::range(1         , 1 + N_INPUTS));
  SubPattern h(X, ub::range(1+N_INPUTS, N_NEURONS+1 ));
  srand48(RAND_SEED);
  init_1d_dk2011_bnf(v,h);

  // set_1d_ordered_connection(); // 一つの解＋ノイズからはじめる．
  // test_symmetric_connections(1); // 対称結合になっているかのチェック．

  init_gnuplot ();
  for (i = 0; i <= N_LEARNING; i++){
  //  for (i = 0; i < 500; i++){
    Wm = W0;
    Wd = W0;

    one_cycle(v,h);
    Wa = Wd - Wm;
    if (i % 10 == 0) {
      sprintf (buf, "t=%d", i);
      plot3d_gnuplot_file (buf, i);
    }
  }
}

void init_1d_dk2011_bnf(SubPattern &v, SubPattern &h)
{
  int i, j;
  int start, end;
  double k = (double)g_PARAMETER_k;

  X[0] = 0.0; /* 0th unit is always active if we consider threshold */
  /* 結合係数の初期化 */
  for (i = 0; i <= N_NEURONS; i++) {
    W(i,i) = 0.0; /* すべての素子は自己結合なし */
  }

  /* 結合係数の初期化：入力層の素子間には結合がない */
  for (i = 1; i <= N_INPUTS; i++){
    for (j =i+1; j <= N_INPUTS; j++){
      W(i,j) = 0.0;   
      W(j,i) = 0.0;
    }
    W(i,0) = g_sigma*(drand48() - 0.5);
    //    W(i,0) = nrand();
    W(0,i) = W(i,0);
  }
  
  /* 結合係数の初期化：第一層と第二層の素子間の結合の初期値は乱数 */
  for (i = 1; i <= N_INPUTS; i++) {
    for (j = N_INPUTS+1; j <= N_NEURONS; j++) {
      // W[i][j] = nrand();
      W(i,j) = g_sigma*(drand48() - 0.5);
      W(j,i) = W(i,j);
    }
  }

  /*
   * 結合係数の初期化：第二層の素子間にも結合を考えない．
   * ボルツマンマシンとしては動かさない 
   */
  for (i = N_INPUTS+1; i <= N_NEURONS; i++) {
    for (j = i+1; j <= N_NEURONS; j++) {
      W(i,j) = 0.0;   
      W(j,i) = 0.0;
    }
    W(i,0) = g_sigma*(drand48() - 0.5);
    W(0,i) = W(i,0);
  }
}

double
compute_energy()
{
  int i,j;
  double e = 0.0;
  for (i = 0; i <= N_NEURONS; i++)
    for (j = i+1; j <= N_NEURONS; j++)
      e += W(i,j)*X[i]*X[j];
  return -e;
}

void one_cycle(SubPattern &v, SubPattern &h)
{
  int i,j;
  int hebb;

  generate_2d_input_pattern(v,h);
  learning_forward_hebb(v,h);
  generate_input_from_upper_layer(v,h); /* check ! */
  learning_forward_antihebb(v,h);
}

void
learning_forward_hebb(SubPattern& v, SubPattern& h)
{
  int i, j, k, l;
  int start,end;
  int argmin;
  double e;
  double min;
  double sum;

  for (j = N_INPUTS+1; j <= N_NEURONS; j++){
    h = all_h[j-(N_INPUTS+1)];
    e = compute_energy();
    q[j] = exp(-e);
    sum += q[j];
  }
  q /= sum;
  // learning start
  for (j = N_INPUTS+1; j <= N_NEURONS; j++){
    // match a pattern on output layer
    h = all_h[j-(N_INPUTS+1)];
    // k-cell on input layer and i-cell on output layer
    for (k=1; k<=N_INPUTS; k++){
      for (l=N_INPUTS+1; l<=N_NEURONS; l++){
	double weight = transition_tbl
	  [static_cast<int>(X[k])]
	  [static_cast<int>(X[l])](q[j]);
	W (k,l) += weight;
	Wd(k,l) += weight;
	W (l,k) = W (k,l);
	Wd(l,k) = Wd(k,l);
	
	// if (static_cast<int>(X[k]) && static_cast<int>(X[l])){
// 	  printf("X[%d]=%lf\t",k,X[k]);
// 	  printf("X[%d]=%lf\t",l,X[l]);
// 	  printf("q[%d]=%lf\n",j,q[j]);
// 	}
      }
    }
    /* learn bias term */
    for (k=1; k<=N_NEURONS; k++){
      double weight = transition_tbl
	[static_cast<int>(X[k])][1](q[j]);
      W (k,0) += weight;
      Wd(k,0) += weight;
      W (0,k) = W (k,0);
      Wd(0,k) = Wd(k,0);
      // if (static_cast<int>(X[k])){
// 	printf("X[%d]=%lf\t",k,X[k]);
// 	printf("bias\t");
// 	printf("q[%d]=%lf\n",j,q[j]);
//       }
    }
  } /* end of for j */
  //assert(false);
  /* q[i] を累積密度関数に変換しておく */
  //  printf("g_counter=%d\n",g_counter);
  for (j = N_INPUTS+2; j <= N_NEURONS; j++){
    //    printf("q[%d]=%lf\t",j,q[j]);
    q[j] = q[j] + q[j-1];
    //    printf("q[%d]=%lf\n",j,q[j]);
  }
}

void
learning_forward_antihebb(SubPattern& v, SubPattern& h)
{
  int i, j, k, l;
  int start,end;
  int argmin;
 
  double e;
  double min;
  double sum;

  // output pattern
  //for (j = N_INPUTS+2; j <= N_NEURONS; j++){
  for (j = N_INPUTS+1; j <= N_NEURONS; j++){
    h = all_h[j-(N_INPUTS+1)];
    e = compute_energy();
    q[j] = exp(-e);
    sum += q[j];
  }
  q/=sum;
  
  for (j = N_INPUTS+1; j <= N_NEURONS; j++){
    h = all_h[j-(N_INPUTS+1)];
    /* k番目の素子（入力層）と l番目の素子の学習 */
    for (k=1; k<=N_INPUTS; k++){
      for (l=N_INPUTS+1; l<=N_NEURONS; l++){
	double weight = transition_tbl
	  [static_cast<int>(X[k])]
	  [static_cast<int>(X[l])](q[j]);
	W (k,l) -= weight;
	Wm(k,l) += weight;
	W (l,k)  = W(k,l);
	Wm(l,k)  = Wm(k,l);
	// if (static_cast<int>(X[k]) &&
	// 	    static_cast<int>(X[l])){
	// 	  printf("forward anti-hebb\n");
	// 	  printf("X[%d]=%lf\t",k,X[k]);
	// 	  printf("X[%d]=%lf\t",l,X[l]);
	// 	  printf("q[%d]=%lf\n",j,q[j]);
	// 	}
      }
    }
    /* 各素子のバイアス項の学習 */
    for (k=1; k<=N_NEURONS; k++){
      double weight = transition_tbl
	[static_cast<int>(X[k])][1](q[j]);
      W (0,k) -= weight;
      Wm(0,k) += weight;
      W (k,0)  = W (0,k);
      Wm(k,0)  = Wm(0,k);
      // if (static_cast<int>(X[k])){
      // 	printf("forward anti-hebb\n");
      // 	printf("X[%d]=%lf\t",k,X[k]);
      // 	printf("bias\t");
      // 	printf("q[%d]=%lf\n",j,q[j]);
      //       }
    }
  }/* end of for j */

  /* q[i] を累積密度関数に変換しておく */
  //printf("anti-hebb mode");
  for (j = N_INPUTS+2; j <= N_NEURONS; j++){
    //printf("q[%d]=%lf\t",j,q[j]);
    q[j] = q[j] + q[j-1];
    //printf("q[%d]=%lf\n",j,q[j]);
  }
  // if (g_counter++ < 5){
//     for (size_t i=0; i<q.size(); i++){
//       printf("q[%d]=%lf\n",i,q[i]);
//     }
//     assert(false);
//   }
}

/* すべての出力パターン q[i] の確率をもとに，各入力素子の活動を更新 */
void generate_input_from_upper_layer(SubPattern &v, SubPattern &h)
{
  int i,k;
  double r;
  r = drand48();
  i = N_INPUTS+1;
  while (q[i]<r){
    i++;
  }
  k = i;
  assert(k >= N_INPUTS+1 && k <= N_NEURONS);

  h = all_h[k-(N_INPUTS+1)];
  for (i = 1; i <= N_INPUTS; i++){
    update_a_unit(i);
  }
}

void
generate_2d_input_pattern(SubPattern &v, SubPattern &h)
{
  int d = (int)(drand48()*(double)N_INPUTS) + 1;
  v = all_v[d-1];
}

double update_a_unit(int k)
{
  int i;
  double u = 0.0;
  double p;

  for (i=0; i<=N_NEURONS; i++) {
    u += W(k,i)*X[i];
  }

  p = 1.0/(1.0 + exp(-u));

  if ( drand48() < p ){ 
    X[k] = 1.0; 
  }
  else{
    X[k] = 0.0;
  }
  return p;
}


const char *g_terminal_suffix[] = {
  "eps",
  "gif",
  "dat",
};

void
init_gnuplot ()
{
  sprintf (buf, "gnuplot -geometry %dx%d", WINDOW_SIZE,WINDOW_SIZE);
  gp = popen(buf,"w");
  switch (g_terminal_idx){
  case 0:
    fprintf(gp, "set terminal postscript eps color \"Times\" 20\n");
    break;
  case 1:
    fprintf(gp, "set terminal gif\n");
    break;
  case 2:
    fprintf(gp, "set terminal x11\n");
    break;
  default:
    break;
  }
  //fprintf(gp, "set cbrange[-2.5:2.5]\n");
  fprintf(gp, "set cbrange[-3.0:3.0]\n");
  fprintf(gp, "set yrange [] reverse\n");
  fprintf(gp, "unset key\n");
  fprintf(gp, "set palette defined (0 \"white\", 1 \"blue\")\n"); 
}

typedef void (*cmdsfunc)();
void _cmds0()
{
  fprintf(gp, "plot '-' matrix with image title 'W '\n");  
  write_multi_images(4,4,W);
}

void _cmds1()
{
  // multi image
  fprintf(gp, "set multiplot layout 1,4\n");
  fprintf(gp, "set cbrange[-2.5:2.5]\n");
  fprintf(gp, "plot '-' matrix with image title 'W '\n");
  write_multi_images(4,4,W);
  fprintf(gp, "set cbrange[-0.01:0.01]\n");
  fprintf(gp, "plot '-' matrix with image title 'Wa'\n");
  write_multi_images(4,4,Wa);
  fprintf(gp, "set cbrange[-0.01:0.01]\n");
  fprintf(gp, "plot '-' matrix with image title 'Wd'\n");
  write_multi_images(4,4,Wd);
  fprintf(gp, "set cbrange[-0.01:0.01]\n");
  fprintf(gp, "plot '-' matrix with image title 'Wm'\n");
  write_multi_images(4,4,Wm);
  fprintf(gp, "unset multiplot\n");
}

void _cmds2()
{
  // data
  fprintf(gp, "set multiplot layout 1,4\n");
  fprintf(gp, "set cbrange[-2.5:2.5]\n");
  fprintf(gp, "plot '-' matrix with image title 'W '\n");
  write_2d_data(W );
  fprintf(gp, "set cbrange[-0.01:0.01]\n");
  fprintf(gp, "plot '-' matrix with image title 'Wa'\n");
  write_2d_data(Wa);
  fprintf(gp, "set cbrange[-0.01:0.01]\n");
  fprintf(gp, "plot '-' matrix with image title 'Wd'\n");
  write_2d_data(Wd);
  fprintf(gp, "set cbrange[-0.01:0.01]\n");
  fprintf(gp, "plot '-' matrix with image title 'Wm'\n");
  write_2d_data(Wm);
  fprintf(gp, "unset multiplot\n");
}

void _cmds3()
{
  // multi images
  fprintf(gp, "set multiplot layout 2,2\n");
  fprintf(gp, "plot '-' matrix with image title 'W '\n");
  write_multi_images(4,4,W);
  fprintf(gp, "plot '-' matrix with image title 'Wa'\n");
  write_multi_images(4,4,Wa);
  fprintf(gp, "plot '-' matrix with image title 'Wd'\n");
  write_multi_images(4,4,Wd);
  fprintf(gp, "plot '-' matrix with image title 'Wm'\n");
  write_multi_images(4,4,Wm);
  fprintf(gp, "unset multiplot\n");
}

void _cmds4()
{
  // data
  fprintf(gp, "set multiplot layout 2,2\n");
  fprintf(gp, "plot '-' matrix with image title 'W '\n");
  write_2d_data(W );
  fprintf(gp, "plot '-' matrix with image title 'Wa'\n");
  write_2d_data(Wa);
  fprintf(gp, "plot '-' matrix with image title 'Wd'\n");
  write_2d_data(Wd);
  fprintf(gp, "plot '-' matrix with image title 'Wm'\n");
  write_2d_data(Wm);
  fprintf(gp, "unset multiplot\n");
}

void _cmds5()
{
  fprintf(gp, "set multiplot layout 2,2\n");
  fprintf(gp, "plot '-' matrix with image title 'W '\n");
  write_2d_all(W );
  fprintf(gp, "plot '-' matrix with image title 'Wa'\n");
  write_2d_all(Wa);
  fprintf(gp, "plot '-' matrix with image title 'Wd'\n");
  write_2d_all(Wd);
  fprintf(gp, "plot '-' matrix with image title 'Wm'\n");
  write_2d_all(Wm);
  fprintf(gp, "unset multiplot\n");
}

void _cmds6()
{
  fprintf(gp, "set multiplot layout 1,4\n");
  fprintf(gp, "plot '-' matrix with image title 'W '\n");
  write_2d_all(W );
  fprintf(gp, "plot '-' matrix with image title 'Wa'\n");
  write_2d_all(Wa);
  fprintf(gp, "plot '-' matrix with image title 'Wd'\n");
  write_2d_all(Wd);
  fprintf(gp, "plot '-' matrix with image title 'Wm'\n");
  write_2d_all(Wm);
  fprintf(gp, "unset multiplot\n");
}

cmdsfunc cmdtbl[]={
  _cmds0,
  _cmds1,
  _cmds2,
  _cmds3,
  _cmds4,
  _cmds5,
  _cmds6,
};

void
plot3d_gnuplot_file(char *title, int t)
{
  static cmdsfunc plot_cmds = cmdtbl[g_cmds_idx];
  double n_matrix = 4.0;
  double x_offset = 0.0;
#define y_offset(n) ((n_matrix - (n))*(1.0/n_matrix))
  // fprintf(gp, "set output 'dir-dat/%s-t%d.%s'\n"
  // 	  ,g_prefix, t, g_terminal_suffix[g_terminal_idx]);
  fprintf(gp, "set output 'dir-dat/%s-t%d.%s'\n"
	  ,g_prefix
	  ,t
	  ,g_terminal_suffix[g_terminal_idx]);
  fprintf(gp, "set rmargin %lf\n", 2.5);
  fprintf(gp, "set lmargin %lf\n", 1.8);
  fprintf(gp, "unset key\n");
  fprintf(gp, "unset xtics\n");
  fprintf(gp, "unset ytics\n");
  fprintf(gp, "set size square\n");
  fprintf(gp, "set title \"%s\"\n", title);
  plot_cmds();
  fflush (gp);
}

void write_multi_images(int nx, int ny, Layer& W)
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
      img[a] = W(i,j);
    }
  }

  //  fprintf(gp, "set terminal postscript enhacned color eps\n");
  for (i=0; i < n_pixels; i++){
    fprintf(gp, "%.8lf ", img[i]);
    if ((i+1)%width == 0){
      fprintf(gp, "\n");
    }
  }
  fprintf(gp,"e\n");
  fprintf(gp,"e\n");
}
