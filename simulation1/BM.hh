#ifndef BM_H_
#define BM_H_

#include "Layer.hh"

typedef double (*annealing_ptr)(double temperature);
typedef Layer (*init_layer_ptr)(int D, double scale);

// using FBM
typedef double (*energy_ptr)(Pattern& s, Layer& W);
typedef double (*denergy_ptr)(size_t i, Pattern& s, Layer& W);
typedef double (*prop_ptr)(int i, Pattern& s, Layer& W, double temperature);
// using GBM
typedef double (*energy_gbm_ptr)(Pattern &v,Pattern &h,SubLayer &L,SubLayer &J,SubLayer &W);
typedef double (*denergy_gbm_ptr)(size_t i, Pattern &v,Pattern &h,SubLayer &L,SubLayer &J,SubLayer &W);

// using RBM. will expanding

extern double g_init_alpha;		// append 2012/01/21
extern size_t g_annealing_schedule_idx;
extern size_t g_boltzmann_neural_field_idx; // append 2012/01/21
extern size_t g_denergy_idx;	// use at fbm,gbm table
extern size_t g_energy_idx ;	// use at fbm,gbm table
extern size_t g_init_layer_idx;
extern size_t g_neighbor_k;	// append 2012/01/21
extern size_t g_prop_idx   ;	// using fbm only
extern unsigned g_seed;


extern annealing_ptr	g_annealing_schedule_tbl[];
extern init_layer_ptr   g_init_layer_tbl[];
// FBM
extern energy_ptr	g_energy_tbl [];
extern denergy_ptr	g_denergy_tbl[];
// GBM
extern energy_gbm_ptr	g_energy_gbm_tbl [];
extern denergy_gbm_ptr	g_denergy_up_gbm_tbl[];
extern denergy_gbm_ptr	g_denergy_dw_gbm_tbl[];
// RBM

template<typename T> void
pmat(T& m, std::ostream& os = std::cout)
{
  using namespace std;
  for (size_t i=0; i<m.size1(); i++){
    ub::matrix_row<T> mr(m,i);
    os << mr << "\n";
  }
}

// General using Functions
void pvec(Pattern v, std::ostream& os);
double trace(Layer& W);
void pprobability(Probability& Pr, std::ostream& os);
Pattern dtob(int n, size_t size);
int btod(const Pattern& x);
double sigma(double x);
double KLd(Probability& p, Probability& q);

double annealing_schedule(double time);
Layer init_layer(int D, double scale);

// General using SubFunctions
double _annealing_schedule0(double time);
double _annealing_schedule1(double time);
double _annealing_schedule2(double time);
double _annealing_schedule3(double time);
double _annealing_schedule4(double time);

Layer _init_layer0(int D, double scale);
Layer _init_layer1(int D, double scale);
Layer _init_layer2(int D, double scale);
Layer _init_layer3(int D, double scale);
Layer _init_layer4(int D, double scale);
Layer _init_layer5(int D, double scale);
Layer _init_layer6(int D, double scale);

// FBM
double energy(Pattern& v, Layer& L);
double denergy(size_t i, Pattern& v, Layer& L);
double prop(int i, Pattern& v, Layer& L, double temperature);
double partition_function(Layer& L, Patterns& all_v);
double boltzmann_distribution(Pattern& v,Layer& L, Patterns& all_v);
void kstep_gibbs(size_t K, Pattern& v, Layer& L, double temperature);
void kstep_gibbs_count_onbit(size_t K, Probability& vf, Layer& L, double temperature, Probability &onbit_count);
void gibbs_sampling(Pattern& v, Layer& L, double temperature, Probability &onbit_count);
// FBM -sub functions
double _energy0(Pattern& v, Layer& L);
double _energy1(Pattern& v, Layer& L);
double _energy2(Pattern& v, Layer& L);
double _energy3(Pattern& v, Layer& L);
double _energy4(Pattern& v, Layer& L);
double _energy5(Pattern& v, Layer& L);
double _energy6(Pattern& v, Layer& L);
double _energy7(Pattern& v, Layer& L);
double _energy8(Pattern& v, Layer& L); // collect code candidate
double _energy9(Pattern& v, Layer& L); // hold water
double _energy10(Pattern& v, Layer& L);// maybe same as _energy9

double _denergy0(size_t i, Pattern& v, Layer& L);
double _denergy1(size_t i, Pattern& v, Layer& L);
double _denergy2(size_t i, Pattern& v, Layer& L);// correct code candidate
double _denergy3(size_t i, Pattern& v, Layer& L);
double _denergy4(size_t i, Pattern& v, Layer& L);// I think this is hold water

double _prop0(int i, Pattern& v, Layer& L, double temperature);
double _prop1(int i, Pattern& v, Layer& L, double temperature);
double _prop2(int i, Pattern& v, Layer& L, double temperature);
double _prop3(int i, Pattern& v, Layer& L, double temperature);
double _prop4(int i, Pattern& v, Layer& L, double temperature);

// GBM
double energy(Pattern &v, Pattern &h, SubLayer &L, SubLayer &J, SubLayer &W);
double denergy_up(size_t j,Pattern &v,Pattern &h,SubLayer &L,SubLayer &J,SubLayer &W);
double denergy_dw(size_t i,Pattern &v,Pattern &h,SubLayer &L,SubLayer &J,SubLayer &W);
double partition_function(SubLayer& L,SubLayer& J,SubLayer& W,Patterns& all_v,Patterns& all_h);
double prob_v(Pattern& v, SubLayer& L, SubLayer& J, SubLayer& W, Patterns& all_v, Patterns& all_h);	// -> boltzmann distribution
// GBM -sub functions
double _energy0(Pattern &v, Pattern &h, SubLayer &L, SubLayer &J, SubLayer &W);
double _denergy_up0(size_t j,Pattern &v,Pattern &h,SubLayer &L,SubLayer &J,SubLayer &W);
double _denergy_dw0(size_t i,Pattern &v,Pattern &h,SubLayer &L,SubLayer &J,SubLayer &W);

// RBM. will extending later ...






#endif	// BM_H_
