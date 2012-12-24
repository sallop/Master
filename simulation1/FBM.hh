#ifndef FBM_H_
#define FBM_H_

#include "BM.hh"
#include "config.hh"
//#include "util.hh"

typedef double (*annealing_ptr)(double temperature);
typedef double (*energy_ptr)(Pattern& s, Layer& W);
typedef double (*denergy_ptr)(size_t i, Pattern& s, Layer& W);
typedef double (*prop_ptr)(int i, Pattern& s, Layer& W, double temperature);
typedef double (*generator_ptr)();
typedef Layer (*init_layer_ptr)(int D, double scale);

extern energy_ptr	g_energy_tbl [];
extern denergy_ptr	g_denergy_tbl[];
extern prop_ptr		g_prop_tbl   [];
extern annealing_ptr	g_annealing_schedule_tbl[];
extern init_layer_ptr   g_init_layer_tbl[];

extern size_t		g_energy_idx;
extern size_t		g_denergy_idx;
extern size_t		g_prop_idx;
extern size_t		g_annealing_schedule_idx;
extern PlotTable	g_plot_table;

double energy_test(Pattern& v, Layer& L, std::ostream &os);
double denergy_test(size_t i, Pattern& v, Layer& L,std::ostream &os);
double prop_test(int i, Pattern& v, Layer& L, double temperature,std::ostream& os);

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
double energy(Pattern& v, Layer& L);

double _denergy0(size_t i, Pattern& v, Layer& L);
double _denergy1(size_t i, Pattern& v, Layer& L);
double _denergy2(size_t i, Pattern& v, Layer& L);// correct code candidate
double _denergy3(size_t i, Pattern& v, Layer& L);
double _denergy4(size_t i, Pattern& v, Layer& L);// I think this is hold water
double denergy(size_t i, Pattern& v, Layer& L);

double _prop0(int i, Pattern& v, Layer& L, double temperature);
double _prop1(int i, Pattern& v, Layer& L, double temperature);
double _prop2(int i, Pattern& v, Layer& L, double temperature);
double _prop3(int i, Pattern& v, Layer& L, double temperature);
double _prop4(int i, Pattern& v, Layer& L, double temperature);
double prop(int i, Pattern& v, Layer& L, double temperature);

double partition_function(Layer& L, Patterns& all_v);
double boltzmann_distribution(Pattern& v,Layer& L, Patterns& all_v);
double KLd(Probability& p, Probability& q);
void kstep_gibbs(size_t K, Pattern& v, Layer& L,double temperature);

double _annealing_schedule0(double temperature);
double _annealing_schedule1(double temperature);
double _annealing_schedule2(double temperature);
double _annealing_schedule3(double temperature);
double _annealing_schedule4(double temperature);
double annealing_schedule(double temperature);

Layer _init_layer0(int D, double scale);
Layer _init_layer1(int D, double scale);
Layer _init_layer2(int D, double scale);
Layer _init_layer3(int D, double scale);
Layer _init_layer4(int D, double scale);
Layer _init_layer5(int D, double scale);
Layer _init_layer6(int D, double scale);
Layer init_layer(int D, double scale);

void kstep_gibbs(size_t K, Pattern& v, Layer& L, double temperature);
void kstep_gibbs_count_onbit(size_t K, Probability& vf, Layer& L,
			     double temperature, Probability &onbit_count);

// gibbs sampling
void gibbs_sampling(Pattern& v, Layer& L,
		    double temperature, Probability &onbit_count);
#endif	// FBM_H_
