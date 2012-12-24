# macro's
CC     = g++
OBJS   = FBM.o fbm_view.o date_main.o kurata_main.o kurata2_main.o util.o BM.o
SRCS   = $(OBJS, .o=.cc)
CFLAGS = -Wall -g
#CFLAGS = -Wall -O3
DEPS   = FBM.hh fbm_view.hh

.cc.o:
	$(CC) $(CFLAGS) -c $<

date : FBM.o date_main.o
	$(CC) -c FBM.cc -g -DDATEIMP
	$(CC) $(CFLAGS) FBM.o date_main.o -o $@


kurata1 : kurata_main.cc FBM.hh BM.hh BM.o FBM.o util.o
	$(CC) $(CFLAGS) $^ -o $@

kurata2 : kurata2_main.cc FBM.hh BM.hh BM.o FBM.o util.o
	$(CC) $(CFLAGS) $^ -o $@

kurata3 : kurata3_main.cc FBM.hh BM.hh BM.o FBM.o util.o
	$(CC) $(CFLAGS) $^ -o $@

# gibbs sampling. No use K parameter
kurata4 : kurata4_main.cc FBM.hh BM.hh BM.o FBM.o util.o
	$(CC) $(CFLAGS) $^ -o $@

kurata5 : kurata5_main.cc FBM.hh BM.hh BM.o FBM.o util.o
	$(CC) $(CFLAGS) $^ -o $@

# k is nothing and gibbs sampling M
kurata6 : kurata6_main.cc FBM.hh BM.hh BM.o FBM.o util.o
	$(CC) $(CFLAGS) $^ -o $@

# kurata impliments is wrong. I confuse temperature with learning rate.
# So this implimeent, separate temperature and learning rate(alpha).
# and, learning rate adapt to simulated annealing.
kurata7 : kurata7_main.cc BM.o FBM.o config.o FBM.hh BM.hh config.hh
	$(CC) $(CFLAGS) $^ -o $@

# this impliments alike kurata7.
# temperature adapt to simulated annealing.
kurata8 : kurata8_main.cc FBM.hh BM.hh BM.o FBM.o util.o
	$(CC) $(CFLAGS) $^ -o $@

# this impliments is derived kurata7.
# change a fantasy particle vf.(Pattern vf -> Patterns vf)
# vf take a t state.
kurata9 : kurata9_main.cc FBM.hh BM.hh BM.o FBM.o config.o
	$(CC) $(CFLAGS) $^ -o $@

# kurata9. dimension D is 4
kurata10 : kurata10_main.cc FBM.hh BM.hh BM.o FBM.o config.o
	$(CC) $(CFLAGS) $^ -o $@

kurata11 : kurata11_main.cc FBM.hh BM.hh BM.o FBM.o config.o
	$(CC) $(CFLAGS) $^ -o $@


# this impliments is derived kurata7.
# change a fantasy particle vf.(Pattern vf -> Patterns vf)
# vf take a t state.
fbm : kurata9_main.cc FBM.hh BM.hh BM.o FBM.o config.o
	$(CC) $(CFLAGS) $^ -o $@

# 
gbm1 : gbm_main.cc BM.o GBM.o config.o BM.hh config.hh
	$(CC) $(CFLAGS) $^ -o $@

# 
gbm2 : gbm2_main.cc BM.o GBM.o config.o BM.hh config.hh
	$(CC) $(CFLAGS) $^ -o $@


# pseudo restricted boltzmann machine
prbm1 : pseudo_rbm_main.cc BM.o GBM.o config.o BM.hh config.hh
	$(CC) $(CFLAGS) $^ -o $@

# 
prbm2 : pseudo_rbm2_main.cc BM.o GBM.o config.o BM.hh config.hh
	$(CC) $(CFLAGS) -O3 $^ -o $@


# topological mapping
topological : topological_mapping.o
	$(CC) $(CFLAGS) $^ -o $@

fbmview : fbm_view.o
	$(CC) $(CFLAGS) $^ -o $@

energy : energy_cmp.cc FBM.cc FBM.hh
	$(CC) $(CFLAGS) $^ -o $@

pltkl : plot_kl.cc BM.o FBM.o config.o
	$(CC) $(CFLAGS) -g $^ -o $@

plgbm : plot_gbm.cc BM.o GBM.o Layer.hh config.o config.hh
	$(CC) $(CFLAGS) $^ -g -o $@

draw_energy : draw_energy.cc BM.o GBM.o
	$(CC) $(CFLAGS) $^ -g -o $@


nf : neural_field.cc
	$(CC) $(CFLAGS) -g $^ -o $@

# boltzmann neural field


# algorithm is using GBM
bnf : bnf_main.cc BM.o GBM.o config.o BM.hh config.hh
	$(CC) $(CFLAGS) $^ -o $@

# algorithm 
bnf2 : bnf2_main.cc BM.o GBM.o config.o BM.hh config.hh
	$(CC) $(CFLAGS) $^ -o $@


bnf_search : bnf_search.cc BM.o GBM.o config.o BM.hh config.hh
	$(CC) $(CFLAGS) $^ -o $@

bnf2_search : bnf2_search.cc BM.o GBM.o config.o BM.hh config.hh
	$(CC) $(CFLAGS) $^ -o $@

# date algorithm impliment 1
bnffbm : bnffbm_main.cc BM.hh BM.o GBM.o config.o
	$(CC) $(CFLAGS) $^ -o $@

# date algorithm impliment 2
bnffbm2 : bnffbm2_main.cc BM.hh BM.o GBM.o config.o
	$(CC) $(CFLAGS) $^ -o $@



util.o   : util.hh
config.o : config.hh Layer.hh
BM.o     : BM.hh Layer.hh
GBM.o    : GBM.hh GBM.cc Layer.hh
FBM.o    : FBM.hh BM.cc Layer.hh

date_main.o : FBM.o
kurata_main.o : FBM.o
energy_cmp.o : FBM.o
neural_field.o : neural_field.cc
plot_kl.o : util.o FBM.o util.hh
topological_mapping.o : topological_mapping.cc


# symmetric and diagonal matrix is 0
kurata2_main.o : kurata2_main.cc

# symmetric and diagonal matrix is self expection. 
# K step gibbs is replaced K=N
kurata3_main.o : kurata3_main.cc

fbm_view.o : fbm_view.cc

boltzmann_neuralfield.o : boltzmann_neuralfield.cc


clean:
	rm *.o