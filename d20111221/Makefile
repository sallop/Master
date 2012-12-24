# macro's
CC     = g++
OBJS   = FBM.o fbm_view.o date_main.o kurata_main.o kurata2_main.o
SRCS   = $(OBJS, .o=.cc)
#SRC_DATEIMP = $(OBJS_DATEIMP, .o=.cc)
#OBJS_DATEIMP = FBM.o date_main.o kurata_main.o
CFLAGS = -Wall -g
DEPS   = FBM.hh fbm_view.hh

.cc.o:
	$(CC) $(CFLAGS) -c $<

date : FBM.o date_main.o
	$(CC) -c FBM.cc -g -DDATEIMP
	$(CC) $(CFLAGS) FBM.o date_main.o -o $@


kurata : FBM.hh FBM.o kurata_main.o
	$(CC) $(CFLAGS) $^ -o $@

kurata2 : FBM.o kurata2_main.o
	$(CC) $(CFLAGS) $^ -o $@

kurata3 : FBM.o kurata3_main.o
	$(CC) $(CFLAGS) $^ -o $@

fbmview : fbm_view.o
	$(CC) $(CFLAGS) $^ -o $@

energy : energy_cmp.cc FBM.cc FBM.hh
	$(CC) $(CFLAGS) $^ -o $@


FBM.o : FBM.cc FBM.hh
date_main.o : date_main.cc
kurata_main.o : kurata_main.cc
energy_cmp.o : energy_cmp.cc FBM.hh

# symmetric and diagonal matrix is 0
kurata2_main.o : kurata2_main.cc

# symmetric and diagonal matrix is self expection. 
# K step gibbs is replaced K=N
kurata3_main.o : kurata3_main.cc

fbm_view.o : fbm_view.cc







clean:
	rm *.o