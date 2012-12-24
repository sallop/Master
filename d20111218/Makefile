# macro's
CC     = g++
OBJS   = FBM.o fbm_view.o date_main.o kurata_main.o kurata2_main.o
SRCS   = $(OBJS, .o=.cc)
CFLAGS = -Wall -g
DEPS   = FBM.hh fbm_view.hh

OBJS_DATEIMP = FBM.o date_main.o kurata_main.o
SRC_DATEIMP = $(OBJS_DATEIMP, .obj=.cc)

date : FBM.o date_main.o
	$(CC) $(CFLAGS) $(LIBS) $^ -o $@

kurata : FBM.o kurata_main.o
	$(CC) $(CFLAGS) $(LIBS) $^ -o $@

kurata2 : FBM.o kurata2_main.o
	$(CC) $(CFLAGS) $(LIBS) $^ -o $@

kurata3 : FBM.o kurata3_main.o
	$(CC) $(CFLAGS) $(LIBS) $^ -o $@




fbmview : fbm_view.o

FBM.o : FBM.hh
date_main.o : date_main.cc
kurata_main.o : kurata_main.cc FBM.hh
kurata2_main.o : kurata2_main.cc FBM.hh # symmetric matrix and diag(W)=0
kurata3_main.o : kurata3_main.cc FBM.hh # K = N impliment

fbm_view.o : fbm_view.cc

%.o:%.cc $(DEPS)
	$(CC) $(CFLAGS) -c $<

clean:
	rm *.o