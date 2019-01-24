LBFGS_INC=/Users/pjweggy/Documents/tools/lbfgs/include
LBFGS_LIB=/Users/pjweggy/Documents/tools/lbfgs/lib
CC = g++
#GLPK = glpk-4.34
#LIBPATH = $(GLPK)/src

#ifneq ($(MAKECMDGOALS), CPLEX)
#CFLAGS = -O3 -Wall -IARG -I.. -I.
#LIBPATH =
#else
#CFLAGS = -O3 -lpthread -Wall -IARG -I.. -I.
#LIBPATH = 
#endif

CFLAGS = -O3 -Wall -I${LBFGS_INC} -L${LBFGS_LIB} -Wl,-rpath,${LBFGS_LIB}
LIBPATH = 

all: PedMix 


Utils.o : Utils.cpp  Utils.h
	$(CC)  $(CFLAGS) -c Utils.cpp

Utils2.o : Utils2.cpp  Utils2.h Utils.h
	$(CC)  $(CFLAGS) -c Utils2.cpp

Utils3.o : Utils3.cpp  Utils3.h Utils.h
	$(CC)  $(CFLAGS) -c Utils3.cpp

Utils4.o : Utils4.cpp  Utils4.h
	$(CC)  $(CFLAGS) -c Utils4.cpp

Utils-basic.o : Utils-basic.cpp  Utils-basic.h
	$(CC)  $(CFLAGS) -c Utils-basic.cpp

UtilsNumerical.o : UtilsNumerical.cpp UtilsNumerical.h
	$(CC)  $(CFLAGS) -c UtilsNumerical.cpp

#UnWeightedGraph.o : UnWeightedGraph.cpp  UnWeightedGraph.h
#	$(CC)  $(CFLAGS) -c UnWeightedGraph.cpp

FFTZ2n.o : FFTZ2n.cpp
	$(CC)  $(CFLAGS) -c FFTZ2n.cpp
    
PedigreeMixLikelihood.o : PedigreeMixLikelihood.cpp PedigreeMixLikelihood.h
	$(CC)  $(CFLAGS) -c PedigreeMixLikelihood.cpp

PedMixTest.o : PedMixTest.cpp PedMixTest.h
	$(CC)  $(CFLAGS) -c PedMixTest.cpp

PedigreeMixIO.o : PedigreeMixIO.cpp PedigreeMixIO.h
	$(CC)  $(CFLAGS) -c PedigreeMixIO.cpp
    
MixModelExplorer.o : MixModelExplorer.cpp MixModelExplorer.h
	$(CC)  $(CFLAGS) -c MixModelExplorer.cpp
        
# Now final link
main.o: main.cpp
	$(CC)  $(CFLAGS) -c main.cpp

PedMix: Utils.o Utils2.o  Utils3.o Utils4.o Utils-basic.o UtilsNumerical.o FFTZ2n.o \
        PedigreeMixLikelihood.o PedMixTest.o PedigreeMixIO.o MixModelExplorer.o  main.o
	$(CC) $(CFLAGS) -o PedMix main.o Utils.o Utils2.o  Utils3.o Utils4.o     \
    Utils-basic.o  UtilsNumerical.o FFTZ2n.o \
    PedigreeMixLikelihood.o PedMixTest.o PedigreeMixIO.o MixModelExplorer.o -llbfgs
#clean
clean:
	$(RM)  *.o PedMix
