MPICXX := mpic++ 
LIBS := -lfftw3 -lmpi 
INCLUDE_PATH := -I../libs/fftw3/include
LIB_PATH     := -L../libs/fftw3/lib
FLAGS        := -std=c++11 -O3 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE

TSpecDyn: main.o CSpecDyn.o MikeFFT.o
	$(MPICXX) $^ -o $@ $(LIB_PATH) $(LIBS) $(FLAGS)
	
main.o: main.cc CSpecDyn.h define.h
	$(MPICXX) $(INCLUDE_PATH) $(FLAGS) -o $@ -c $<

CSpecDyn.o: CSpecDyn.cc CSpecDyn.h define.h
	$(MPICXX) $(INCLUDE_PATH) $(FLAGS) -o $@ -c $<
	
MikeFFT.o: MikeFFT.cc MikeFFT.h
	$(MPICXX) $(INCLUDE_PATH) $(FLAGS) -o $@ -c $<

clean: 
	rm -f CSpecDyn CSpecDyn_* *.o
