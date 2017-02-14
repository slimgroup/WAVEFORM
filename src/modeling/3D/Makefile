CC=mex
MEXOPTS='/home/slim/curtd/mexopts/icc_glnxa64.xml'
CFLG='$(CFLAGS) -fopenmp -O3 -std=c99 -fPIC -funsafe-math-optimizations -march=native -ffast-math -fassociative-math'  
CFLG1='$(CLFAGS) -fopenmp -O3 -std=c99 -fPIC '
LDFLG='$(LDFLAGS) -fopenmp'
MEXFLAGS= -largeArrayDims -lmwblas


all: Helm3dmvp_forw_mex Helm3dmvp_forw_deriv_mex Helm3dmvp_adj_mex Helm3dmvp_adj_deriv_mex Helm3dmvp_forw_wnc_mex Helm3dmvp_adj_wnc_mex

Helm3d : Helm3d_27pt_cpp.cpp

Helm3dmvp_forw_mex : Helm3dmvp_mex.c
	$(CC) $^ -output $@ CFLAGS=$(CFLG) LDFLAGS=$(LDFLG) 

Helm3dmvp_forw_wnc_mex : Helm3dmvp_mex.c
	$(CC) $(MEXFLAGS) $^ -DWNCMPLX -output $@ CFLAGS=$(CFLG) LDFLAGS=$(LDFLG) 

Helm3dmvp_forw_deriv_mex : Helm3dmvp_mex.c
	$(CC) $(MEXFLAGS) $^ -DDERIV -output $@ CFLAGS=$(CFLG) LDFLAGS=$(LDFLG)	

Helm3dmvp_adj_mex : Helm3dmvp_mex.c
	$(CC) $(MEXFLAGS) $^ -DADJ -output $@ CFLAGS=$(CFLG) LDFLAGS=$(LDFLG)	

Helm3dmvp_adj_wnc_mex : Helm3dmvp_mex.c
	$(CC) $(MEXFLAGS) $^ -DADJ -DWNCMPLX -output $@ CFLAGS=$(CFLG) LDFLAGS=$(LDFLG)	

Helm3dmvp_adj_deriv_mex : Helm3dmvp_mex.c
	$(CC) $(MEXFLAGS) $^ -DADJ -DDERIV -output $@ CFLAGS=$(CFLG) LDFLAGS=$(LDFLG)
