RM = /bin/rm -f

CC =  mpicc 
FC =  mpif90 
FL =  mpif90 

F77 = mpif90 
F90 = mpif90 
ICC = mpicc 

AR = ar vr
RANLIB  = echo

#DOUBLE VERSUS SINGLE PRECISION LU FACTORIZATION:
#   -D_sp                   (single precision)
#   -D_dp                   (double precision)
EXT=
OPT_PRE = -cpp $(EXT)

#OPTIONS COMPILE
 
OPTF =  -O3 -fopenmp $(OPT_PRE)
OPTC =  -O3 -fopenmp $(OPT_PRE)
OPTL =  -O3 -fopenmp $(OPT_PRE)
