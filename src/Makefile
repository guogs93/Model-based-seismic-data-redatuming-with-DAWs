include ./Makefile.inc

EXEDIR = ../bin

.SUFFIXES: .f .o
.f.o:
	$(F77) $(OPTF) $(INC) -c $*.f

.SUFFIXES: .f90 .o
.f90.o:
	$(F90) $(OPTF) $(INC) -c $*.f90

.SUFFIXES: .c .o
.c.o:
	$(CC) $(OPTC) $(INC) -c $*.c

all: compile

compile: main_redatum_DAW

clean:
	$(RM) *.o *.mod *~
	
SUB_v1 = libfourier.o kwsinc.o main_redatum_DAW.o minresModule.o sub_redatum_DAW.o sublinearminres.o \
         subacqui.o subreadpar.o submodeling2.o subresiduals.o subname.o para_range.o subwavefieldresample.o

main_redatum_DAW : $(SUB_v1)
	$(F90) $(OPTF) $(SUB_v1) -o $(EXEDIR)/main_redatum_DAW
