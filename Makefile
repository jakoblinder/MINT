
# -fbounds-check causes a weird error due to non-lazy evaluation of bolean in gfortran
F77= gfortran -Wall -fno-automatic -ffixed-line-length-none # -fbounds-check

## For debugging uncomment the following
DEBUG= -O -ggdb -ffpe-trap=invalid,zero,overflow -fcheck=all -g -O -finit-real=nan -fbacktrace

FF=$(F77) $(OPT) $(DEBUG)


INCLUDE =


vpath %.f ../
vpath %.c ../
vpath %.h ../

%.o: %.f $(INCLUDE)
	$(FF) -c $<


gauss: integrate_functions.o mint-integrator.o functions.o cernroutines.o random.o
	$(FF) $^ -o $@


clean:
	rm -f *.o *.mod gauss


integrate_functions.o: functions.o
integrate_functions.o: mint-integrator.o
