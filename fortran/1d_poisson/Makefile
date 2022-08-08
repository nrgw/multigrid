FC=gfortran
FFLAGS = -O2
EXEC = multigrid
SRC = multigrid_data.f90 multigrid_functions.f90 main.f90 

OBJ = ${SRC:.f90=.o}

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $<

$(EXEC): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

clean:
	rm *.mod *.o $(EXEC)
