
FC = mpifort
FCFLAGS = -cpp -dM -Wall -fPIC -fcheck=all -ffree-line-length-0 -O3 -Wtabs -fopenmp

BASE = \
	base/kind.f90 \
	base/parallel.f90 \
	base/errors.f90 \
	base/basinpar.f90 \
	base/decomposition.f90 \
	base/data_types.f90

DATA = \
	data/grid.f90 \
	data/ocean.f90

TOOLS = \
	tools/io.f90

CONTROL = \
	control/init_data.f90

clean:
	rm model *.mod *.o

compile:
	$(FC) $(FCFLAGS) -o model $(BASE) $(DATA) $(TOOLS) $(CONTROL) model.f90 

test:
	$(FC) -o model $(BASE) $(DATA) model.f90 