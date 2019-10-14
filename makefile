
FC = mpifort
FCFLAGS = -cpp -dM -Wall -fPIC -fcheck=all -ffree-line-length-0 -O3 -Wtabs -fopenmp

SHARED = \
	shared/kind.f90 \
	shared/mpp.f90 \
	shared/errors.f90 \
	shared/constants.f90 \
	shared/basinpar.f90

CORE = \
	core/decomposition.f90 \
	core/data_types.f90 \
	core/grid.f90 \
	core/ocean.f90

TOOLS = \
	tools/io.f90

CONTROL = \
	control/init_data.f90

all: compile

clean:
	rm model *.mod *.o

compile:
	$(FC) $(FCFLAGS) -o model $(SHARED) $(CORE) $(TOOLS) $(CONTROL) model.f90 