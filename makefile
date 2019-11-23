## no implicit rules
.SUFFIXES: 

## definitions
FC = mpifort
FCFLAGS = -cpp -dM -Wall -fPIC -fcheck=all -ffree-line-length-0 -O3 -Wtabs -fopenmp

## sources for the program
SHARED = \
	shared/kind.f90 \
	shared/mpp/mpp.f90 \
	shared/mpp/hilbert_curve.f90 \
	shared/errors.f90 \
	shared/constants.f90 \
	shared/basinpar.f90

CORE = \
	core/decomposition.f90 \
	core/data_types.f90 \
	core/kernel_interface.f90 \
	core/grid.f90 \
	core/ocean.f90

TOOLS = \
	tools/io.f90

PHYSICS = \
	physics/velocity.f90

CONTROL = \
	control/init_data.f90 \
	control/output.f90 \
	control/ocean_model.f90

## main and clean targets
model: $(subst .f90,.o, $(SHARED) $(CORE)) #$(TOOLS) $(PHYSICS) $(CONTROL))
	$(FC) $(FCFLAGS) -o $@ $+

.PHONY: clean
clean:
	-rm -rf *.o *.mod

## compilation rules
%.o %.mod: %.f90
	$(FC) $(FCFLAGS) -c -o $*.o $<

## .o -> .mod of the modules it uses
#main.o: one.mod
#one@proc.o: two.mod
#two@proc.o: one.mod

## .o of a submodule -> .smod of its parent module
#one@proc.o: one.smod
#two@proc.o: two.smod

