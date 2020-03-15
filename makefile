## no implicit rules
.SUFFIXES: 

## definitions
FC = mpifort
FCFLAGS = -cpp -dM -Wall -fPIC -fcheck=all -ffree-line-length-0 -O3 -Wtabs -fopenmp -I ./

## sources for the program
SHARED = \
	shared/kind.f90 \
	shared/constants.f90 \
	shared/mpp/mpp.f90 \
	shared/configs/basinpar.f90 \
	shared/errors.f90 \
	shared/debug.f90 \
	shared/mpp/hilbert_curve.f90

LEGACY = \
	legacy/service/input_output_data.f90 \
	legacy/service/rw_ctl_file.f90

CORE = \
	core/math_tools.f90 \
	core/decomposition.f90 \
	core/data_types.f90 \
	shared/mpp/sync.f90 \
	core/kernel_interface.f90 \
	core/grid.f90 \
	core/ocean.f90

TOOLS = \
	tools/io.f90

SERVICE = \
	service/gridcon.f90 \
	service/grid_parameters.f90 \
	service/basinpar.f90

# Kernel Layer
PHYSICS = \
	#physics/velocity.f90

# Parallel System Layer
INTERFACE = \
	#interface/ocean_interface.f90

# Algorithm Layer
CONTROL = \
	control/init_data.f90 \
	control/output.f90

## main and clean targets
model: $(subst .f90,.o, $(SHARED) $(LEGACY) $(CORE) $(TOOLS) $(SERVICE) $(PHYSICS) $(INTERFACE) $(CONTROL) model.f90)
	$(FC) $(FCFLAGS) -o $@ $+

.PHONY: clean
clean:
	find . -name "*.o" -delete
	find . -name "*.mod" -delete

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

