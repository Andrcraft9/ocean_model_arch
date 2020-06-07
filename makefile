## no implicit rules
.SUFFIXES: 

# Definitions and Options:
FCGCC = -fopenmp -cpp -dM -ffree-line-length-0 -I./ -Imacros/
FCINTEL = -assume byterecl -openmp -cpp -I./ -Imacros/
# Debug: (from book Introduction to Programming with Fortran)
# ?: -fPIC
FCGCC_DEBUG = -g -O -fcheck=all -finit-real=nan -Warray-temporaries -fbacktrace -pedantic-errors -Wunderflow -ffpe-trap=zero,overflow,underflow
FCINTEL_DEBUG = -g -O -finit-real=nan
# Fast:
# ?: --ffast-math -auto -stack_temp
FCGCC_FAST = -O3
FCINTEL_FAST = -O3

################# USER SECTION BEGIN ##########################################
### Default compiler (GCC or Intel / Debug or Production):
FCD = mpif90 $(FCGCC) $(FCGCC_DEBUG)
#FCD = mpif90 $(FCGCC) $(FCGCC_FAST)
#FCD = mpiifort $(FCINTEL) $(FCINTEL_DEBUG)
#FCD = mpiifort $(FCINTEL) $(FCINTEL_FAST)

### Profiler compiler:
FCPROF = /home/andr/lib/tau-2.29/x86_64/bin/tau_f90.sh -O3
######## TAU GUIDE ########
# Need to:
# Configure Tau with:
# ./configure -c++=mpicxx -cc=mpicc -fortran=mpif90 -mpi -openmp -opari -iowrapper -bfd=download -dwarf=download -otf=download -unwind=download -pdt=/home/andr/lib/pdtoolkit-3.25.1
# Set in .bashrc:
# export TAU_MAKEFILE=/home/andr/lib/tau-2.29/x86_64/lib/Makefile.tau-mpi-pdt-openmp-opari
# export PATH=$PATH:/home/andr/lib/tau-2.29/x86_64/bin
# Set in env:
# export TAU_OPTIONS='-optTauSelectFile=select.tau -optVerbose -optCompInst'
# export TAU_OMPT_SUPPORT_LEVEL=full
# export TAU_OMPT_RESOLVE_ADDRESS_EAGERLY=1
# export TAU_PROFILE=1
# export TAU_TRACE=1 
# Compile:
# make profiler
# Run:
# mpirun -n 1 model
# Pack profiler:
# paraprof --pack tau.ppk
# Pack trace
# tau_treemerge.pl
# tau2slog2 tau.trctau.edf -o tau.slog2
################# USER SECTION END ############################################

## sources for the program
SHARED = \
	shared/kind.f90 \
	shared/system.f90 \
	shared/kernel_runtime.f90 \
	shared/constants.f90 \
	shared/mpp/mpp.f90 \
	shared/configs/basinpar.f90 \
	shared/configs/sw.f90 \
	shared/errors.f90 \
	shared/debug.f90 \
	shared/mpp/hilbert_curve.f90

LEGACY = \
	legacy/service/input_output_data.f90 \
	legacy/service/rw_ctl_file.f90 \
	legacy/service/read_write_parameters.f90 \
	legacy/service/time_tools.f90

CORE = \
	core/math_tools.f90 \
	core/decomposition.f90 \
	core/data_types.f90 \
	shared/mpp/sync.f90 \
	core/kernel_interface.f90 \
	core/grid.f90 \
	core/ocean.f90

TOOLS = \
	tools/io.f90 \
	tools/time_manager.f90

SERVICE = \
	service/gridcon.f90 \
	service/grid_parameters.f90 \
	service/basinpar_construction.f90

# Kernel Layer
PHYSICS = \
	kernel/shallow_water/depth.f90 \
	kernel/shallow_water/vel_ssh.f90 \
	kernel/shallow_water/mixing.f90

# Parallel System Layer
INTERFACE = \
	interface/shallow_water/sw_interface.f90

# Algorithm Layer
CONTROL = \
	control/init_data.f90 \
	control/output.f90 \
	control/shallow_water/shallow_water.f90

## main and clean targets
model: $(subst .f90,.o, $(SHARED) $(LEGACY) $(CORE) $(TOOLS) $(SERVICE) $(PHYSICS) $(INTERFACE) $(CONTROL) model.f90)
	$(FCD) -o $@ $+

.PHONY: clean
clean:
	find . -name "*.o" -delete
	find . -name "*.mod" -delete
	find . -name "*.inst.f90" -delete
	find . -name "*.pomp.f90" -delete

## compilation rules
%.o %.mod: %.f90
	$(FCD) -c -o $*.o $<

# Profiler TAU
profiler:
	 $(FCPROF) -o model $(SHARED) $(LEGACY) $(CORE) $(TOOLS) $(SERVICE) $(PHYSICS) $(INTERFACE) $(CONTROL) model.f90 

pack_profiler:
	paraprof --pack LOG_TAU.ppk
	rm profile.*

pack_trace:
	tau_treemerge.pl
	tau2slog2 tau.trc tau.edf -o LOG_TAU.slog2
	rm *.trc
	rm *.edf

## .o -> .mod of the modules it uses
#main.o: one.mod
#one@proc.o: two.mod
#two@proc.o: one.mod

## .o of a submodule -> .smod of its parent module
#one@proc.o: one.smod
#two@proc.o: two.smod

