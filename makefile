### Pls export TAU and Nvidia HPC libs
###
###

## no implicit rules
.SUFFIXES: 

# Definitions and Options:
FCGCC = -fopenmp -cpp -dM -ffree-line-length-0 -I./ -Imacros/
FCINTEL = -assume byterecl -qopenmp -fpp -allow nofpp_comments -I./ -Imacros/
# Debug: (from book Introduction to Programming with Fortran)
# ?: -fPIC -Warray-temporaries
FCGCC_DEBUG = -g -O -fcheck=all -finit-real=nan  -fbacktrace -pedantic-errors -Wunderflow -ffpe-trap=zero,overflow,underflow
FCINTEL_DEBUG = -g -O -check all -traceback
# Fast:
# ?: --ffast-math -auto -stack_temp
FCGCC_FAST = -Ofast
FCINTEL_FAST = -O3

### PGI, GPU, CUDA
FCPGI = -mp -cpp -dM -cuda -gpu=cc60 -I./ -Imacros/
FCPGI_DEBUG = -g -C -gopt -Mneginfo=all
FCPGI_FAST = -fast -Minfo=all -Mnoopenmp -Mcache_align -Mprefetch -Mipa=inline,reshape
FCPGI_FAST = -fast
PROF_CUDA = nvprof
# PROF_CUDA = /opt/nvidia/hpc_sdk/Linux_x86_64/2021/cuda/bin/nvprof

################# USER SECTION BEGIN ##########################################
### Default compiler (GCC or Intel / Debug or Production):
#FCD = mpif90 $(FCGCC) $(FCGCC_DEBUG)
#FCD = mpif90 $(FCGCC) $(FCGCC_FAST)
#FCD = mpiifort $(FCINTEL) $(FCINTEL_DEBUG)
#FCD = mpiifort $(FCINTEL) $(FCINTEL_FAST)
### GPU compilers
#FCD = /opt/nvidia/hpc_sdk/Linux_x86_64/21.2/comm_libs/mpi/bin/mpif90 $(FCPGI) $(FCPGI_DEBUG)
FCD = /opt/nvidia/hpc_sdk/Linux_x86_64/21.2/comm_libs/mpi/bin/mpif90 $(FCPGI) $(FCPGI_FAST)
#FCD = nvfortran $(FCPGI) $(FCPGI_DEBUG)
# FCD = /opt/nvidia/hpc_sdk/Linux_x86_64/2021/compilers/bin/nvfortran $(FCPGI) $(FCPGI_FAST)

### Profiler compiler:
FCPROF = tau_f90.sh $(FCGCC) $(FCGCC_FAST)
# FCPROF = /home/andr/lib/tau-2.29/x86_64/bin/tau_f90.sh $(FCGCC) $(FCGCC_FAST)
# FCPROF = /data4t/chaplygin/lib/tau-2.28/x86_64/bin/tau_f90.sh $(FCINTEL) $(FCINTEL_FAST)

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
# export TAU_TRACK_SIGNALS=1
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
	shared/errors.f90 \
	shared/debug.f90 \
	shared/mpp/hilbert_curve.f90

LEGACY = \
	legacy/service/input_output_data.f90 \
	legacy/service/rw_ctl_file.f90 \
	legacy/service/read_write_parameters.f90 \
	legacy/service/time_tools.f90

CONFIGS = \
	configs/basinpar.f90 \
	configs/sw.f90 \
	configs/parallel.f90 \
	configs/cmd.f90

CORE = \
	core/math_tools.f90 \
	core/decomposition.f90 \
	core/data_types.f90 \
	shared/mpp/sync.f90 \
	core/grid.f90 \
	core/ocean.f90 \
	core/kernel_interface.f90

TOOLS = \
	tools/io.f90 \
	tools/time_manager.f90

# Kernel Layer
PHYSICS = \
	kernel/shallow_water/depth.f90 \
	kernel/shallow_water/vel_ssh.f90 \
	kernel/shallow_water/mixing.f90 \
	kernel/service/grid_parameters.f90 \
	kernel/service/grid_kernels.f90 \
	kernel/tracer/leapfrog_tracer.f90

# Parallel System Layer
INTERFACE = \
	interface/shallow_water/sw_interface.f90 \
	interface/service/grid_interface.f90 \
	interface/tracer/tracer_interface.f90

# Algorithm Layer
SERVICE = \
	service/gridcon.f90 \
	service/basinpar_construction.f90

GPU = \
	gpu/kernel/depth_gpu.f90 \
	gpu/kernel/mixing_gpu.f90 \
	gpu/kernel/vel_ssh_gpu.f90 \
	gpu/kernel/leapfrog_tracer_gpu.f90 \
	gpu/interface/sw_interface_gpu.f90 \
	gpu/interface/tracer_interface_gpu.f90

CONTROL = \
	control/init_data.f90 \
	control/output.f90 \
	control/shallow_water/shallow_water.f90 \
	control/preprocess.f90 \
	control/tracer.f90

## main and clean targets
model: $(subst .f90,.o, $(SHARED) $(LEGACY) $(CONFIGS) $(CORE) $(TOOLS)  $(PHYSICS) $(INTERFACE) $(SERVICE) $(GPU) $(CONTROL) model.f90)
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
	 $(FCPROF) -o model $(SHARED) $(LEGACY) $(CONFIGS) $(CORE) $(TOOLS) $(PHYSICS) $(INTERFACE) $(SERVICE) $(GPU) $(CONTROL) model.f90 

pack_profiler:
	paraprof --pack LOG_TAU.ppk
	rm profile.*

pack_trace:
	tau_treemerge.pl
	tau2slog2 tau.trc tau.edf -o LOG_TAU.slog2
	rm *.trc
	rm *.edf

run_nv_mpi:
	/opt/nvidia/hpc_sdk/Linux_x86_64/21.2/comm_libs/mpi/bin/mpirun --bind-to none -n 1 ./model

#--print-gpu-trace -o gpu.nvprof
run_nv_prof:
	/opt/nvidia/hpc_sdk/Linux_x86_64/21.2/comm_libs/mpi/bin/mpirun --bind-to none -n 1 nvprof  ./model

## .o -> .mod of the modules it uses
#main.o: one.mod
#one@proc.o: two.mod
#two@proc.o: one.mod

## .o of a submodule -> .smod of its parent module
#one@proc.o: one.smod
#two@proc.o: two.smod

