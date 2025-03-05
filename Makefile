
#Compiler options
PRECOMP = gcc
PRCOMP_FLAGS = -traditional-cpp -cpp -E -P 
FC = mpif90
COMPILER_FLAGS = -g -fbacktrace -fexternal-blas

# SERIAL/PARALLEL
#PARALLEL_FLAG = -DSERIAL
PARALLEL_FLAG = -DPARALLEL

# Location of LIBSTELL
LIBSTELL_INC = -I$(STELLOPT_PATH)LIBSTELL/Release
LIBSTELL_LIB = ~/bin/libstell.a

# NETCDF (old versions use nc-config)
NETCDF_INC = $(shell nf-config --fflags)
NETCDF_LIB = $(shell nf-config --flibs)

# Pick one of the following for your machine
# BLAS/LAPACK (MACOS)
BLASLAPACK_LIB = -framework Accelerate # MacOS
# BLAS/LAPACK (MKL)
#BLASLAPACK_LIB = -I${MKLROOT}/include ${MKLROOT}/lib/libmkl_intel_lp64.a ${MKLROOT}/lib/libmkl_core.a ${MKLROOT}/lib/libmkl_sequential.a -lpthread -lm

############### DO NOT EDIT BELOW #####################################

PRCOMP_FLAGS += $(PARALLEL_FLAG)
GIT_REPO      = $(shell git ls-remote --get-url)
GIT_BRANCH    = $(shell git rev-parse --abbrev-ref HEAD)
GIT_VERSION   = $(shell git describe --abbrev=4 --dirty --always --tags)
GIT_HASH      = $(shell git show -s --format=%H)
TIME          = $(shell date +"%d.%m.%Y %H:%M:%S")
PRCOMP_FLAGS += -DGIT_VERSION_EXT="'$(GIT_VERSION)'"
PRCOMP_FLAGS += -DGIT_HASH_EXT="'$(GIT_HASH)'"
PRCOMP_FLAGS += -DGIT_BRANCH_EXT="'$(GIT_BRANCH)'"
PRCOMP_FLAGS += -DGIT_REPO_EXT="'$(GIT_REPO)'"
PRCOMP_FLAGS += -DBUILT_ON_EXT="'$(TIME)'"

FFLAGS= ${LIBSTELL_INC} ${NETCDF_INC}
LDFLAGS=${LIBSTELL_LIB} ${NETCDF_LIB} ${BLASLAPACK_LIB}

OBJ=fitpack.o Fourier_lib_convolve.o 

%.o: %.f
	$(PRECOMP) $(PRCOMP_FLAGS) $^ > temp.f
	$(FC) $(COMPILER_FLAGS) -c -o $@ temp.f $(FFLAGS)

xmetric: metric_element_create_ver8.46.o
	$(FC) $(COMPILER_FLAGS) -o $@ $^ $(LDFLAGS)

xstgap: stellgap_ver5.o $(OBJ) 
	$(FC) $(COMPILER_FLAGS) -o $@ $^ $(LDFLAGS)

xstgap_snd_ver6: stellgap_soundwave_lagrng_ver6.o $(OBJ) 
	$(FC) $(COMPILER_FLAGS) -o $@ $^ $(LDFLAGS)

xstgap_snd_ver7: stellgap_soundwave_lagrng_ver7.o $(OBJ) 
	$(FC) $(COMPILER_FLAGS) -o $@ $^ $(LDFLAGS)

clean:
	@rm -rf *.o xmetric xstgap xstgap_snd_ver6 xstgap_snd_ver7 temp.f

help:
	@echo 'BUILD OPTIONS:'
	@echo '   xmetric xstgap xstgap_snd_ver6 xstgap_snd_ver7'
	@echo 'Directories and flags for build.'
	@echo '--------------------------------'
	@echo LIBSTELL_INC is $(LIBSTELL_INC)
	@echo NETCDF_INC is $(NETCDF_INC)
	@echo LIBSTELL_LIB is $(LIBSTELL_LIB)
	@echo NETCDF_LIB is $(NETCDF_LIB)
	@echo BLASLAPACK_LIB is $(BLASLAPACK_LIB)
	@echo '--------------------------------'
	@echo Precompiler is $(PRECOMP)
	@echo Precompiler flags are $(PRCOMP_FLAGS)
	@echo '--------------------------------'
	@echo Compiler $(FC)
	@echo Compiler flags are $(COMPILER_FLAGS)
	@echo '--------------------------------'
#	@echo `python --version` `which python`
#	@echo '--------------------------------'