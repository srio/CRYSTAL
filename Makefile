# --------------------------------------------------------------
#
# Usage:
#	make            [to create the executable diff_pat]
#	make install	 [to install  it in the correct directory]
#	make clean	 [to clean the source directory]
#
#
# --------------------------------------------------------------
#   in ESRF/NICE: 
#     . /scisoft/ESRF_sw/opteron2/set_environment
#
# System definition (architecture)
#

ARCH := $(shell /scisoft/xop2.4/xop_uname)
FC = gfortran -ffree-line-length-none
#FC = g95 -ffree-line-length-huge
FFLAGS =   -static
#FFLAGS =   
COMPILEOPT=-cpp -D_COMPILE4NIX

#in ESRF/NICE: 
#   source /scisoft/ESRF_sw/opteron2/set_environment.tcsh
#   use  FC=gfortran
#


#*********** PLATFORM DEPENDENT VARIABLES ******************
#
#-------------------------- libraries for sun5 (expgh) 
#
ifeq ($(ARCH),sun5)
endif
#
#-------------------------- libraries for linux (expglio)
#
ifeq ($(ARCH),linux)
endif
#
#-------------------------- libraries for darwin
#
ifeq ($(ARCH),darwin)
endif
#*********** END OF PLATFORM DEPENDENT VARIABLES ******************

all: clean diff_pat compliance

compliance: compliance.F90 shadow_globaldefinitions.F90 stringio.F90 shadow_math.F90 elasticity.F90
	$(FC) $(FFLAGS) $(COMPILEOPT) -c shadow_globaldefinitions.F90
	$(FC) $(FFLAGS) -c shadow_math.F90
	$(FC) $(FFLAGS) -c stringio.F90
	$(FC) $(FFLAGS) -c elasticity.F90
	$(FC) $(FFLAGS) -c compliance.F90
	$(FC) $(FFLAGS) -o compliance shadow_globaldefinitions.o stringio.o shadow_math.o elasticity.o compliance.o

diff_pat: diff_pat.F90 crystal3.F90 shadow_globaldefinitions.F90 stringio.F90 shadow_math.F90 elasticity.F90
	$(FC) $(FFLAGS) $(COMPILEOPT) -c shadow_globaldefinitions.F90
	$(FC) $(FFLAGS) -c shadow_math.F90
	$(FC) $(FFLAGS) -c stringio.F90
	$(FC) $(FFLAGS) -c elasticity.F90
	$(FC) $(FFLAGS) -c crystal3.F90
	$(FC) $(FFLAGS) -c diff_pat.F90
	$(FC) $(FFLAGS) -o diff_pat shadow_globaldefinitions.o stringio.o shadow_math.o elasticity.o crystal3.o diff_pat.o   

clean: 
	/bin/rm -f *.mod *.o diff_pat diff_pat.dat diff_pat.par *.mod compliance  compliance

install:
	cp diff_pat ./../../bin.$(ARCH)/diff_pat
	cp diff_pat ./../../../xop2.3/bin.$(ARCH)/diff_pat
	cp compliance ./../../bin.$(ARCH)/compliance
	cp compliance ./../../../xop2.3/bin.$(ARCH)/compliance
