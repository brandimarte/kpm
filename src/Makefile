#  *******************************************************************  #
#  KPM Fortran Code 2014                                                #
#                                                                       #
#  By Pedro Brandimarte (brandimarte@gmail.com) and                     #
#     Eric de Castro e Andrade (eandrade@ift.unesp.br)                  #
#                                                                       #
#  Copyright (c), All Rights Reserved                                   #
#                                                                       #
#  This program is free software. You can redistribute it and/or        #
#  modify it under the terms of the GNU General Public License          #
#  (version 3 or later) as published by the Free Software Foundation    #
#  <http://fsf.org/>.                                                   #
#                                                                       #
#  This program is distributed in the hope that it will be useful, but  #
#  WITHOUT ANY WARRANTY, without even the implied warranty of           #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU     #
#  General Public License for more details (file 'LICENSE_GPL'          #
#  distributed along with this program or at                            #
#  <http://www.gnu.org/licenses/gpl.html>).                             #
#  *******************************************************************  #
#                               Makefile                                #
#  *******************************************************************  #
#  Description: src directory makefile to build KPM code.               #
#                                                                       #
#  Written by Pedro Brandimarte, Dec 2014.                              #
#  Instituto de Fisica                                                  #
#  Universidade de Sao Paulo                                            #
#  e-mail: brandimarte@gmail.com                                        #
#  ***************************** HISTORY *****************************  #
#  Original version:    December 2014                                   #
#  *******************************************************************  #

DIR = .
KPM_DIR = ..
include $(KPM_DIR)/make.inc

# If not set, use default variable names.
EXEC ?= kpm_exec

# Check for MKL directory and give user hint how to set it.
ifeq ($(MKL),)
    ifeq ($(MKLROOT),)
        vars_sh  := $(shell which icc 2>/dev/null | perl -pe            \
         's|/bin.*/icc|/mkl/bin/mklvars.sh|')
        vars_csh := $(shell which icc 2>/dev/null | perl -pe            \
         's|/bin.*/icc|/mkl/bin/mklvars.csh|')
        ifeq ($(vars_sh),)
            vars_sh  := /opt/intel/composerxe/mkl/bin/mklvars.sh
            vars_csh := /opt/intel/composerxe/mkl/bin/mklvars.csh
        endif
        $(error Set $$MKLROOT, preferably in your environment, e.g.,    \
         run "source $(vars_sh) intel64" in ~/.bashrc, or               \
         "source $(vars_csh) intel64" in ~/.cshrc)
    endif

    ifeq ($(wildcard $(MKLROOT)),)
        $(error $$MKLROOT=$(MKLROOT) does not exist. Please set         \
         $$MKLROOT to where MKL is installed.)
    endif
endif

#  *******************************************************************  #

# This block tells make to consider only these suffixes in its operation.
.SUFFIXES:
.SUFFIXES: .f .F .o .a .f90 .F90 .cpp

# Rules to compile object files.
.F90.o:
	$(FC) $(FFLAGS) -c $(INCFLAGS) $(FPPFLAGS) $<

.F90:
	make $*.o
	$(FC) $(FFLAGS) -o $(INCFLAGS) $(FPPFLAGS) $* $*.o $(LDLIBS) 

#  *******************************************************************  #

# All libraries.
LDLIBS = $(MKL_LIBS) $(SCALAPACK_LIBS) $(OTHER_LIBS) $(MPI_LIBS)

# All includes.
INCFLAGS = $(MKL_INCLUDE) $(MPI_INCLUDE)

# Preprocessor definitions or flags.
FPPFLAGS = $(FPPFLAGS_MPI)

# Remove command.
RM = /bin/rm -f

# All source files.
SRCS = precision.F90 parallel.F90 io.F90 string.F90 options.F90         \
	random.F90 lattice.F90 init.F90 hsparse.F90 hstd.F90 hadt.F90   \
	hlm.F90 moment.F90 kernel.F90 dct.F90 output.F90 end.F90 kpm.F90

# All objects.
OBJS = $(SRCS:.F90=.o)

#  *******************************************************************  #

.PHONY: all clean cleanall

.DEFAULT_GOAL := all

all: default

default: what $(EXEC)

#  *******************************************************************  #

FDF=libfdf.a
$(FDF): 
	( cd fdf ; $(MAKE) "FC=$(FC)" "LDFLAGS=$(LDFLAGS)" \
	"VPATH=$(DIR)/fdf" module )

# Routines using fdf calls.
#init.o options.o: $(FDF)

#  *******************************************************************  #

# Creates MKL_DFTI modules:
MKLDFT=libmkl_dfti.a
$(MKLDFT):
	@echo "$(MKL)"
	@cp -f $(MKL)/include/mkl_dfti.f90 $(DIR)/MKLDFT
	@chmod a+x $(DIR)/MKLDFT/*
	@chmod u+w $(DIR)/MKLDFT/*
	( cd MKLDFT ; $(MAKE) "FC=$(FC)" "FFLAGS=$(FFLAGS)" \
	"LDFLAGS=$(LDFLAGS)" "FPPFLAGS=$(FPPFLAGS)" \
	"INCFLAGS=$(INCFLAGS)" "AR=$(AR)" "RANLIB=$(RANLIB)" \
	"VPATH=$(DIR)" module )

# Routines using mkl_dfti calls.
dct.o : $(MKLDFT)

#  *******************************************************************  #

what:
	@echo	
	@echo "Compilation architecture to be used: ${KPM_ARCH}"
	@echo
	@echo "Hit ^C to abort..."
	@echo
	@sleep 2

$(EXEC): $(FDF) $(MKLDFT) $(OBJS)
	$(FC) -o $(KPM_DIR)/$@ $(LDFLAGS) $(OBJS) $(MKLDFT) $(FDF)      \
	$(LDLIBS)

#  *******************************************************************  #

clean:
	@echo "==> Cleaning KPM object and executable files"
	$(RM) *~ \#~ .\#* *.o *.mod *.log $(EXEC) core a.out

cleanall: clean
	@echo "==> Cleaning libraries compilations"
	$(RM) *.a
	( cd fdf ; $(MAKE) clean )
	@if [ -d MKLDFT ] ; then ( cd MKLDFT && $(MAKE) clean ) ; fi

#  *******************************************************************  #
