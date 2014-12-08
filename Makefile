#  *******************************************************************  #
#  KPM Fortran Code 2014                                                #
#                                                                       #
#  Written by Eric de Castro e Andrade (eandrade@ift.unesp.br) and      #
#             Pedro Brandimarte (brandimarte@gmail.com).                #
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
#  Description: makefile to build KPM code.                             #
#                                                                       #
#  Written by Pedro Brandimarte, Dec 2014.                              #
#  Instituto de Fisica                                                  #
#  Universidade de Sao Paulo                                            #
#  e-mail: brandimarte@gmail.com                                        #
#  ***************************** HISTORY *****************************  #
#  Original version:    December 2014                                   #
#  *******************************************************************  #

.SUFFIXES: .f .F .o .a .f90 .F90

# Executable name.
EXEC = kpmdev

# Compilation architecture.
KPM_ARCH = x86_64-mac-openmpi-mkl

# Fortran compiler.
FC = mpif90
FFLAGS = -O2 -xHost -traceback -fpe0 -heap-arrays 1024 \
         -check assume,bounds,format,output_conversion,pointers,stack,uninit \
         -warn unused,truncated_source,uncalled,declarations,usage
LDFLAGS = -static-intel

# BLAS and LAPACK with MKL.
MKL           = /opt/intel/composerxe/mkl
MKL_LIBS      = $(MKL)/lib/libmkl_intel_lp64.a \
	        $(MKL)/lib/libmkl_sequential.a \
	        $(MKL)/lib/libmkl_core.a
SCALAPACK_LIBS = /opt/local/lib/libscalapack.a
MKL_INCLUDE    = -I$(MKL)/include

# MPI
MPI_LIBS     = -L/opt/local/lib -lmpi
MPI_INCLUDE  = -I/opt/local/include
FPPFLAGS_MPI = -DMPI

#  *******************************************************************  #
#  Do not change the following lines.                                   #
#  *******************************************************************  #

# All libraries.
LDLIBS = $(MKL_LIBS) $(SCALAPACK_LIBS) -lpthread -lm $(MPI_LIBS)

# All includes.
INCFLAGS = $(MKL_INCLUDE) $(MPI_INCLUDE)

# Preprocessor definitions or flags.
FPPFLAGS = $(FPPFLAGS_MPI)

# Remove command.
RM = /bin/rm -f

# All source files.
SRCS = precision.F90 parallel.F90 io.F90 options.F90 init.F90 end.F90 \
	kpm.F90

# All objects.
OBJS = $(SRCS:.F90=.o)
#OBJS = precision.o parallel.o io.o init.o end.o kpm.o

VPATH=.

#  *******************************************************************  #

# This block tells make to consider only these suffixes in its operation.
.SUFFIXES:
.SUFFIXES: .f .F .o .a .f90 .F90 .cpp

.F90.o:
	$(FC) $(FFLAGS) -c $(INCFLAGS) $(FPPFLAGS) $<

.F90:
	make $*.o
	$(FC) $(FFLAGS) -o $(INCFLAGS) $(FPPFLAGS) $* $*.o $(LDLIBS) 

#  *******************************************************************  #

default: what $(EXEC)

#  *******************************************************************  #

FDF=libfdf.a
$(FDF): 
	(cd fdf ; $(MAKE) "FC=$(FC)" "LDFLAGS=$(LDFLAGS)" \
	"VPATH=$(VPATH)/fdf" module)

# Routines using fdf calls.
init.o: $(FDF)

#  *******************************************************************  #

what:
	@echo	
	@echo "Compilation architecture to be used: ${KPM_ARCH}"
	@echo
	@echo "Hit ^C to abort..."
	@echo
	@sleep 2

$(EXEC): $(FDF) $(OBJS)
	$(FC) -o $@ $(LDFLAGS) $^ $(LDLIBS)

#  *******************************************************************  #

clean:
	@echo "==> Cleaning object, library, and executable files"
	$(RM) *~ \#~ .\#* *.o *.a *.mod $(EXEC) core a.out
	(cd fdf ; $(MAKE) clean)

#  *******************************************************************  #
