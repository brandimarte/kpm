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
#                               makefile                                #
#  *******************************************************************  #
#  Description: makefile to build the compiler dependent MKL interface  #
#  for computing discrete Fourier transform through the fast Fourier    #
#  transform algorithm on distributed-memory architectures (clusters).  #
#                                                                       #
#  Written by Pedro Brandimarte, Jan 2015.                              #
#  Instituto de Fisica                                                  #
#  Universidade de Sao Paulo                                            #
#  e-mail: brandimarte@gmail.com                                        #
#  ***************************** HISTORY *****************************  #
#  Original version:    January 2015                                    #
#  *******************************************************************  #

# For VPATH operation.
INCFLAGS:=-I$(VPATH) $(INCFLAGS)

# Remove command.
RM = /bin/rm -f

#  *******************************************************************  #

mkl_dfti.o: 
	$(FC) $(FFLAGS) -c $(INCFLAGS) $(FPPFLAGS) $(VPATH)/$*.f90

#  *******************************************************************  #

# Include copying operations in rule to make
# sure that they are always performed.
default: module

#  *******************************************************************  #

module: libmkl_dfti.a
	cp libmkl_dfti.a ..
	@cp -p *.*d ..
#	@cp *.mod ..
#	@touch module_built

#  *******************************************************************  #

check_flags:
	@echo "In MKLDFT, INCFLAGS is: $(INCFLAGS)"

#  *******************************************************************  #

libmkl_dfti.a: check_flags mkl_dfti.o
	@$(AR) cru libmkl_dfti.a mkl_dfti.o
	-$(RANLIB) libmkl_dfti.a

#  *******************************************************************  #

clean:
	$(RM) *.o *.a *.mod *~ \#~ .\#*

#  *******************************************************************  #
