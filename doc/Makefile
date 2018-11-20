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
#  Description: makefile for Doxygen documentation.                     #
#                                                                       #
#  Written by Pedro Brandimarte, Dec 2014.                              #
#  Instituto de Fisica                                                  #
#  Universidade de Sao Paulo                                            #
#  e-mail: brandimarte@gmail.com                                        #
#  ***************************** HISTORY *****************************  #
#  Original version:    December 2014                                   #
#  *******************************************************************  #

# Remove command.
RM = /bin/rm -f

.DEFAULT_GOAL := all

default: all

all:
	doxygen doxyfile

clean:
	@echo "==> Cleaning documentation files"
	$(RM) output_err
	( cd html && $(RM) *.png *.html *.css *.js && \
	cd search && $(RM) *.png *.html *.css *.js )
	( cd latex && $(MAKE) clean )

cleanall:
	@echo "==> Cleaning ALL documentation files"
	( $(RM) -r html latex output_err )

