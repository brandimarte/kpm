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
#  Description: main makefile for building KPM code.                    #
#                                                                       #
#  Written by Pedro Brandimarte, Dec 2014.                              #
#  Instituto de Fisica                                                  #
#  Universidade de Sao Paulo                                            #
#  e-mail: brandimarte@gmail.com                                        #
#  ***************************** HISTORY *****************************  #
#  Original version:    December 2014                                   #
#  *******************************************************************  #

KPM_DIR = .
include $(KPM_DIR)/make.inc

#  *******************************************************************  #

.PHONY: all doc build clean cleanall

.DEFAULT_GOAL := all

default: all

all: what doc build

#  *******************************************************************  #

what:
	@echo "#  ***************************************************  #"
	@echo "#  KPM Fortran Code 2014                                #"
	@echo "#                                                       #"
	@echo "#    Pedro Brandimarte (brandimarte@gmail.com)          #"
	@echo "#    Eric de Castro e Andrade (eandrade@ift.unesp.br)   #"
	@echo "#                                                       #"
	@echo "#  Copyright (c), All Rights Reserved                   #"
	@echo "#                                                       #"
	@echo "#  This program is free software. You can redistribute  #"
	@echo "#  it and/or modify it under the terms of the GNU       #"
	@echo "#  General Public License (version 3 or later) as       #"
	@echo "#  published by the Free Software Foundation            #"
	@echo "#  <http://fsf.org/>.                                   #"
	@echo "#  ***************************************************  #"
	@echo ""

#  *******************************************************************  #

# Doxygen documentation.
doc: 
	( cd doc && $(MAKE) )

#  *******************************************************************  #

# KPM code.
build:
	( cd src && $(MAKE) )

#  *******************************************************************  #

clean:
	( cd doc && $(MAKE) clean )
	( cd src && $(MAKE) clean )

cleanall:
	( cd doc && $(MAKE) cleanall )
	( cd src && $(MAKE) cleanall )

