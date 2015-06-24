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

.PHONY: all build clean cleanall

.DEFAULT_GOAL := all

default: all

all: build

build:
	@echo "#  ***************************************************  #"
	@echo "#  KPM Fortran Code 2014                                #"
	@echo "#                                                       #"
	@echo "#  Written by:                                          #"
	@echo "#                                                       #"
	@echo "#    Eric de Castro e Andrade (eandrade@ift.unesp.br)   #"
	@echo "#    Pedro Brandimarte (brandimarte@gmail.com)          #"
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
	( cd src && $(MAKE) )

clean:
	( cd src && $(MAKE) clean )

cleanall:
	( cd src && $(MAKE) cleanall )

