#########################################################################
#                                                                       #
# TAPSim - an atom probe data simulation program                        #
#                                                                       #
# Copyright (C) 2011 Christian Oberdorfer                               #
#                                                                       #
# This program is free software: you can redistribute it and/or modify  #
# it under the terms of the GNU General Public License as published by  #
# the Free Software Foundation, either version 3 of the License, or any #
# any later version.                                                    #
#                                                                       #
# This program is distributed in the hope that it will be useful,       #
# but WITHOUT ANY WARRANTY; without even the implied warranty of        #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
# GNU General Public License for more details.                          #
#                                                                       #
# You should have received a copy of the GNU General Public License     #
# along with this program.  If not, see 'http://www.gnu.org/licenses'.  #
#                                                                       #
#########################################################################

####### Compiler, tools and options

COPY     := cp -rf
MOVE     := mv -f
DELETE   := rm -rf

MAKE_DIR := mkdir -p

####### Targets

LIB_TARGETS := libpredicates.a \
               libtetgen.a \
               libvector.a

BIN_DIR   := ./bin

.PHONY: tapsim meshgen $(LIB_TARGETS)

first: all

####### Build rules

all: tapsim meshgen $(LIB_TARGETS)
	$(MAKE_DIR) $(BIN_DIR)
	$(MOVE) ./tapsim/tapsim $(BIN_DIR)
	$(MOVE) ./meshgen/meshgen $(BIN_DIR)

clean:
	$(MAKE) --directory=./predicates clean
	$(MAKE) --directory=./tetgen clean
	$(MAKE) --directory=./vector clean
	$(MAKE) --directory=./tapsim clean
	$(MAKE) --directory=./meshgen clean
	$(DELETE) $(BIN_DIR)/tapsim
	$(DELETE) $(BIN_DIR)/meshgen

####### Distributed compilation

tapsim: $(LIB_TARGETS)
	$(MAKE) --directory=./tapsim all

meshgen: libtetgen.a
	$(MAKE) --directory=./meshgen all

libpredicates.a:
	$(MAKE) --directory=./predicates all

libtetgen.a:
	$(MAKE) --directory=./tetgen all

libvector.a:
	$(MAKE) --directory=./vector all


