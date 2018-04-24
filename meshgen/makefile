#########################################################################
#                                                                       #
# MeshGen - a mesh generator utility for TapSim                         #
#                                                                       #
# Copyright (C) 2012 Christian Oberdorfer                               #
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

DEFINES  :=

CXX      := g++
CXXFLAGS := -g -O0 -Wall $(DEFINES)
INCPATH  :=  -I. -I/usr/include -I../tetgen

LINK     := g++
LNKFLAGS :=
LIBPATH  := -L/usr/lib64

COPY     := cp -rf
MOVE     := mv -f
DELETE   := rm -rf

MAKE_DIR := mkdir -p

####### Files

OBJECTS := meshgen.o

LIBS    := -ltetgen

vpath libtetgen.a ../tetgen

####### Build rules

all: meshgen

meshgen: meshgen.o $(LIBS)
	$(LINK) $(LNKFLAGS) $(LIBPATH) -o $@ $^

clean:
	$(DELETE) $(OBJECTS)
	$(DELETE) ./meshgen

####### Compile

meshgen.o: meshgen.cpp
	$(CXX) $(CXXFLAGS) $(INCPATH) -c $< -o $@
