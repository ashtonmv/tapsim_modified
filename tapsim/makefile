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

DEFINES  := #-DVERBOSE #-DDEMO_MESH_MODE

CXX      := g++
CXXFLAGS := -g -O3 -Wall $(DEFINES)
INCPATH  :=  -I. -I/usr/include -I../predicates -I../tetgen -I../vector

LINK     := g++
LNKFLAGS :=
LIBPATH  := -L/usr/lib64

COPY     := cp -rf
MOVE     := mv -f
DELETE   := rm -rf

MAKE_DIR := mkdir -p

####### Files

OBJECTS := tapsim.o \
	   process.o \
	   commands.o \
	   system_3d.o \
	   configuration.o \
	   geometry_3d.o \
	   grid_3d.o \
	   surface_3d.o \
	   fieldline_3d.o \
	   trajectory_3d.o \
	   info.o \
	   file_io.o \
	   file_util.o \
	   logResults.o \
	   logTrajectory.o \
	   logGrid.o \
	   logSurface.o \
	   utils.o \
	   debug.o

LIBS    := -lpredicates \
	   -ltetgen \
	   -lpthread

vpath libtetgen.a ../tetgen
vpath libvector.a ../vector
vpath libpredicates.a ../predicates
vpath libpthread.so /usr/lib64

####### Build rules

all: tapsim

tapsim: $(OBJECTS) $(LIBS)
	$(LINK) $(LNKFLAGS) $(LIBPATH) -o $@ $^

clean:
	$(DELETE) $(OBJECTS)
	$(DELETE) ./tapsim

####### Compile

tapsim.o: tapsim.cpp version.h process.h commands.h system_3d.h configuration.h surface_3d.h trajectory_3d.h info.h file_io.h debug.h
	$(CXX) $(CXXFLAGS) $(INCPATH) -c -o $@ $<
process.o: process.cpp process.h system_3d.h configuration.h surface_3d.h trajectory_3d.h info.h file_io.h debug.h
	$(CXX) $(CXXFLAGS) $(INCPATH) -c -o $@ $<
commands.o: commands.cpp commands.h configuration.h
	$(CXX) $(CXXFLAGS) $(INCPATH) -c -o $@ $<
system_3d.o: system_3d.cpp system_3d.h configuration.h geometry_3d.h grid_3d.h
	$(CXX) $(CXXFLAGS) $(INCPATH) -c -o $@ $<
configuration.o: configuration.cpp configuration.h
	$(CXX) $(CXXFLAGS) $(INCPATH) -c -o $@ $<
geometry_3d.o: geometry_3d.cpp geometry_3d.h configuration.h utils.h debug.h
	$(CXX) $(CXXFLAGS) $(INCPATH) -c -o $@ $<
grid_3d.o: grid_3d.cpp grid_3d.h configuration.h info.h utils.h debug.h
	$(CXX) $(CXXFLAGS) $(INCPATH) -c -o $@ $<
surface_3d.o: surface_3d.cpp surface_3d.h configuration.h
	$(CXX) $(CXXFLAGS) $(INCPATH) -c -o $@ $<
fieldline_3d.o: fieldline_3d.cpp fieldline_3d.h configuration.h geometry_3d.h grid_3d.h
	$(CXX) $(CXXFLAGS) $(INCPATH) -c -o $@ $<
trajectory_3d.o: trajectory_3d.cpp configuration.h trajectory_3d.h debug.h
	$(CXX) $(CXXFLAGS) $(INCPATH) -c -o $@ $<
info.o: info.cpp info.h
	$(CXX) $(CXXFLAGS) $(INCPATH) -c -o $@ $<
file_io.o: file_io.cpp file_io.h configuration.h file_util.h
	$(CXX) $(CXXFLAGS) $(INCPATH) -c -o $@ $<
file_util.o: file_util.cpp file_util.h
	$(CXX) $(CXXFLAGS) $(INCPATH) -c -o $@ $<
logResults.o: logResults.cpp logResults.h configuration.h
	$(CXX) $(CXXFLAGS) $(INCPATH) -c -o $@ $<
logTrajectory.o: logTrajectory.cpp logTrajectory.h configuration.h
	$(CXX) $(CXXFLAGS) $(INCPATH) -c -o $@ $<
logGrid.o: logGrid.cpp logGrid.h configuration.h
	$(CXX) $(CXXFLAGS) $(INCPATH) -c -o $@ $<
logSurface.o: logSurface.cpp logSurface.h trajectory_3d.h logTrajectory.h configuration.h
	$(CXX) $(CXXFLAGS) $(INCPATH) -c -o $@ $<
utils.o: utils.cpp utils.h
	$(CXX) $(CXXFLAGS) $(INCPATH) -c -o $@ $<
debug.o: debug.cpp debug.h system_3d.h configuration.h surface_3d.h process.h
	$(CXX) $(CXXFLAGS) $(INCPATH) -c -o $@ $<
