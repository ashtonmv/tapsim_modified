#########################################################################
#                                                                       #
# Vector class template library                                         #
#                                                                       #
# Copyright (C) 2011 Christian Oberdorfer                               #
#                                                                       #
# This library is free software: you can redistribute it and/or modify  #
# it under the terms of the GNU General Public License as published by  #
# the Free Software Foundation, either version 3 of the License, or any #
# any later version.                                                    #
#                                                                       #
# This library is distributed in the hope that it will be useful,       #
# but WITHOUT ANY WARRANTY; without even the implied warranty of        #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
# GNU General Public License for more details.                          #
#                                                                       #
# You should have received a copy of the GNU General Public License     #
# along with this library.  If not, see 'http://www.gnu.org/licenses'.  #
#                                                                       #
#########################################################################

####### Compiler, tools and options

CXX           = g++
CXXFLAGS      = -O3 -Wall $(DEFINES)
INCPATH       = -I. -I/usr/include

AR       := ar rv

COPY     := cp -rf
MOVE     := mv -f
DELETE   := rm -rf

MAKE_DIR := mkdir -p

######## tasks

GCH_OBJECTS := staticVector.h.gch \
               staticVector_tuple.h.gch \
               mathVector.h.gch \
               mathVector_2d.h.gch \
               mathVector_3d.h.gch \
               dynamicVector.h.gch

OBJECTS := staticVector_tuple.o \
           mathVector_2d.o \
           mathVector_3d.o \
           dynamicVector.o

first: all

####### Build rules

all: libvector.a $(GCH_OBJECTS)

libvector.a: $(OBJECTS)
	$(AR) $@ $^

clean:
	$(DELETE) $(GCH_OBJECTS) $(OBJECTS)
	$(DELETE) libvector.a

####### Compile

staticVector.h.gch: staticVector.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) $<
staticVector_tuple.h.gch: staticVector_tuple.h staticVector.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) $<
mathVector.h.gch: mathVector.h staticVector.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) $<
mathVector_2d.h.gch: mathVector_2d.h mathVector.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) $<
mathVector_3d.h.gch: mathVector_3d.h mathVector.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) $<
dynamicVector.h.gch: dynamicVector.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) $<

mathVector_2d.o: mathVector_2d.cpp
	$(CXX) -c $(CXXFLAGS) -fno-implicit-templates $(INCPATH) $< -o $@
mathVector_3d.o: mathVector_3d.cpp
	$(CXX) -c $(CXXFLAGS) -fno-implicit-templates $(INCPATH) $< -o $@
staticVector_tuple.o: staticVector_tuple.cpp
	$(CXX) -c $(CXXFLAGS) -fno-implicit-templates $(INCPATH) $< -o $@
dynamicVector.o: dynamicVector.cpp
	$(CXX) -c $(CXXFLAGS) -fno-implicit-templates $(INCPATH) $< -o $@