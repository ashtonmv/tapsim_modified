####### Compiler, tools and options

DEFINES  := -DLINUX

CC       := gcc
CFLAGS   := -O3 -w $(DEFINES)
INCPATH  := -I. -I/usr/include

AR       := ar rv
DELETE   := rm -rf

OBJECTS := predicates_single.o \
           predicates_double.o

first: all

####### Build rules

all: libpredicates.a

libpredicates.a: $(OBJECTS)
	$(AR) $@ $^
clean:
	$(DELETE) $(OBJECTS)
	$(DELETE) libpredicates.a

####### Compile

predicates_single.o: predicates.c
	$(CC) -c $(CFLAGS) -DSINGLE $(INCPATH) -o $@ $<
predicates_double.o: predicates.c
	$(CC) -c $(CFLAGS) $(INCPATH) -o $@ $<
