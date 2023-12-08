########################################
# Makefile for PX425 2021 Assignment 4 #
########################################

# C compiler and flags
CC=mpiicc 
CFLAGS=-O3 -Wall -Wextra

# Command to use for linking and executable
LD=$(CC)
EXE=rfim
LDFLAGS= -lpng -lm

OBJECTS=mt19937ar.o grid.o comms.o makePNG.o rfim.o

# Default build target
perc : $(OBJECTS)
	$(LD) $(CFLAGS) -o $(EXE) $(OBJECTS) $(LDFLAGS)


# Purge build files and executable
clean :
	rm -f *.o *.mod $(EXE)

# Rules for building object files
%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<
