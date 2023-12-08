########################################
# Makefile for PX425 2021 Assignment 5 #
########################################

# C compiler and flags
CC=gcc -O0 -march=native -fopenmp
CFLAGS=-Wall -Wextra 

# Command to use for linking and executable
LD=gcc 
EXE=wifi
LDFLAGS=-lm

all: wifi

# Default build target
wifi : wifi.o mt19937ar.o
	$(CC) $(CFLAGS) -o $(EXE) wifi.o mt19937ar.o $(LDFLAGS)

# Purge build files and executable
clean :
	rm -rf *.o *.mod ./$(EXE) 

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<
