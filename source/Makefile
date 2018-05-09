# Makefile for PoMiN

CC=gcc
CFLAGS= -std=gnu99 

PROGRAM=pomin.c
EXECPATH=/usr/local/bin/
MANPATH=/usr/local/man/
DOCDIR= ./doc
README=./README

# Get the OS
UNAME := $(shell uname)

# Set OS-specific flags
ifeq ($(UNAME), Darwin)
CFLAGS += -DDARWIN
else
CC=gcc
LIBFLAGS = -lrt
#CC=icc          OR run "make <target> CC=icc" from the command line.
endif

# Set host-specific flags
HOST := $(shell hostname)
ifeq ($(HOST), godel.ph.utexas.edu)
CC=/usr/local/gnu-devtools/bin/gcc
LIBFLAGS += -L/usr/local/gnu-devtools/lib64
endif

############################################################
# TARGETS
############################################################
.PHONY: double quad conv install docs clean header

# Compile with double precision.
double: $(PROGRAM)
	$(CC) $(CFLAGS) $(PROGRAM) -lm $(LIBFLAGS) -o pomin

# Compile with quad precision.
quad: $(PROGRAM)
	$(CC) $(CFLAGS) $(PROGRAM) -lquadmath $(LIBFLAGS) -DQUAD -o pomin

# Compile convergence calculation file.
conv: ./validation/conv/convergence.c
	$(CC) $(CFLAGS) -o ./validation/conv/convergence.exe ./validation/conv/convergence.c -lquadmath

# This is only setup for systems that conform to the Filesystem Hierarchy Standard (FHS) (i.e., Linux distributions).
install: $(DOCDIR)/pomin.1 double quad
	sudo cp pomin $(EXECPATH)
	sudo cp pomin.1 $(MANPATH)/man1

# Create man page from the README file.
docs: $(README)
	pandoc $(README) -s -o README.html -o $(DOCDIR)/pomin.1
	cd $(DOCDIR) && pdflatex pomin.tex && bibtex pomin.aux &&\
		pdflatex pomin.tex && pdflatex pomin.tex

# Remove some files generated through usual testing procedures.
clean:
	git clean -dXf

# Create a header for n particles. To execute this target run: "make header N=n"
header:
	@header="description,start time,end time,timestep,iterations,courant number,gravitational constant,speed of light" ; \
	i=1 ; while [[ $$i -le ${N} ]] ; do \
		header=$$header",mass_$$i,qx_$$i,qy_$$i,qz_$$i,px_$$i,py_$$i,pz_$$i" ; \
		((i = i + 1)) ; \
	done ; echo $$header
