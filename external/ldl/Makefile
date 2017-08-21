#------------------------------------------------------------------------------
# LDL Makefile
#------------------------------------------------------------------------------

SUITESPARSE ?= $(realpath $(CURDIR)/..)
export SUITESPARSE

default: all

include ../SuiteSparse_config/SuiteSparse_config.mk

# Compile all C code, including the C-callable routine.
all:
	( cd lib    ; $(MAKE) )

# compile just the C-callable libraries
library:
	( cd lib    ; $(MAKE) )

# compile the static libraries only
static:
	( cd lib    ; $(MAKE) static )

# remove object files, but keep the compiled programs and library archives
clean:
	( cd lib    ; $(MAKE) clean )

# clean, and then remove compiled programs and library archives
purge:
	( cd lib    ; $(MAKE) purge )

distclean: purge

ccode: library

lib: library

# install LDL
install:
	( cd lib ; $(MAKE) install )

# uninstall LDL
uninstall:
	( cd lib ; $(MAKE) uninstall )

