#------------------------------------------------------------------------------
# AMD Makefile
#------------------------------------------------------------------------------

SUITESPARSE ?= $(realpath $(CURDIR)/..)
export SUITESPARSE

default: all

include ../SuiteSparse_config/SuiteSparse_config.mk

demos: all

# Compile all C code.  Do not compile the FORTRAN versions.
all:
	( cd lib    ; $(MAKE) )

# compile just the C-callable libraries (not Demos)
library:
	( cd lib    ; $(MAKE) )

# compile the static libraries only
static:
	( cd lib    ; $(MAKE) static )

# compile the FORTRAN libraries and demo programs (not compiled by "make all")
fortran:
	( cd lib    ; $(MAKE) fortran )

# remove object files, but keep the compiled programs and library archives
clean:
	( cd lib    ; $(MAKE) clean )

# clean, and then remove compiled programs and library archives
purge:
	( cd lib    ; $(MAKE) purge )

distclean: purge

ccode: library

lib: library

# install AMD
install:
	( cd lib  ; $(MAKE) install )

# uninstall AMD
uninstall:
	( cd lib  ; $(MAKE) uninstall )

