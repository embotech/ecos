#-------------------------------------------------------------------------------
# Makefile for the LDL library
#-------------------------------------------------------------------------------

default: all

include ../../ecos.mk

INCLUDES = -I../SuiteSparse_config -Iinclude

C = $(CC) $(CFLAGS) $(INCLUDES)

all: libldl.a

#-------------------------------------------------------------------------------
# the ldl library:
#-------------------------------------------------------------------------------

libldl.a: src/ldl.c include/ldl.h ../../include/glblopts.h ../../include/ecos.h
	$(C) -c src/ldl.c -o ldl.o
	$(C) -DLDL_LONG -c src/ldl.c -o ldll.o
	$(ARCHIVE)  libldl.a ldl.o ldll.o
	- $(RANLIB) libldl.a

purge:
	- $(RM) libldl.a
clean:
	- $(RM) $(CLEAN)
