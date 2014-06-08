#-------------------------------------------------------------------------------
# Makefile for the LDL library
#-------------------------------------------------------------------------------

default: all

include ../../ecos.mk

INCLUDES = -I../SuiteSparse_config -Iinclude

C = $(CC) $(CFLAGS) $(INCLUDES)

all: libldl.so

#-------------------------------------------------------------------------------
# the ldl library:
#-------------------------------------------------------------------------------

libldl.so: src/ldl.c include/ldl.h
	$(C) -fPIC -DLDL_LONG -c src/ldl.c -o ldl.o
	$(C) -shared -o libldl.so ldl.o

purge: 
	- $(RM) libldl.so
clean:
	- $(RM) $(CLEAN)


