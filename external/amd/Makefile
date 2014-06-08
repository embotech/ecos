#-------------------------------------------------------------------------------
# AMD Makefile for compiling on Unix systems (for GNU make only)
#-------------------------------------------------------------------------------

default: libamd.so

include ../../ecos.mk

C = $(CC) $(CFLAGS) -Iinclude -I../SuiteSparse_config

#-------------------------------------------------------------------------------
# source files
#-------------------------------------------------------------------------------

AMD = amd_aat amd_1 amd_2 amd_dump amd_postorder amd_post_tree amd_defaults \
	amd_order amd_control amd_info amd_valid amd_preprocess

CONFIG = ../SuiteSparse_config/SuiteSparse_config.h

INC = include/amd.h include/amd_internal.h $(CONFIG)

#-------------------------------------------------------------------------------
# object files for each version
#-------------------------------------------------------------------------------

#AMDI = $(addsuffix .o, $(subst amd_,amd_i_,$(AMD)))
AMDL = $(addsuffix .o, $(subst amd_,amd_l_,$(AMD)))

#-------------------------------------------------------------------------------
# compile each int and long routine (with no real/complex version)
#-------------------------------------------------------------------------------

amd_global.o: src/amd_global.c $(INC)
	$(C) -fPIC -c $< -o $@

#amd_i_%.o: src/amd_%.c $(INC)
#	$(C) -fPIC -DDINT -c $< -o $@

amd_l_%.o: src/amd_%.c $(INC)
	$(C) -fPIC -DDLONG -c $< -o $@

#-------------------------------------------------------------------------------
# Create the libamd.a library (C versions only)
#-------------------------------------------------------------------------------

libamd.so: amd_global.o $(AMDI) $(AMDL)
	$(C) -shared -o libamd.so $^

#-------------------------------------------------------------------------------
# Remove all but the files in the original distribution
#-------------------------------------------------------------------------------

clean:
	- $(RM) $(CLEAN)

purge: 
	- $(RM) libamd.so libamdf77.so

