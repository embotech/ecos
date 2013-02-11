# Makefile configuration for ECOS

## Intel C Compiler
#CC = icc
#CFLAGS = -O3 -m64 -strict-ansi
#LIBS = -lm

## GNU C Compiler
#CC = gcc
CFLAGS = -O3 -DLDL_LONG -DDLONG

UNAME := $(shell uname)

ifeq ($(UNAME), Linux)
# we're on a linux system, use accurate timer provided by clock_gettime()
LIBS = -lm -lrt
else
# we're on apple, no need to link rt library
LIBS = -lm
endif


## AR and RANLIB FOR GENERATING LIBRARIES
AR = ar
ARFLAGS = rv
ARCHIVE = $(AR) $(ARFLAGS)
RANLIB = ranlib

## WHICH FILES TO CLEAN UP
CLEAN = *.o *.obj *.ln *.bb *.bbg *.da *.tcov *.gcov gmon.out *.bak *.d *.gcda *.gcno
