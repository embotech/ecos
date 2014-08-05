# Makefile configuration for ECOS

## Intel C Compiler
#CC = icc
#CFLAGS = -O3 -m64 -Wall -strict-ansi -DLDL_LONG -DDLONG
#LIBS = -lm

## GNU C Compiler
#CC = gcc
CFLAGS = -O2 -Wall -DLDL_LONG -DDLONG -Wextra -fPIC #-ansi -ipo

UNAME := $(shell uname)

ifeq (CYGWIN, $(findstring CYGWIN, $(UNAME)))
ISWINDOWS := 1
else ifeq (MINGW, $(findstring MINGW, $(UNAME)))
ISWINDOWS := 1
else ifeq (MSYS, $(findstring MSYS, $(UNAME)))
ISWINDOWS := 1
else
ISWINDOWS := 0
endif

ifeq ($(UNAME), Darwin)
# we're on apple, no need to link rt library
LIBS = -lm
# shared library has extension .dylib
SHAREDNAME = libecos.dylib
else ifeq ($(ISWINDOWS), 1)
# we're on windows (cygwin or msys)
LIBS = -lm
# shared library has extension .dll
SHAREDNAME = libecos.dll
else
# we're on a linux system, use accurate timer provided by clock_gettime()
LIBS = -lm -lrt
# shared library has extension .so
SHAREDNAME = libecos.so
endif


## AR and RANLIB FOR GENERATING LIBRARIES
AR = ar
ARFLAGS = rcs
ARCHIVE = $(AR) $(ARFLAGS)
RANLIB = ranlib

## WHICH FILES TO CLEAN UP
CLEAN = *.o *.obj *.ln *.bb *.bbg *.da *.tcov *.gcov gmon.out *.bak *.d *.gcda *.gcno libecos.a libecos.so libecos.dylib libecos.dll ecos_bb_test ecostester ecostester.exe
