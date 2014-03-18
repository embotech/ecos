# Makefile for ECOS
# Configuration of make process in ecos.mk

include ecos.mk
C = $(CC) $(CFLAGS) -Iinclude -Iexternal/ldl/include -Iexternal/amd/include -Iexternal/SuiteSparse_config 
TEST_INCLUDES = -Itest -Itest/quadratic

# Compile all C code, including the C-callable routine
all: ldl amd ecos demo	

# build Tim Davis' sparse LDL package
ldl: 
	( cd external/ldl    ; $(MAKE) )
	$(AR) -x external/ldl/libldl.a
	
# build Tim Davis' AMD package
amd: 
	( cd external/amd    ; $(MAKE) )
	$(AR) -x external/amd/libamd.a

# build ECOS
ecos: ecos.o kkt.o cone.o spla.o timer.o preproc.o splamm.o equil.o
	$(ARCHIVE) libecos.a *.o
	- $(RANLIB) libecos.a

ecos.o: src/ecos.c include/ecos.h
	$(C) -c src/ecos.c -o ecos.o

kkt.o: src/kkt.c include/kkt.h
	$(C) -c src/kkt.c -o kkt.o

cone.o: src/cone.c include/cone.h
	$(C) -c src/cone.c -o cone.o

preproc.o: src/preproc.c
	$(C) -c src/preproc.c -o preproc.o

spla.o: src/spla.c include/spla.h
	$(C) -c src/spla.c -o spla.o

splamm.o: src/splamm.c include/splamm.h
	$(C) -c src/splamm.c -o splamm.o

timer.o: src/timer.c include/timer.h
	$(C) -c src/timer.c -o timer.o

equil.o: src/equil.c include/equil.h
	$(C) -c src/equil.c -o equil.o

# ECOS demo
demo: ldl amd ecos src/runecos.c 
	$(C) -o runecos src/runecos.c libecos.a $(LIBS)
	echo ECOS successfully built. Type ./runecos to run demo problem.
	
# ECOS tester
QUADRATIC_TEST_OBJS = qcml_utils.o norm.o sq_norm.o sum_sq.o quad_over_lin.o
test: ldl amd ecos test/ecostester.c $(QUADRATIC_TEST_OBJS)
	$(C) $(TEST_INCLUDES) -o ecostester test/ecostester.c libecos.a $(LIBS) $(QUADRATIC_TEST_OBJS)

qcml_utils.o: test/quadratic/qcml_utils.c test/quadratic/qcml_utils.h
	$(C) $(TEST_INCLUDES) -c test/quadratic/qcml_utils.c -o $@

norm.o: test/quadratic/norm/norm.c test/quadratic/norm/norm.h
	$(C) $(TEST_INCLUDES) -c test/quadratic/norm/norm.c -o $@

quad_over_lin.o: test/quadratic/quad_over_lin/quad_over_lin.c test/quadratic/quad_over_lin/quad_over_lin.h
	$(C) $(TEST_INCLUDES) -c test/quadratic/quad_over_lin/quad_over_lin.c -o $@

sq_norm.o: test/quadratic/sq_norm/sq_norm.c test/quadratic/sq_norm/sq_norm.h
	$(C) $(TEST_INCLUDES) -c test/quadratic/sq_norm/sq_norm.c -o $@

sum_sq.o: test/quadratic/sum_sq/sum_sq.c test/quadratic/sum_sq/sum_sq.h
	$(C) $(TEST_INCLUDES) -c test/quadratic/sum_sq/sum_sq.c -o $@

# remove object files, but keep the compiled programs and library archives
clean:
	( cd external/ldl    ; $(MAKE) clean )
	( cd external/amd    ; $(MAKE) clean )
	- $(RM) $(CLEAN)

# clean, and then remove compiled programs and library archives
purge: clean
	( cd external/ldl    ; $(MAKE) purge )
	( cd external/amd    ; $(MAKE) purge )	
	- $(RM) libecos.a runecos

