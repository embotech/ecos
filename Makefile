# Makefile for ECOS
# Configuration of make process in ecos.mk

include ecos.mk
C = $(CC) $(CFLAGS) -Iinclude -Iexternal/ldl/include -Iexternal/amd/include -Iexternal/SuiteSparse_config
TEST_INCLUDES = -Itest -Itest/generated_tests

# Compile all C code, including the C-callable routine
.PHONY: all
all: ldl amd ecos demo

# build Tim Davis' sparse LDL package
.PHONY: ldl
ldl:
	( cd external/ldl    ; $(MAKE) )
	$(AR) -x external/ldl/libldl.a

# build Tim Davis' AMD package
.PHONY: amd
amd:
	( cd external/amd    ; $(MAKE) )
	$(AR) -x external/amd/libamd.a

# build ECOS
ECOS_OBJS = ecos.o kkt.o cone.o spla.o timer.o preproc.o splamm.o equil.o
.PHONY: ecos
ecos: $(ECOS_OBJS)
	$(ARCHIVE) libecos.a $(ECOS_OBJS) amd_*.o ldl*.o
	- $(RANLIB) libecos.a

.PHONY: ecos_bb
ecos_bb: ldl amd ecos ecos_bb/bb_test.c
	$(C) -o ecos_bb_test ecos_bb/bb_test.c libecos.a $(LIBS)

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
.PHONY: demo
demo: ldl amd ecos src/runecos.c
	$(C) -o runecos src/runecos.c libecos.a $(LIBS)
	echo ECOS successfully built. Type ./runecos to run demo problem.

# Shared library
shared: ldl amd ecos
	$(C) -shared -o $(SHAREDNAME) ecos.o kkt.o cone.o preproc.o spla.o splamm.o timer.o equil.o -lldl -lamd -Lexternal/amd/ -Lexternal/ldl/ $(LIBS)

# ECOS tester
TEST_OBJS = qcml_utils.o norm.o sq_norm.o sum_sq.o quad_over_lin.o inv_pos.o
.PHONY: test
test: ldl amd ecos test/ecostester.c $(TEST_OBJS)
	$(C) $(TEST_INCLUDES) -o ecostester test/ecostester.c libecos.a $(LIBS) $(TEST_OBJS)

qcml_utils.o: test/generated_tests/qcml_utils.c test/generated_tests/qcml_utils.h
	$(C) $(TEST_INCLUDES) -c test/generated_tests/qcml_utils.c -o $@

norm.o: test/generated_tests/norm/norm.c test/generated_tests/norm/norm.h
	$(C) $(TEST_INCLUDES) -c test/generated_tests/norm/norm.c -o $@

quad_over_lin.o: test/generated_tests/quad_over_lin/quad_over_lin.c test/generated_tests/quad_over_lin/quad_over_lin.h
	$(C) $(TEST_INCLUDES) -c test/generated_tests/quad_over_lin/quad_over_lin.c -o $@

sq_norm.o: test/generated_tests/sq_norm/sq_norm.c test/generated_tests/sq_norm/sq_norm.h
	$(C) $(TEST_INCLUDES) -c test/generated_tests/sq_norm/sq_norm.c -o $@

sum_sq.o: test/generated_tests/sum_sq/sum_sq.c test/generated_tests/sum_sq/sum_sq.h
	$(C) $(TEST_INCLUDES) -c test/generated_tests/sum_sq/sum_sq.c -o $@

inv_pos.o: test/generated_tests/inv_pos/inv_pos.c test/generated_tests/inv_pos/inv_pos.h
	$(C) $(TEST_INCLUDES) -c test/generated_tests/inv_pos/inv_pos.c -o $@

# remove object files, but keep the compiled programs and library archives
.PHONY: clean
clean:
	( cd external/ldl    ; $(MAKE) clean )
	( cd external/amd    ; $(MAKE) clean )
	- $(RM) $(CLEAN)

# clean, and then remove compiled programs and library archives
.PHONY: purge
purge: clean
	( cd external/ldl    ; $(MAKE) purge )
	( cd external/amd    ; $(MAKE) purge )
	- $(RM) libecos.a runecos
