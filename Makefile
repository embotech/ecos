# Makefile for ECOS
# Configuration of make process in ecos.mk

include ecos.mk
CFLAGS += -Iinclude -Iexternal/ldl/include -Iexternal/amd/include -Iexternal/SuiteSparse_config
TEST_INCLUDES = -Itest -Itest/generated_tests

# Compile all C code, including the C-callable routine
.PHONY: all
all: ecos ecos_bb demo

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
ECOS_OBJS = ecos.o kkt.o cone.o spla.o ctrlc.o timer.o preproc.o splamm.o equil.o
.PHONY: ecos
ecos: $(ECOS_OBJS) ldl amd
	$(ARCHIVE) libecos.a $(ECOS_OBJS) amd_*.o ldl*.o
	- $(RANLIB) libecos.a

# build ECOS branch-and-bound
ECOS_BB_OBJS = $(ECOS_OBJS) ecos_bb_preproc.o ecos_bb.o
.PHONY: ecos_bb
ecos_bb: $(ECOS_BB_OBJS) ldl amd
	$(ARCHIVE) libecos_bb.a $(ECOS_BB_OBJS) amd_*.o ldl*.o
	- $(RANLIB) libecos_bb.a

ecos_bb.o: ecos_bb/ecos_bb.c include/ecos_bb.h
	$(CC) $(CFLAGS) -c ecos_bb/ecos_bb.c -o ecos_bb.o

ecos_bb_preproc.o: ecos_bb/ecos_bb_preproc.c include/ecos_bb.h
	$(CC) $(CFLAGS) -c ecos_bb/ecos_bb_preproc.c -o ecos_bb_preproc.o

ecos.o: src/ecos.c include/ecos.h
	$(CC) $(CFLAGS) -c src/ecos.c -o ecos.o

kkt.o: src/kkt.c include/kkt.h
	$(CC) $(CFLAGS) -c src/kkt.c -o kkt.o

cone.o: src/cone.c include/cone.h
	$(CC) $(CFLAGS) -c src/cone.c -o cone.o

preproc.o: src/preproc.c
	$(CC) $(CFLAGS) -c src/preproc.c -o preproc.o

spla.o: src/spla.c include/spla.h
	$(CC) $(CFLAGS) -c src/spla.c -o spla.o

splamm.o: src/splamm.c include/splamm.h
	$(CC) $(CFLAGS) -c src/splamm.c -o splamm.o

ctrlc.o: src/ctrlc.c include/ctrlc.h
	$(CC) $(CFLAGS) -c src/ctrlc.c -o ctrlc.o

timer.o: src/timer.c include/timer.h
	$(CC) $(CFLAGS) -c src/timer.c -o timer.o

equil.o: src/equil.c include/equil.h
	$(CC) $(CFLAGS) -c src/equil.c -o equil.o

# ECOS demo
.PHONY: demo
demo: ldl amd ecos src/runecos.c
	$(CC) $(CFLAGS) -o runecos src/runecos.c libecos.a $(LIBS)
	echo ECOS successfully built. Type ./runecos to run demo problem.

# Shared library
shared: ldl amd ecos
	$(CC) $(CFLAGS) -shared -o $(SHAREDNAME) $(ECOS_OBJS) -lldl -lamd -Lexternal/amd/ -Lexternal/ldl/ $(LIBS)

# ECOS tester
TEST_OBJS = qcml_utils.o norm.o sq_norm.o sum_sq.o quad_over_lin.o inv_pos.o
.PHONY: test
test: ldl amd ecos test/ecostester.c $(TEST_OBJS)
	$(CC) $(CFLAGS) $(TEST_INCLUDES) -o ecostester test/ecostester.c libecos.a $(LIBS) $(TEST_OBJS)

qcml_utils.o: test/generated_tests/qcml_utils.c test/generated_tests/qcml_utils.h
	$(CC) $(CFLAGS) $(TEST_INCLUDES) -c test/generated_tests/qcml_utils.c -o $@

norm.o: test/generated_tests/norm/norm.c test/generated_tests/norm/norm.h
	$(CC) $(CFLAGS) $(TEST_INCLUDES) -c test/generated_tests/norm/norm.c -o $@

quad_over_lin.o: test/generated_tests/quad_over_lin/quad_over_lin.c test/generated_tests/quad_over_lin/quad_over_lin.h
	$(CC) $(CFLAGS) $(TEST_INCLUDES) -c test/generated_tests/quad_over_lin/quad_over_lin.c -o $@

sq_norm.o: test/generated_tests/sq_norm/sq_norm.c test/generated_tests/sq_norm/sq_norm.h
	$(CC) $(CFLAGS) $(TEST_INCLUDES) -c test/generated_tests/sq_norm/sq_norm.c -o $@

sum_sq.o: test/generated_tests/sum_sq/sum_sq.c test/generated_tests/sum_sq/sum_sq.h
	$(CC) $(CFLAGS) $(TEST_INCLUDES) -c test/generated_tests/sum_sq/sum_sq.c -o $@

inv_pos.o: test/generated_tests/inv_pos/inv_pos.c test/generated_tests/inv_pos/inv_pos.h
	$(CC) $(CFLAGS) $(TEST_INCLUDES) -c test/generated_tests/inv_pos/inv_pos.c -o $@

sqrt.o: test/generated_tests/sqrt/sqrt.c test/generated_tests/sqrt/sqrt.h
	$(CC) $(CFLAGS) $(TEST_INCLUDES) -c test/generated_tests/sqrt/sqrt.c -o $@

ecos_bb_test: ecos_bb
	$(CC) $(CFLAGS) -L. -o ecos_bb_test test/bb_test.c -lecos_bb $(LIBS)

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
	- $(RM) libecos.a libecos_bb.a runecos
