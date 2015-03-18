# Makefile for ECOS

# Disable implicit rules such as 'make ecos' when 'ecos.o' exists
.SUFFIXES:

# Configuration of make process in ecos.mk
include ecos.mk
CFLAGS += -Iinclude -Iexternal/ldl/include -Iexternal/amd/include -Iexternal/SuiteSparse_config 
TEST_INCLUDES = -Itest -Itest/generated

# Compile all C code, including the C-callable routine
.PHONY: all
all: libecos.a libecos_bb.a runecos runecosexp

# build Tim Davis' sparse LDL package
$(LDL):
	( cd external/ldl    ; $(MAKE) )
	$(AR) -x external/ldl/libldl.a

# build Tim Davis' AMD package
$(AMD):
	( cd external/amd    ; $(MAKE) )
	$(AR) -x external/amd/libamd.a

# build ECOS
ECOS_OBJS = ecos.o kkt.o cone.o spla.o ctrlc.o timer.o preproc.o splamm.o equil.o expcone.o wright_omega.o
libecos.a: $(ECOS_OBJS) $(LDL) $(AMD)
	$(ARCHIVE) $@ $^
	- $(RANLIB) $@

# build ECOS branch-and-bound
ECOS_BB_OBJS = $(ECOS_OBJS) ecos_bb_preproc.o ecos_bb.o
libecos_bb.a: $(ECOS_BB_OBJS) $(LDL) $(AMD)
	$(ARCHIVE) $@ $^
	- $(RANLIB) $@

%.o : src/%.c
	$(CC) $(CFLAGS) -c $< -o $@

%.o : ecos_bb/%.c
	$(CC) $(CFLAGS) -c $< -o $@

ecos_bb.o           : include/ecos_bb.h
ecos_bb_preproc.o   : include/ecos_bb.h
ecos.o              : include/ecos.h include/glblopts.h include/cone.h
kkt.o               : include/kkt.h include/glblopts.h include/cone.h include/ecos.h
cone.o              : include/cone.h include/glblopts.h include/ecos.h
preproc.o           : include/glblopts.h include/cone.h include/ecos.h
spla.o              : include/spla.h include/glblopts.h include/cone.h include/ecos.h
splamm.o            : include/splamm.h include/glblopts.h include/cone.h include/ecos.h
ctrlc.o             : include/ctrlc.h include/glblopts.h include/cone.h include/ecos.h
timer.o             : include/timer.h include/glblopts.h include/cone.h include/ecos.h
equil.o             : include/equil.h include/glblopts.h include/cone.h include/ecos.h
expcone.o           : include/expcone.h  
wright_omega.o      : include/wright_omega.h


# ECOS demo
.PHONY: demo
demo: runecos
runecos: src/runecos.c libecos.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)
	echo ECOS successfully built. Type ./runecos to run demo problem.

runecosexp: src/runecos_exp.c libecos.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)
	echo ECOS-Exp successfully built. Type ./runecosexp to run demo problem.

# Shared library
.PHONY: shared
shared: $(SHAREDNAME)
$(SHAREDNAME): $(LDL) $(AMD) $(ECOS_OBJS)
	$(CC) $(CFLAGS) -shared -o $@ $^ $(LDFLAGS)

# ECOS tester
TEST_OBJS = qcml_utils.o norm.o sq_norm.o sum_sq.o quad_over_lin.o inv_pos.o
.PHONY: test
test: ecostester ecos_bb_test
ecostester: test/ecostester.c $(TEST_OBJS) libecos.a
	$(CC) $(CFLAGS) $(TEST_INCLUDES) -o $@ $^ $(LDFLAGS)

ecos_bb_test: test/bb_test.c libecos_bb.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

%.o: test/generated/%.c test/generated/%.h
	$(CC) $(CFLAGS) $(TEST_INCLUDES) -c $< -o $@

%.o: test/generated/*/%.c test/generated/*/%.h
	$(CC) $(CFLAGS) $(TEST_INCLUDES) -c $< -o $@


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
	- $(RM) libecos.a libecos_bb.a runecos runecosexp
