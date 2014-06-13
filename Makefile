# Makefile for ECOS
# Configuration of make process in ecos.mk
# ECOS JNILib
PACKAGE=com.verizon.jecos
PACKAGE_PATH=$(subst .,/,$(PACKAGE))

SRC=src/main
SRC_JAVA=$(SRC)/java
SRC_C=$(SRC)/native
RESOURCES=$(SRC)/main/resources

TARGET_C=target/
LIB_PATH=$(RESOURCES)/lib

GENERATED_HEADERS=include/com_verizon_jecos_NativeECOS.h
#ant javah generates the NativeECOS header, implement the NativeECOS.c driver
GENERATED_SOURCES=${SRC_C}/NativeECOS.c

include ecos.mk
#jni headers specific for MacOSX 10.9
C = $(CC) $(CFLAGS) -I/System/Library/Frameworks/JavaVM.framework/Headers/ -Iinclude -Iexternal/ldl/include -Iexternal/amd/include -Iexternal/SuiteSparse_config 
TEST_INCLUDES = -Itest -Itest/quadratic

# Compile all C code, including the C-callable routine
all: ldl amd ecos demo

# build Tim Davis' sparse LDL package
ldl:
	( cd external/ldl    ; $(MAKE) )

# build Tim Davis' AMD package
amd:
	( cd external/amd    ; $(MAKE) )

# build ECOS
ecos: jniecos.o ecos.o kkt.o cone.o spla.o timer.o preproc.o splamm.o equil.o
	$(C) -L./external/amd -L./external/ldl -shared -o libecos.so jniecos.o ecos.o kkt.o cone.o spla.o timer.o preproc.o splamm.o equil.o -lldl -lamd
	cp libecos.so libecos.jnilib

jniecos.o: ${GENERATED_SOURCES} ${GENERATED_HEADERS}
	$(C) -fPIC -c ${GENERATED_SOURCES} -o jniecos.o

ecos.o: ${SRC_C}/ecos.c include/ecos.h
	$(C) -fPIC -c ${SRC_C}/ecos.c -o ecos.o

kkt.o: ${SRC_C}/kkt.c include/kkt.h
	$(C) -fPIC -c ${SRC_C}/kkt.c -o kkt.o

cone.o: ${SRC_C}/cone.c include/cone.h
	$(C) -fPIC -c ${SRC_C}/cone.c -o cone.o

preproc.o: ${SRC_C}/preproc.c
	$(C) -fPIC -c ${SRC_C}/preproc.c -o preproc.o

spla.o: ${SRC_C}/spla.c include/spla.h
	$(C) -fPIC -c ${SRC_C}/spla.c -o spla.o

splamm.o: ${SRC_C}/splamm.c include/splamm.h
	$(C) -fPIC -c ${SRC_C}/splamm.c -o splamm.o

timer.o: ${SRC_C}/timer.c include/timer.h
	$(C) -fPIC -c ${SRC_C}/timer.c -o timer.o

equil.o: ${SRC_C}/equil.c include/equil.h
	$(C) -fPIC -c ${SRC_C}/equil.c -o equil.o

# ECOS demo
demo: ldl amd ecos ${SRC_C}/runecos.c 
	$(C) -L. -L./external/amd -L./external/ldl -o runecos ${SRC_C}/runecos.c -lecos -lamd -ldl $(LIBS)
	echo ECOS successfully built. Type ./runecos to run demo problem.

# ECOS tester
TEST_OBJS = qcml_utils.o norm.o sq_norm.o sum_sq.o quad_over_lin.o inv_pos.o sqrt.o
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

sqrt.o: test/generated_tests/sqrt/sqrt.c test/generated_tests/sqrt/sqrt.h
	$(C) $(TEST_INCLUDES) -c test/generated_tests/sqrt/sqrt.c -o $@

# remove object files, but keep the compiled programs and library archives
clean:
	( cd external/ldl    ; $(MAKE) clean )
	( cd external/amd    ; $(MAKE) clean )
	- $(RM) $(CLEAN)

# clean, and then remove compiled programs and library archives
purge: clean
	( cd external/ldl    ; $(MAKE) purge )
	( cd external/amd    ; $(MAKE) purge )	
	- $(RM) libecos.so libecos.jnilib runecos
