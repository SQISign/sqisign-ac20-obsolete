CC=clang
CFLAGS=-Wall -Wextra -pedantic -std=gnu99 -I./include -I/usr/local/include
CFLAGSLINK=-lpari -lm -L/usr/local/lib/ 
DEBUG_FLAGS=-g
BENCH_FLAGS=-DNDEBUG -O3 -Os -march=native -mtune=native

default: all

# "Library" object files depended upon by all executables
SRC6983    =fp2.c constants.c precomputed.c steps_tunecycles.c
SRC6983_ASM=fp.s 
SRC6983_GMP=fp.c
SRC    =toolbox.c ideal.c klpt.c idiso.c rng.c poly.c mont.c tedwards.c steps.c\
	isogenies.c isomorphism.c two_walks.c mitm.c sqisign.c #verif.c
SRC_ASM=uint.s
SRC_GMP=uint.c

OBJ6983    =$(SRC6983:%.c=build/obj/p6983/%.o)
OBJ6983_ASM=$(SRC6983_ASM:%.s=build/obj/p6983/%.o)
OBJ6983_GMP=$(SRC6983_GMP:%.c=build/obj/p6983/%_gmp.o)
OBJ        =$(SRC:%.c=build/obj/%.o)
OBJ_ASM    =$(SRC_ASM:%.s=build/obj/%.o)
OBJ_GMP    =$(SRC_GMP:%.c=build/obj/%_gmp.o)

LIB    =$(OBJ6983) $(OBJ) $(OBJ6983_ASM) $(OBJ_ASM)
LIB_GMP=$(OBJ6983) $(OBJ) $(OBJ6983_GMP) $(OBJ_GMP)

# Benchmarks
BENCHS_GMP=$(patsubst bench/%.c,build/bench_%_gmp,$(wildcard bench/*.c))
BENCHS=$(patsubst bench/%.c,build/bench_%,$(wildcard bench/*.c))

$(BENCHS): build/bench_%: bench/%.c $(SRC:%=src/%) $(SRC_ASM:%=src/%)\
$(SRC6983:%=src/p6983/%) $(SRC6983_ASM:%=src/p6983/%)
	@mkdir -p $(@D)
	$(CC) $^ $(CFLAGS) $(BENCH_FLAGS) $(CFLAGSLINK) -o $@

$(BENCHS_GMP): build/bench_%_gmp: bench/%.c $(SRC:%=src/%) $(SRC_GMP:%=src/%)\
$(SRC6983:%=src/p6983/%) $(SRC6983_GMP:%=src/p6983/%)
	@mkdir -p $(@D)
	$(CC) $^ $(CFLAGS) $(BENCH_FLAGS) $(CFLAGSLINK) -lgmp -o $@

# Tests
TESTS_GMP=$(patsubst test/%.c,build/test_%_gmp,$(wildcard test/*.c))
TESTS=$(patsubst test/%.c,build/test_%,$(wildcard test/*.c))

$(TESTS): build/test_%: test/%.c $(LIB)
	$(CC) $< $(LIB) $(CFLAGS) $(DEBUG_FLAGS) $(CFLAGSLINK) -o $@

$(TESTS_GMP): build/test_%_gmp: test/%.c $(LIB_GMP)
	$(CC) $< $(LIB_GMP) $(CFLAGS) $(DEBUG_FLAGS) $(CFLAGSLINK) -lgmp -o $@

# Additional executables
EXES=build/precomp

build/precomp: src/precomp.c $(LIB_GMP)
	$(CC) $< $(LIB_GMP) $(CFLAGS) $(DEBUG_FLAGS) $(CFLAGSLINK) -lgmp -o $@

# Velusqrt Tuning

build/tunecycles_6983: src/tunecycles.c src/isogenies.c src/mont.c src/p6983/fp2.c src/uint.s\
src/p6983/fp.s src/p6983/constants.c src/rng.c src/poly.c src/steps.c src/steps_default.c
	@mkdir -p $(@D)
	$(CC) $^ $(CFLAGS) $(BENCH_FLAGS) $(CFLAGSLINK) -o $@

src/p6983/tunecycles.out: build/tunecycles_6983
	# 8 minutes on 1.9GHz Kaby Lake
	time ./$< > $@

tune: src/tune2c src/p6983/tunecycles.out
	./src/tune2c < src/p6983/tunecycles.out > src/p6983/steps_tunecycles.c

# Object files
$(OBJ6983) $(OBJ): build/obj/%.o: src/%.c
	@mkdir -p $(@D)
	$(CC) $< $(CFLAGS) $(DEBUG_FLAGS) -c -o $@
$(OBJ6983_ASM) $(OBJ_ASM): build/obj/%.o: src/%.s
	@mkdir -p $(@D)
	$(CC) $< -c -o $@
$(OBJ6983_GMP) $(OBJ_GMP): build/obj/%_gmp.o: src/%.c
	@mkdir -p $(@D)
	$(CC) $< $(CFLAGS) $(DEBUG_FLAGS) -c -o $@

# Run tests
$(TESTS:%=%_run): %_run: %
	@echo
	./$^

$(TESTS_GMP:%=%_run): %_run: %
	@echo
	./$^

check_asm: $(TESTS:%=%_run)
check_gmp: $(TESTS_GMP:%=%_run)
check: check_asm check_gmp

# Run benchmarks
$(BENCHS:build/%=%.tsv): %.tsv: build/%
	./$^ >> $@

$(BENCHS_GMP:build/%=%.tsv): %.tsv: build/%
	./$^ >> $@

benchmark_asm: $(BENCHS:build/%=%.tsv)
benchmark_gmp: $(BENCHS_GMP:build/%=%.tsv)
benchmark: benchmark_asm

# Phony targets
benchs: $(BENCHS) $(BENCHS_GMP)
tests: $(TESTS) $(TESTS_GMP)
asm: $(BENCHS) $(TESTS)
gmp: $(BENCHS_GMP) $(TESTS_GMP)
all: $(EXES) asm gmp

distclean:
	rm -r build

.PHONY: distclean all gmp asm tests benchs benchmark benchmark_asm benchmark_gmp\
check check_gmp check_asm tune default
