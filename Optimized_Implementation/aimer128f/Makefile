# SPDX-License-Identifier: MIT

CC ?= gcc
CFLAGS += -I. -O3 -Wall -Wextra -march=native -fomit-frame-pointer
NISTFLAGS = -Wno-sign-compare -Wno-unused-but-set-variable -Wno-unused-parameter -Wno-unused-result
AVX2FLAGS = -mavx2 -mpclmul

EXECUTABLE_TESTAIM2   = tests/test_aim2
EXECUTABLE_DEVEL      = tests/devel
EXECUTABLE_TESTSIGN   = tests/test_sign
EXECUTABLE_TESTSPEED  = tests/test_speed
EXECUTABLE_MICROBENCH = bench/microbench
EXECUTABLE_KAT        = PQCgenKAT_sign

RNG   = common/rng.o
SHA3  = common/aes.o common/keccak_avx2/fips202.o common/keccak_avx2/KeccakP-1600-AVX2.s common/keccak4x/KeccakP-1600-times4-SIMD256.o
AIMER = aim2.o field.o hash.o sign.o tree.o
BENCH = bench/bench.o bench/bench_aim2.o bench/bench_field.o

.PHONY: all

all: $(EXECUTABLE_TESTAIM2) $(EXECUTABLE_TESTSIGN) $(EXECUTABLE_TESTSPEED) $(EXECUTABLE_KAT) $(EXECUTABLE_MICROBENCH)

dev: $(EXECUTABLE_DEVEL)

.c.o:
	$(CC) -c $(CFLAGS) $(AVX2FLAGS) $< -o $@

$(EXECUTABLE_DEVEL): $(EXECUTABLE_DEVEL).c aim2.c field.c hash.c common/keccak_avx2/fips202.c common/keccak_avx2/KeccakP-1600-AVX2.s common/keccak4x/KeccakP-1600-times4-SIMD256.c
	$(CC) $(CFLAGS) $(AVX2FLAGS) $^ -o $@

$(EXECUTABLE_TESTAIM2): $(EXECUTABLE_TESTAIM2).o $(AIMER) $(SHA3) $(RNG)
	$(CC) $(CFLAGS) $(AVX2FLAGS) $^ -o $@

$(EXECUTABLE_TESTSIGN): $(EXECUTABLE_TESTSIGN).o $(AIMER) $(SHA3) $(RNG)
	$(CC) $(CFLAGS) $(AVX2FLAGS) $^ -o $@

$(EXECUTABLE_TESTSPEED): $(EXECUTABLE_TESTSPEED).o $(AIMER) $(SHA3) $(RNG) common/speed_print.o common/cpucycles.o
	$(CC) $(CFLAGS) $(AVX2FLAGS) $^ -o $@

$(EXECUTABLE_KAT): common/$(EXECUTABLE_KAT).o $(AIMER) $(SHA3) $(RNG)
	$(CC) $(CFLAGS) $(AVX2FLAGS) $^ -o $@

$(EXECUTABLE_MICROBENCH): $(EXECUTABLE_MICROBENCH).o $(AIMER) $(SHA3) $(RNG) $(BENCH)
	$(CC) $(CFLAGS) $(AVX2FLAGS) $^ -o $@

clean:
	rm -f $(wildcard *.o) $(wildcard tests/*.o) $(BENCH) $(EXECUTABLE_TESTAIM2) $(EXECUTABLE_TESTSIGN) $(EXECUTABLE_TESTSPEED) $(EXECUTABLE_KAT) $(EXECUTABLE_MICROBENCH) PQCsignKAT_*.rsp PQCsignKAT_*.req
