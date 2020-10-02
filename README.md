# SQISign: compact post-quantum signatures from quaternions and isogenies

This code implements the isogeny-based signature scheme SQISign.

(C) 2020, The SQISign team. MIT license.

## Dependencies

The code depends on the latest stable version of the [PARI/GP
library](http://pari.math.u-bordeaux.fr/), 2.11.4.

The code has an optional dependency on [GMP](https://gmplib.org/),
which is also an optional dependency of PARI/GP and is typically
installed along with it.

## Supported platforms

The code compiles and runs on Linux and MacOS.

It contains two implementations of the low-level arithmetic functions:

- One based on handwritten assembly for the x86-64 platform,
  starting from Broadwell architectures.

- One based on GMP.

By default, both versions are compiled and tested.

## Compile

To compile and test the code, run

```
make
make check
```

The tests typically take 2-3 minutes.

To only compile and test the assembly version:

```
make asm
make check_asm
```

To only compile and test the GMP version:

```
make gmp
make check_gmp
```

## Tuning

The algorithms have already been tuned on a modern processors, and the
default values should be fine. However, should you wish to fine-tune
for your machine, run:

```
make tune
make
```

Be patient, tuning takes several minutes.

## Run benchmarks

To run benchmarks type

```
make benchmark
```

Allow a few minutes for the benchmarks to complete.  This produces
files named `bench_xxx.tsv` containing timing information on the
various parts of the signature.

By default, only benchmarks for the assembly version are produced. If
you want benchmarks for the GMP version, run

```
make benchmark_gmp
```

### Timings

The following are timings (medians) obtained running the benchmarks
above on an Intel Core i7-6700 CPU @ 3.40GHz with Turbo Boost
disabled.

|        | Mcycles |    ms |
|:-------|--------:|------:|
| keygen |   1,959 |   575 |
| sign   |   7,767 | 2,279 |
| verify |     142 |    42 |
