---
title : "Quick start guide"
authors:
- Guido Cossu
---

Please send all pull requests to the `develop` branch.

## Requirements

### Required libraries

* [GMP](https://gmplib.org/) is the GNU Multiple Precision Library (RHMC support).
* [MPFR](http://www.mpfr.org/) is a C library for multiple-precision floating-point computations with correct rounding (RHMC support). 
* [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page): bootstrapping GRID downloads and uses for internal dense matrix (non-QCD operations) the Eigen library.

Grid optionally uses:

* [HDF5](https://support.hdfgroup.org/HDF5/) for structured data I/O

* [LIME](http://usqcd-software.github.io/c-lime/) for ILDG and SciDAC file format support.

* [FFTW](http://www.fftw.org/) either generic version or via the Intel MKL library.

* [LAPACK](http://www.netlib.org/lapack/) either generic version or Intel MKL library.


### Compilers

Grid requires a C++ compiler with support for C++17.

There have historically been issues with certain GCC versions
(see e.g. issue [#100](https://github.com/paboyle/Grid/issues/100));
these are not expected to be observed
with current-generation GCC versions supporting C++17.

## Quick start
First, start by cloning the repository:

``` bash
git clone https://github.com/paboyle/Grid.git
```

Then enter the cloned directory and set up the build system:

```bash
cd Grid
./bootstrap.sh
```

Now you can execute the `configure` script to generate makefiles as in this example (here from a build directory):

``` bash
mkdir build
cd build
../configure --enable-precision=double --enable-simd=AVX --enable-comms=mpi-auto --prefix=<path>
```

where:

``` bash
  --enable-precision=single|double
```

sets the **default precision**. Since this is largely a benchmarking convenience, it is anticipated that the default precision may be removed in future implementations, and that explicit type selection be made at all points. Naturally, most code will be type templated in any case.::

``` bash
   --enable-simd=GEN|SSE4|AVX|AVXFMA|AVXFMA4|AVX2|AVX512|NEONv8|QPX
```

sets the **SIMD architecture**, 

``` bash
   --enable-comms=mpi|none
```

selects whether to use MPI communication (mpi) or no communication (none). 

``` bash
   --prefix=<path>
```

should be passed the prefix path where you want to install Grid. 

Other options are detailed in the next section, you can also use 

```bash
   configure --help
```

to display them. 

Like with any other program using GNU autotool, the 

```bash
   CXX, CXXFLAGS, LDFLAGS, ... 
```

environment variables can be modified to customise the build.

Finally, you can build and install Grid:

``` bash
make
make install   #this is optional
```

To minimise the build time, only the tests at the root of the `tests` directory are built by default. If you want to build tests in the sub-directory `<subdir>` you can execute:

``` bash
make -C tests/<subdir> tests
```
If you want to build all the tests at once just use `make tests`.

## Build configuration options

A full list of configurations options is available with the `./configure --help` command: 

- `--prefix=<path>`: installation prefix for Grid.
- `--enable-Nc={2|3|4|...}`:
  Set the default number of colours for lattice objects.
  For most datastructures this can be overridden in template parameters;
  for some components such as the HMC this is hardcoded and not overridable.
- `--with-gmp=<path>`: look for GMP in the UNIX prefix `<path>`
- `--with-mpfr=<path>`: look for MPFR in the UNIX prefix `<path>`
- `--with-fftw=<path>`: look for FFTW in the UNIX prefix `<path>`
- `--enable-lapack[=<path>]`:
  enable LAPACK support in Lanczos eigensolver.
  A UNIX prefix containing the library can be specified (optional).
- `--enable-mkl[=<path>]`:
  use Intel MKL for FFT (and LAPACK if enabled) routines.
  A UNIX prefix containing the library can be specified (optional).
- `--enable-simd=<code>`: set up Grid for the SIMD target `<code>` (default: `GEN`).
  A list of possible SIMD targets is detailed in
  [Architecture targets](simd-targets.md).
- `--enable-gen-simd-width=<size>`:
  select the size (in bytes) of the generic SIMD vector type
  (default: 64 bytes).
- `--enable-comms=<comm>`:
  Use `<comm>` for message passing (default: `none`).
  A list of possible communications targets is detailed in
  [Communication interfaces](communications-interfaces.md).
- `--enable-rng={sitmo|ranlux48|mt19937}`:
  choose the RNG (default: `sitmo `).
  See also [Random number generators](../api-description/random-number-generators.md)
- `--disable-timers`: disable system dependent high-resolution timers.
- `--enable-chroma`: enable Chroma regression tests.
- `--enable-doxygen-doc`:
  enable the Doxygen documentation generation
  (build with `make doxygen-doc`)
- `--disable-fermion-reps`, `--disable-fermion-instantiations`:
  disable building instantiations of higher fermion representations
  (two-index symmetric antisymmetric, adjoint),
  or all fermion representations
  (including the fundamental)
  respectively,
  in order to reduce the compilation time.

More details on the *Getting started* menu entries on the left. 

This document was updated in March 2026.
