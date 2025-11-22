#!/bin/bash

../configure \
  --prefix=/home/dp208/dp208/dc-bonn2/SwanSea/SourceCodes/external_lib/prefix_grid_202410 \
  --enable-doxygen-doc \
  --enable-comms=mpi \
  --enable-simd=GPU \
  --enable-shm=nvlink \
  --enable-accelerator=cuda \
  --enable-gen-simd-width=64 \
  --disable-gparity \
  --disable-fermion-reps \
  --enable-llr \
  --with-lime=/home/dp208/dp208/dc-bonn2/SwanSea/SourceCodes/external_lib/prefix_grid_202410 \
  --with-gmp=/home/dp208/dp208/dc-bonn2/SwanSea/SourceCodes/external_lib/prefix_grid_202410 \
  --with-mpfr=/home/dp208/dp208/dc-bonn2/SwanSea/SourceCodes/external_lib/prefix_grid_202410 \
  --disable-unified \
  CXX=nvcc \
  LDFLAGS="-cudart shared -lcublas" \
  CXXFLAGS="-ccbin mpicxx -gencode arch=compute_80,code=sm_80 -std=c++17 -cudart shared --diag-suppress 177,550,611"



