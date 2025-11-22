#!/usr/bin/bash
tput bold;
echo "! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !"
echo "!                                                                       !"
echo "!     Code to configure Grid-Telos                                      !"
echo "!     $scrfipt_file_name                                            !"
echo "!     [Author]: Frederic Bonnet November 2024                           !"
echo "!     [usage]: getConfigure-Grid-Telos.sh                               !"
echo "!     [example]: getConfigure-Grid-Telos-Devel.sh                       !"
echo "!                                                                       !"
echo "! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !"
tput sgr0;
#colors
red="tput setaf 1"  ;green="tput setaf 2"  ;yellow="tput setaf 3"
blue="tput setaf 4" ;magenta="tput setaf 5";cyan="tput setaf 6"
white="tput setaf 7";bold=""               ;reset_colors="tput sgr0"
##########################
# MAresNostrum setup
##########################
#module module load cuda/12.3 openmpi/4.1.5-cuda12.3 ucx/1.15.0-cuda12.3 gcc/9.3.0;


source /etc/profile.d/modules.sh ;
module load /mnt/lustre/tursafs1/home/y07/shared/tursa-modules/setup-env ;
module load cuda/12.3 openmpi/4.1.5-cuda12.3 ucx/1.15.0-cuda12.3 gcc/9.3.0;

module list;


# Start of the script
prefix="/home/dp208/dp208/dc-bonn2/SwanSea/SourceCodes/external_lib/prefix_grid_202410"

:'
../configure \
    --prefix=${prefix} \
    --enable-comms=mpi \
    --enable-simd=GPU \
    --enable-shm=nvlink \
    --enable-accelerator=cuda \
    --enable-gen-simd-width=64 \
    --disable-unified \
    --with-lime=${prefix} \
    --with-gmp=${prefix} \
    --with-mpfr=${prefix} \
    --disable-fermion-reps \
    --disable-gparity \
    --enable-Sp \
    --enable-Nc=4 \
    --enable-llr \
    CXX=nvcc \
    LDFLAGS="-cudart shared -lcublas" \
    CXXFLAGS="-ccbin mpicxx -gencode arch=compute_80,code=sm_80 -std=c++17 -cudart shared --diag-suppress 177,550,611"
'

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
--enable-Sp \
--enable-Nc=4 \
--with-lime=/home/dp208/dp208/dc-bonn2/SwanSea/SourceCodes/external_lib/prefix_grid_202410 \
--with-gmp=/home/dp208/dp208/dc-bonn2/SwanSea/SourceCodes/external_lib/prefix_grid_202410 \
--with-mpfr=/home/dp208/dp208/dc-bonn2/SwanSea/SourceCodes/external_lib/prefix_grid_202410 \
--disable-unified \
CXX=nvcc \
LDFLAGS="-cudart shared -lcublas"  \
CXXFLAGS="-ccbin mpicxx -gencode arch=compute_80,code=sm_80 -std=c++17 -cudart shared --diag-suppress 177,550,611"

#-------------------------------------------------------------------------------
#End of the script
echo
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
$cyan; echo `date`; $blue;
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
echo "-                  getConfigure-Grid-Telos-Devel.sh Done.               -"
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
$white; $reset_colors;
#exit
#-------------------------------------------------------------------------------
