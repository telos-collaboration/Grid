#!/bin/bash
#SBATCH --job-name=Run_Test_llr_contrained_action_nodes001_mpi01-01-01-04
#SBATCH --output=/home/dp208/dp208/dc-bonn2/SwanSea/SourceCodes/Grid-Telos-Devel/%x.out
#SBATCH --error=/home/dp208/dp208/dc-bonn2/SwanSea/SourceCodes/Grid-Telos-Devel/%x.err
#SBATCH --time=02:00:00
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1    # nodes * ntasks
## GPU only
#SBATCH --gpus-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --qos=dev
#SBATCH --account=dp208
#-------------------------------------------------------------------------------
# Getting the common code setup and variables: 
#---> Accelerator type (cpu/gpu)             : gpu
#---> Simulation size in consideration       : small
#---> Machine name that we are working on is : tursa
#-------------------------------------------------------------------------------
#---> this is a Grid_DWF_run_gpu job run
#-------------------------------------------------------------------------------
# Start of the batch body
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Module loads and compiler version
#-------------------------------------------------------------------------------
module load cuda/12.3 openmpi/4.1.5-cuda12.3 ucx/1.15.0-cuda12.3 gcc/9.3.0; module list;

# Check some versions
ucx_info -v
nvcc --version
which mpirun
#-------------------------------------------------------------------------------
# The path structure
#-------------------------------------------------------------------------------
machine_name="tursa"
sourcecode_dir=/home/dp208/dp208/dc-bonn2/SwanSea/SourceCodes
Bench_Grid_HiRep_dir=$sourcecode_dir/Bench_Grid_HiRep

# Application paths
grid_dwf_telos_dir=/home/dp208/dp208/dc-bonn2/SwanSea/SourceCodes/Grid-Telos-Devel/Grid
grid_dwf_telos_build_dir=$grid_dwf_telos_dir/build-Sp2n-withLLR
path_to_test_llr=$grid_dwf_telos_build_dir/tests/llr/
#Extending the library path
prefix=/home/dp208/dp208/dc-bonn2/SwanSea/SourceCodes/external_lib/prefix_grid_202410
#-------------------------------------------------------------------------------
# Export path and library paths
#-------------------------------------------------------------------------------
#Extending the library path
export PREFIX_HOME=$prefix
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PREFIX_HOME/lib

echo "$LD_LIBRARY_PATH"

ls -al "$PREFIX_HOME/lib"
#-------------------------------------------------------------------------------
# Probing the file systems and getting some info
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Compiling the testcode
#-------------------------------------------------------------------------------
#cd "${path_to_test_llr}";
#make check;
#-------------------------------------------------------------------------------
# Variable exports
#-------------------------------------------------------------------------------
# OpenMP
export OMP_NUM_THREADS=8
# MPI
export OMPI_MCA_btl=^uct,openib
export OMPI_MCA_pml=ucx
export OMPI_MCA_io=romio321
export OMPI_MCA_btl_openib_allow_ib=true
export OMPI_MCA_btl_openib_device_type=infiniband
export OMPI_MCA_btl_openib_if_exclude=mlx5_1,mlx5_2,mlx5_3
# UCX
export UCX_TLS=gdr_copy,rc,rc_x,sm,cuda_copy,cuda_ipc
export UCX_RNDV_THRESH=16384
export UCX_RNDV_SCHEME=put_zcopy
export UCX_IB_GPU_DIRECT_RDMA=yes
export UCX_MEMTYPE_CACHE=n
# GRID
export GRID_ALLOC_NCACHE_SMALL=16
export GRID_ALLOC_NCACHE_LARGE=2
export GRID_ALLOC_NCACHE_HUGE=0
#-------------------------------------------------------------------------------
# Output variable.
#-------------------------------------------------------------------------------
LatticeRuns_dir=${sourcecode_dir}/LatticeRuns
path_to_run=${sourcecode_dir}/Grid-Telos-Devel
job_name=Run_Test_llr_contrained_action_nodes001_mpi01-01-01-04
#-------------------------------------------------------------------------------
# Variable list for grid command line argument list
#-------------------------------------------------------------------------------
# Computed on the benchmark batch creation
MPI=1.1.1.4

# Extracted from the configuration

VOL=8.8.8.8
BETA=2.4         # 6.9   beta

# Hardcoded variables

TRAJECTORIES=40 #500 #20         #100000
THERMALIZATIONS=20 #10
MDsteps=40
trajL=1
SAVEFREQ=100 #10

# Extracting the checkpoint from the lattice data
STARTTRAJ=0
echo "#-------------------------------------------------------------------------------"
printf "STARTTRAJ              : "; printf '%s'"$STARTTRAJ"; printf "\n";
echo "#-------------------------------------------------------------------------------"
echo "starting from hot start no configurations to be read in."
echo "#-------------------------------------------------------------------------------"
#-------------------------------------------------------------------------------
# Wrapper scripts Getting the gpu select script
#-------------------------------------------------------------------------------
wrapper_script=${Bench_Grid_HiRep_dir}/doc/BKeeper/gpu-mpi-wrapper-new-Tursa.sh
chmod a+x ${wrapper_script}
#-------------------------------------------------------------------------------
# Launching mechanism
#-------------------------------------------------------------------------------
# run! #########################################################################
device_mem=23000
shm=8192

#  --enable-llr \

mpirun -np ${SLURM_NTASKS} \
"${grid_dwf_telos_build_dir}"/tests/llr/Test_llr_contrained_action \
  --beta ${BETA} \
  --starttraj ${STARTTRAJ} \
  --grid ${VOL} \
  --enable-llr \
  --Trajectories ${TRAJECTORIES} \
  --Thermalizations ${THERMALIZATIONS} \
  --nsteps ${MDsteps} \
  --tlen ${trajL} \
  --savefreq ${SAVEFREQ}

################################################################################
#  --savefreq ${SAVEFREQ} > ./dwf_trials_verybigR1/hmc_${SLURM_JOB_ID}.out
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Finishing up
#-------------------------------------------------------------------------------
#End of the script
echo
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
echo `date`
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
echo "- ${path_to_run}/${job_name}.sh Done. -"
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
# srun --account={account_name} --partition={partition} --time=00:30:00 --nodes=1 --gres=gpu:4 --pty bash
##SBATCH --ntasks-per-socket=4
##SBATCH --mem=494000
#export MPICH_GPU_SUPPORT_ENABLED=1
#export UCX_TLS=self,sm,rc,ud
#export OMPI_MCA_PML="ucx"
#export OMPI_MCA_osc="ucx"
#-------------------------------------------------------------------------------
