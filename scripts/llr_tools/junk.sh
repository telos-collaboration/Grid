

        //LLRGaugeActionR *p_LLRaction_o =
        //        new Grid::LLRGaugeAction<Grid::SpWilsonGaugeActionR, Grid::PeriodicGimplR>(s_llrparams_in, beta);
        //p_LLRaction_o->Sa(&Grid::ConfigurationBase<Field>.)



../configure \
    --prefix=${prefix} \
    --enable-comms=mpi-auto \
    --enable-unified=no \
    --enable-shm=nvlink \
    --enable-accelerator=cuda \
    --enable-gen-simd-width=64 \
    --enable-simd=GPU \
    --enable-accelerator-cshift \
    --with-lime=${prefix} \
    --with-gmp=${prefix} \
    --with-mpfr=${prefix} \
    --disable-fermion-reps \
    --disable-gparity \
    --enable-Sp \
    --enable-Nc=4 \
    --enable-LLR \
    CXX=nvcc \
    LDFLAGS="-cudart shared -lcublas" \
    CXXFLAGS="-ccbin mpicxx -gencode arch=compute_90,code=sm_90 -std=c++17 -cudart shared --diag-suppress 177,550,611";


    --no-create \
    --no-recursion \



:'
mpirun -np ${SLURM_NTASKS} \
  --map-by numa \
  -x LD_LIBRARY_PATH \
  --bind-to none \
  "$wrapper_script" "${grid_dwf_telos_build_dir}"/HMC/LLR_HMC_WilsonGauge  \
  --StartingType HotStart \
  --starttraj $STARTTRAJ \
  --beta ${BETA} \
  --tlen ${TLEN} \
  --grid ${VOL} \
  --dwf_mass ${DWF_MASS} \
  --mobius_b ${MOBIUS_B} \
  --mobius_c ${MOBIUS_C} \
  --Ls ${Ls} \
  --savefreq ${SAVEFREQ} \
  --fermionmass ${MASS} \
  --nsteps ${NSTEPS} \
  --mpi ${MPI} \
  --shm ${shm} \
  --device-mem ${device_mem} \
  --cnfg_dir "./dwf_trials_verybigR1" \
  --accelerator-threads 8 \
  --Trajectories ${TRAJECTORIES} \
  --Thermalizations 10000 \
  --savefreq ${SAVEFREQ}
# > ./hmc_${SLURM_JOB_ID}.out
'
