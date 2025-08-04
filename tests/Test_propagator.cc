#include <Grid/Grid.h>
using namespace Grid;

void testDwfMatrix_momentum(GridCartesian &GRID, GridRedBlackCartesian &RBGRID, SpFundamentalRepresentation::LatticeField &Umu, const std::vector<int>& seeds, const int vol, const Coordinate latt_size, GridSerialRNG sRNG, GridParallelRNG pRNG)
{
  ////////////////////////////////////////////////////
  // Dwf matrix
  ////////////////////////////////////////////////////
    std::cout << "****************************************" << std::endl;
    std::cout << "Testing Fourier representation of Ddwf" << std::endl;
    std::cout << "****************************************" << std::endl;

    const int Ls = 16;
    const int sdir = 0;
    RealD mass = 0.01;
    RealD M5 = 1.0;
    Gamma G5(Gamma::Algebra::Gamma5);

    GridCartesian *FGrid = SpaceTimeGrid::makeFiveDimGrid(Ls, &GRID);
    GridRedBlackCartesian *FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, &GRID);
    
    typedef DomainWallFermion<SpWilsonImplD> FermionAction;
    
    std::cout << "Making Ddwf" << std::endl;
    FermionAction Ddwf(Umu, *FGrid, *FrbGrid, GRID, RBGRID, mass, M5);
    
    typedef typename FermionAction::FermionField SpLatticeFermionD;

    GridParallelRNG RNG5(FGrid);
    RNG5.SeedFixedIntegers(seeds);
    SpLatticeFermionD src5(FGrid);
    gaussian(RNG5, src5);
    SpLatticeFermionD src5_p(FGrid);
    LatticeFermionD result5(FGrid);
    SpLatticeFermionD ref5(FGrid);
    LatticeFermionD tmp5(FGrid);
    
    /////////////////////////////////////////////////////////////////
    // result5 is the non pert operator in 4d mom space
    /////////////////////////////////////////////////////////////////
    Ddwf.M(src5, tmp5);
    ref5 = tmp5;
    /*
    FFT theFFT5(FGrid);
    
    theFFT5.FFT_dim(result5, tmp5, 1, FFT::forward);
    tmp5 = result5;
    
    theFFT5.FFT_dim(result5, tmp5, 2, FFT::forward);
    tmp5 = result5;
    theFFT5.FFT_dim(result5, tmp5, 3, FFT::forward);
    tmp5 = result5;
    theFFT5.FFT_dim(result5, tmp5, 4, FFT::forward);
    result5 = result5 * ComplexD(::sqrt(1.0 / vol), 0.0);

    std::cout << "Fourier xformed Ddwf" << std::endl;

    tmp5 = src5;
    theFFT5.FFT_dim(src5_p, tmp5, 1, FFT::forward);
    tmp5 = src5_p;
    theFFT5.FFT_dim(src5_p, tmp5, 2, FFT::forward);
    tmp5 = src5_p;
    theFFT5.FFT_dim(src5_p, tmp5, 3, FFT::forward);
    tmp5 = src5_p;
    theFFT5.FFT_dim(src5_p, tmp5, 4, FFT::forward);
    src5_p = src5_p * ComplexD(::sqrt(1.0 / vol), 0.0);


    std::cout << "Fourier xformed src5" << std::endl;

    /////////////////////////////////////////////////////////////////
    // work out the predicted from Fourier
    /////////////////////////////////////////////////////////////////
    Gamma::Algebra Gmu[] = {
        Gamma::Algebra::GammaX,
        Gamma::Algebra::GammaY,
        Gamma::Algebra::GammaZ,
        Gamma::Algebra::GammaT,
        Gamma::Algebra::Gamma5};

    SpLatticeFermionD Kinetic(FGrid);
    Kinetic = Zero();
    LatticeComplexD kmu(FGrid);
    LatticeInteger scoor(FGrid);
    LatticeComplexD sk(FGrid);
    sk = Zero();
    LatticeComplexD sk2(FGrid);
    sk2 = Zero();
    LatticeComplexD W(FGrid);
    W = Zero();
    LatticeComplexD one(FGrid);
    one = ComplexD(1.0, 0.0);
    ComplexD ci(0.0, 1.0);

    for (int mu = 0; mu < Nd; mu++)
    {
        LatticeCoordinate(kmu, mu + 1);

        RealD TwoPiL = M_PI * 2.0 / latt_size[mu];

        kmu = TwoPiL * kmu;

        sk2 = sk2 + 2.0 * sin(kmu * 0.5) * sin(kmu * 0.5);
        sk = sk + sin(kmu) * sin(kmu);

        Kinetic = Kinetic + sin(kmu) * ci * (Gamma(Gmu[mu]) * src5_p);
    }

    W = one - M5 + sk2;
    Kinetic = Kinetic + W * src5_p;

    LatticeCoordinate(scoor, sdir);

    tmp5 = Cshift(src5_p, sdir, +1);
    tmp5 = (tmp5 - G5 * tmp5) * 0.5;
    tmp5 = where(scoor == Integer(Ls - 1), mass * tmp5, -tmp5);
    Kinetic = Kinetic + tmp5;

    tmp5 = Cshift(src5_p, sdir, -1);
    tmp5 = (tmp5 + G5 * tmp5) * 0.5;
    tmp5 = where(scoor == Integer(0), mass * tmp5, -tmp5);
    Kinetic = Kinetic + tmp5;

    std::cout << "Momentum space Ddwf  " << norm2(Kinetic) << std::endl;
    std::cout << "Stencil Ddwf         " << norm2(result5) << std::endl;

    result5 = result5 - Kinetic;
    std::cout << "diff " << norm2(result5) << std::endl;
    //assert(norm2(result5) < 1.0e-4);
       
    */
}


/*
void testDwfMatrix_space(GridCartesian &GRID, GridRedBlackCartesian &RBGRID, SpFundamentalRepresentation::LatticeField &Umu, const std::vector<int>& seeds, const int vol, const Coordinate latt_size, GridSerialRNG sRNG, GridParallelRNG pRNG)
{
  
    ////////////////////////////////////////////////////
  // Dwf prop
  ////////////////////////////////////////////////////
  {
    std::cout<<"****************************************"<<std::endl;
    std::cout << "Testing Ddwf Ht Mom space 4d propagator \n";
    std::cout<<"****************************************"<<std::endl;
    typedef DomainWallFermion<SpWilsonImplD> FermionAction;
    typedef typename FermionAction::FermionField SpLatticeFermionD;
    SpLatticeFermionD    src(&GRID); gaussian(pRNG,src);
    SpLatticeFermionD    tmp(&GRID);
    SpLatticeFermionD    ref(&GRID);
    SpLatticeFermionD    diff(&GRID);
    
    Coordinate point(4,0);
    src=Zero();
    SpinColourVectorD ferm; gaussian(sRNG,ferm);
    pokeSite(ferm,src,point);

    const int Ls=32;
    GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,&GRID);
    GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,&GRID);

    RealD mass=0.01;
    RealD M5  =0.8;
    DomainWallFermion<SpWilsonImplD> Ddwf(Umu,*FGrid,*FrbGrid,GRID,RBGRID,mass,M5);

    // Momentum space prop
    std::cout << " Solving by FFT and Feynman rules" <<std::endl;
    bool fiveD = false; //calculate 4d free propagator
    Ddwf.FreePropagator(src,ref,mass,fiveD) ;

    Gamma G5(Gamma::Algebra::Gamma5);

    SpLatticeFermionD    src5(FGrid); src5=Zero();
    SpLatticeFermionD    tmp5(FGrid); 
    SpLatticeFermionD    result5(FGrid); result5=Zero();
    SpLatticeFermionD    result4(&GRID); 
    const int sdir=0;

    ////////////////////////////////////////////////////////////////////////
    // Domain wall physical field source
    ////////////////////////////////////////////////////////////////////////    
    tmp =   (src + G5*src)*0.5;      InsertSlice(tmp,src5,   0,sdir);
    tmp =   (src - G5*src)*0.5;      InsertSlice(tmp,src5,Ls-1,sdir);
    
    ////////////////////////////////////////////////////////////////////////
    // Conjugate gradient on normal equations system
    ////////////////////////////////////////////////////////////////////////
    std::cout << " Solving by Conjugate Gradient (CGNE)" <<std::endl;
    Ddwf.Mdag(src5,tmp5);
    src5=tmp5;
    MdagMLinearOperator<DomainWallFermion<SpWilsonImplD>,SpLatticeFermionD> HermOp(Ddwf);
    ConjugateGradient<SpLatticeFermionD> CG(1.0e-16,10000);
    CG(HermOp,src5,result5);
    
    ////////////////////////////////////////////////////////////////////////
    // Domain wall physical field propagator
    ////////////////////////////////////////////////////////////////////////

    
    ExtractSlice(tmp,result5,0   ,sdir);   result4 =         (tmp-G5*tmp)*0.5;
    ExtractSlice(tmp,result5,Ls-1,sdir);   result4 = result4+(tmp+G5*tmp)*0.5;
    
    std::cout << " Taking difference" <<std::endl;
    std::cout << "Ddwf result4 "<<norm2(result4)<<std::endl;
    std::cout << "Ddwf ref     "<<norm2(ref)<<std::endl;
    
    diff = ref - result4;
    std::cout << "result - ref     "<<norm2(diff)<<std::endl;
    //assert(norm2(diff)<1.0e-4);

  }
      
}
*/

int main(int argc, char **argv)
{
    Grid_init(&argc, &argv);

    int threads = GridThread::GetThreads();
    std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;

    Coordinate latt_size = GridDefaultLatt();
    Coordinate simd_layout = GridDefaultSimd(Nd, vComplexD::Nsimd());
    Coordinate mpi_layout = GridDefaultMpi();

    int vol = 1;
    for (int d = 0; d < latt_size.size(); d++)
    {
        vol = vol * latt_size[d];
    }

    GridCartesian GRID(latt_size, simd_layout, mpi_layout);
    GridRedBlackCartesian RBGRID(&GRID);

    std::vector<int> seeds({1, 2, 3, 4});
    GridSerialRNG sRNG;
    sRNG.SeedFixedIntegers(seeds); // naughty seeding
    GridParallelRNG pRNG(&GRID);
    pRNG.SeedFixedIntegers(seeds);

    SpFundamentalRepresentation::LatticeField Umu(&GRID);
    //LatticeGaugeFieldD Umu(&GRID);
    SpFundamentalRepresentation::LatticeField Urnd(&GRID);
    //LatticeGaugeField Urnd(&GRID);
    LatticeColourMatrix xform(&GRID); // Gauge xform

    //LatticeGaugeField   Utmp(&GRID);
    SpFundamentalRepresentation::LatticeField Utmp(&GRID);
    Sp<Nc>::ColdConfiguration(pRNG, Umu); // Unit gauge
    Urnd = Umu;
    Sp<Nc>::RandomGaugeTransform<PeriodicGimplR>(pRNG, Urnd, xform);
    
    //Gauge fix the randomly transformed field 
    RealD alpha=0.1; //step size
    FourierAcceleratedGaugeFixer<PeriodicGimplR>::SteepestDescentGaugeFix(Urnd,xform,alpha,10000,1.0e-22, 1.0e-22,false);
    
    Utmp = Urnd - Umu;
    std::cout << " Check the output gauge transformation matrices applied to the original field produce the xformed field "<< norm2(Utmp) << " (expect 0)" << std::endl;
    
    // Call the new function here
    testDwfMatrix_momentum(GRID, RBGRID, Umu, seeds, vol, latt_size, sRNG, pRNG);
    std::cout << "\n\n\n\n" << std::endl;
    /*
    testDwfMatrix_space(GRID, RBGRID, Umu, seeds, vol, latt_size, sRNG, pRNG);
    std::cout << "\n\n\n\n" << std::endl;
    testDwfMatrix_space(GRID, RBGRID, Urnd, seeds, vol, latt_size, sRNG, pRNG);
     std::cout << "\n\n\n\n" << std::endl;
    testDwfMatrix_momentum(GRID, RBGRID, Urnd, seeds, vol, latt_size, sRNG, pRNG);
    */
    Grid_finalize();
    return 0;
}

