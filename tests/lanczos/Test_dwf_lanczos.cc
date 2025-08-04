#include <Grid/Grid.h>

using namespace std;
using namespace Grid;

void run(){
  typedef typename MobiusFermionD::FermionField FermionField;
  
  const int Ls=8;
  Coordinate latt_size  ({12,12,12,16});
  //Coordinate latt_size  ({8,8,8,8});
  //std::vector<vComplexD> Simd_vec({1,1,1,2});
  //std::vector<vComplexD> Mpi_vec({1,1,1,2});
  std::vector<int> Simd_vec({vComplexD::Nsimd(),1});
  std::vector<int> Mpi_vec ({1,1,1,4});
  //GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(latt_size, GridDefaultSimd(Nd,vComplexD::Nsimd()),GridDefaultMpi());
  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(latt_size, GridDefaultSimd(Nd,vComplexD::Nsimd()),Mpi_vec);
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);
//  printf("UGrid=%p UrbGrid=%p FGrid=%p FrbGrid=%p\n",UGrid,UrbGrid,FGrid,FrbGrid);
  
  GridCartesian* UGridF = SpaceTimeGrid::makeFourDimGrid(latt_size, GridDefaultSimd(Nd, vComplexD::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian* UrbGridF = SpaceTimeGrid::makeFourDimRedBlackGrid(UGridF);
  GridCartesian* FGridF = SpaceTimeGrid::makeFiveDimGrid(Ls, UGridF);
  GridRedBlackCartesian* FrbGridF = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, UGridF);


  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  GridParallelRNG          RNG5(FGridF);  RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG          RNG4(UGridF);  RNG4.SeedFixedIntegers(seeds4);
  GridParallelRNG          RNG5rb(FrbGridF);  RNG5.SeedFixedIntegers(seeds5);

  LatticeGaugeField Umu(UGrid); 
  SU<Nc>::HotConfiguration(RNG4, Umu);
  FieldMetaData header;
  //std::string file("./ckpoint_EODWF_lat.4");

  int precision32 = 0;
  int tworow      = 0;
  //NerscIO::writeConfiguration(Umu,file,tworow,precision32);
  //NerscIO::readConfiguration(Umu,header,file);

  LatticeGaugeFieldD UmuF(UGridF); 
  precisionChange(UmuF, Umu);
  //static MobiusFermionF* getAction(LatticeGaugeFieldF &Umu, GridCartesian* FGrid, GridRedBlackCartesian* FrbGrid, GridCartesian* UGrid, GridRedBlackCartesian* UrbGrid)
	  
  RealD mass=0.001;
  RealD M5=1.8;
  RealD mob_b=1.5;
  std::vector<Complex> boundary = {1,1,1,-1};
  MobiusFermionD::ImplParams Params(boundary);
  

  std::cout << GridLogMessage << "mass "<<mass<<std::endl;
  std::cout << GridLogMessage << "M5 "<<M5<<std::endl;
  std::cout << GridLogMessage << "mob_b "<<mob_b<<std::endl;
  MobiusFermionD action(UmuF,*FGridF,*FrbGridF,*UGridF,*UrbGridF,mass,M5,mob_b,mob_b-1.,Params);
  
  //Action *action = Setup<Action>::getAction(UmuF,FGridF,FrbGridF,UGridF,UrbGridF);
  //InvG5LinearOperator<MobiusFermionD,FermionField> HermOp(action, -2.);
   
  //MdagMLinearOperator<MobiusFermionD,FermionField> HermOp(action);
  //SchurDiagTwoOperator<MobiusFermionD,FermionField> HermOp(action);
  SchurDiagOneOperator<MobiusFermionD,FermionField> HermOp(action);
  //Gamma5HermitianLinearOperator<MobiusFermionD,FermionField> HermOp(action);
  /*
  RealD tol = 1e-10;
  int cg_MaxIt = 5000;

  ConjugateGradient<FermionField>          FineCG(tol,cg_MaxIt);
  FermionField f_src_e(FrbGridF); f_src_e=1.0;
  FermionField InvHermOp(FrbGridF); InvHermOp = Zero(); 
  FineCG(HermOp,f_src_e,InvHermOp);
  */
  const int Nstop = 10;
  const int Nk = 16;
  const int Np = 50;
  const int Nm = Nk+Np;
  const int MaxIt= 20000;
  RealD resid = 1.0e-6;
  std::cout << GridLogMessage << "Nstop "<<Nstop<<std::endl;
  std::cout << GridLogMessage << "Nk "<<Nk<<std::endl;
  std::cout << GridLogMessage << "Np "<<Np<<std::endl;
  std::cout << GridLogMessage << "resid "<<resid<<std::endl;

  std::vector<double> Coeffs { 0.,-1.};
  Polynomial<FermionField> PolyX(Coeffs);
  Chebyshev<FermionField> Cheby(0.001,4.5,301);
  //Chebyshev<FermionField> Cheby(0.0000006,5.5,4001);
  //Chebyshev<FermionField> Cheby(0.006,5.5,4001);
  //Chebyshev<FermionField> Cheby(5.5,0.006,4001);
  //Chebyshev<FermionField> Cheby(0.003,60,201);

  std::cout << GridLogMessage << "Cheby(0.0000006,5.5,4001) "<<std::endl;

  //Chebyshev<FermionField> Cheby2(3.00,7,50);

  FunctionHermOp<FermionField> OpCheby(Cheby, HermOp);
  //FunctionHermOp<FermionField> OpCheby2(Cheby2, HermOp);
  PlainHermOp<FermionField> Op     (HermOp);
  
  ImplicitlyRestartedLanczos<FermionField> IRL(OpCheby,Op,Nk,Nk,Nm,resid,MaxIt);
  //ImplicitlyRestartedLanczos<FermionField> IRL2(OpCheby,Op,Nk,Nk,Nm,resid,MaxIt); 
  std::vector<RealD>          eval(Nm);
  FermionField    src(FrbGridF); 
  gaussian(RNG5rb,src);
  std::vector<FermionField> evec(Nm,FrbGridF);
  for(int i=0;i<1;i++){
    std::cout << GridLogMessage <<i<<" / "<< Nm<< " grid pointer "<<evec[i].Grid()<<std::endl;
  };

  int Nconv;
  IRL.calc(eval,evec,src,Nconv);
  // Now eval[0] holds the lightest eigenvalue
  std::cout << GridLogMessage << "Smallest eigenvalue: " << eval[0] << std::endl;
  std::cout << eval << std::endl;
  
  std::cout << "ChebyIRL" << std::endl;
  /* 
  std::vector<RealD>          eval2(Nm);
  std::vector<FermionField> evec2(Nm,FrbGridF);
  int Nconv2;
  IRL2.calc(eval2,evec2,src,Nconv2);
  // Now eval[0] holds the lightest eigenvalue
  std::cout << GridLogMessage << "Smallest eigenvalue: " << eval2[0] << std::endl;
  std::cout << eval2 << std::endl;
  */
}
  
int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);
 
  run();
  
  Grid_finalize();
}
