#include <Grid/Grid.h>

int main(int argc, char **argv) {
  using namespace Grid;

  
  typedef Representations<SpFundamentalRepresentation,
                          SpTwoIndexAntiSymmetricRepresentation>
      TheRepresentations;
  
  Grid_init(&argc, &argv);
  typedef GenericSpHMCRunnerHirep<TheRepresentations, MinimumNorm2>
      HMCWrapper;
  typedef SpWilsonTwoIndexAntiSymmetricImplR FermionImplPolicy;
  typedef SpWilsonTwoIndexAntiSymmetricFermionD FermionAction;
  typedef typename FermionAction::FermionField FermionField;

  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  HMCWrapper TheHMC;

  TheHMC.Resources.AddFourDimGrid("gauge");

  // Checkpointer definition
  CheckpointerParameters CPparams;
  CPparams.config_prefix = "ckpoint_lat";
  CPparams.rng_prefix = "ckpoint_rng";
  CPparams.saveInterval = 100;
  CPparams.format = "IEEE64BIG";

  TheHMC.Resources.LoadNerscCheckpointer(CPparams);

  RNGModuleParameters RNGpar;
  RNGpar.serial_seeds = "1 2 3 4 5";
  RNGpar.parallel_seeds = "6 7 8 9 10";
  TheHMC.Resources.SetRNGSeeds(RNGpar);

  // Construct observables
  typedef PlaquetteMod<HMCWrapper::ImplPolicy> PlaqObs;
  TheHMC.Resources.AddObservable<PlaqObs>();

  RealD beta = 6.75;
  std::vector<RealD> hasenbusch({ 0.7 });
  SpWilsonGaugeActionR Waction(beta);


  auto GridPtr = TheHMC.Resources.GetCartesian();
  auto GridRBPtr = TheHMC.Resources.GetRBCartesian();

  SpTwoIndexAntiSymmetricRepresentation::LatticeField U(GridPtr);
 
  // LatticeGaugeField U(GridPtr);

  RealD mass = -1.046;

  std::vector<Complex> boundary = {-1, -1, -1, -1};
  FermionAction::ImplParams bc(boundary);
  FermionAction FermOp(U, *GridPtr, *GridRBPtr, mass, bc);

  ConjugateGradient<FermionField> CG(1.0e-8, 2000, false);

  ActionLevel<HMCWrapper::Field, TheRepresentations> Level1(1);
  ActionLevel<HMCWrapper::Field, TheRepresentations> Level2(4);

  std::cout << GridLogMessage << "mass " << mass << std::endl;

  // Hasenbusch
  std::vector<RealD> light_den;
  std::vector<RealD> light_num;

  int n_hasenbusch = hasenbusch.size();
  light_den.push_back(mass);
  for(int h=0;h<n_hasenbusch;h++){
    light_den.push_back(hasenbusch[h]);
    light_num.push_back(hasenbusch[h]);
  }

  std::vector<FermionAction *> Numerators;
  std::vector<FermionAction *> Denominators;
  std::vector<TwoFlavourEvenOddRatioPseudoFermionAction<FermionImplPolicy> *> Quotients;

  for(int h=0;h<n_hasenbusch+1;h++){
    std::cout << GridLogMessage << " 2f quotient Action  "<< light_num[h] << " / " << light_den[h]<< std::endl;
    Numerators.push_back  (new FermionAction(U, *GridPtr, *GridRBPtr, light_num[h], bc));
    Denominators.push_back(new FermionAction(U, *GridPtr, *GridRBPtr, light_den[h], bc));
    Quotients.push_back   (new TwoFlavourEvenOddRatioPseudoFermionAction<FermionImplPolicy>(*Numerators[h],*Denominators[h],CG,CG));
  }

  for(int h=0;h<n_hasenbusch+1;h++){
    Level1.push_back(Quotients[h]);
  }
  Level2.push_back(&Waction);
  
  TheHMC.TheAction.push_back(Level1);
  TheHMC.TheAction.push_back(Level2);

  TheHMC.Parameters.MD.MDsteps = 16;
  TheHMC.Parameters.MD.trajL = 1.0;

  TheHMC.ReadCommandLine(argc, argv);
  TheHMC.Run();

  Grid_finalize();
}
