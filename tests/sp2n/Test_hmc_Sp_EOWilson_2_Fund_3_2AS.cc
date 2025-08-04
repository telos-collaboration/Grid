#include <Grid/Grid.h>

struct hmc_params {
  std::string save_path;
  int starttraj;
  double fund_mass;
  double AS_mass;
  int nsteps;
  int savefreq;
};



hmc_params ReadCommandLineHMC(int argc, char** argv) {
  hmc_params HMCParams;
  if (Grid::GridCmdOptionExists(argv, argv + argc, "--savepath")) {
     HMCParams.save_path = Grid::GridCmdOptionPayload(argv, argv + argc, "--savepath");
  } else {
    std::cout << Grid::GridLogError << "--savepath must be specified" << std::endl;
  }
  if (Grid::GridCmdOptionExists(argv, argv + argc, "--starttraj")) {
     HMCParams.starttraj = std::stod(Grid::GridCmdOptionPayload(argv, argv + argc, "--starttraj"));
  } else {
    std::cout << Grid::GridLogError << "--starttraj must be specified" << std::endl;
  }
  if (Grid::GridCmdOptionExists(argv, argv + argc, "--fundmass")) {
     HMCParams.fund_mass = std::stod(Grid::GridCmdOptionPayload(argv, argv + argc, "--fundmass"));
  } else {
    std::cout << Grid::GridLogError << "--fundmass must be specified" << std::endl;
  }
  if (Grid::GridCmdOptionExists(argv, argv + argc, "--ASmass")) {
     HMCParams.AS_mass = std::stod(Grid::GridCmdOptionPayload(argv, argv + argc, "--ASmass"));
  } else {
    std::cout << Grid::GridLogError << "--ASmass must be specified" << std::endl;
  }
  if (Grid::GridCmdOptionExists(argv, argv + argc, "--nsteps")) {
     HMCParams.nsteps = std::stod(Grid::GridCmdOptionPayload(argv, argv + argc, "--nsteps"));
  } else {
    std::cout << Grid::GridLogError << "--nsteps must be specified" << std::endl;
  }
    if (Grid::GridCmdOptionExists(argv, argv + argc, "--savefreq")) {
     HMCParams.savefreq = std::stod(Grid::GridCmdOptionPayload(argv, argv + argc, "--savefreq"));
  } else {
    std::cout << Grid::GridLogError << "--savefreq must be specified" << std::endl;
  }


  return HMCParams;
}

int main(int argc, char **argv) {
  using namespace Grid;

  hmc_params HMCParams = ReadCommandLineHMC(argc, argv);

  typedef Representations< SpFundamentalRepresentation, SpTwoIndexAntiSymmetricRepresentation > TheRepresentations;

  Grid_init(&argc, &argv);

  typedef GenericSpHMCRunnerHirep<TheRepresentations, MinimumNorm2> HMCWrapper;

  typedef SpWilsonTwoIndexAntiSymmetricImplR TwoIndexFermionImplPolicy;
  typedef SpWilsonTwoIndexAntiSymmetricFermionD TwoIndexFermionAction;
  typedef typename TwoIndexFermionAction::FermionField TwoIndexFermionField;

  typedef SpWilsonImplR FundFermionImplPolicy;                                    // ok
  typedef SpWilsonFermionD FundFermionAction;                                     // ok
  typedef typename FundFermionAction::FermionField FundFermionField;

  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  HMCparameters HMCparams2;
  HMCparams2.StartTrajectory  = HMCParams.starttraj;

  HMCWrapper TheHMC(HMCparams2);

  TheHMC.Resources.AddFourDimGrid("gauge");

  // Checkpointer definition
  CheckpointerParameters CPparams;
  CPparams.config_prefix = "./" + HMCParams.save_path +"ckpoint_lat";
  CPparams.rng_prefix = "./" + HMCParams.save_path + "ckpoint_rng";
  CPparams.saveInterval = HMCParams.savefreq;
  CPparams.format = "IEEE64BIG";

  TheHMC.Resources.LoadNerscCheckpointer(CPparams);

  RNGModuleParameters RNGpar;
  RNGpar.serial_seeds = "1 2 3 4 5";
  RNGpar.parallel_seeds = "6 7 8 9 10";
  TheHMC.Resources.SetRNGSeeds(RNGpar);

  // Construct observables
  typedef PlaquetteMod<HMCWrapper::ImplPolicy> PlaqObs;
  TheHMC.Resources.AddObservable<PlaqObs>();

//      typedef PolyakovMod<HMCWrapper::ImplPolicy> PolyakovObs;
//      TheHMC.Resources.AddObservable<PolyakovObs>();



  RealD beta = 6.5 ;

  SpWilsonGaugeActionR Waction(beta);

  auto GridPtr = TheHMC.Resources.GetCartesian();
  auto GridRBPtr = TheHMC.Resources.GetRBCartesian();

  SpFundamentalRepresentation::LatticeField fundU(GridPtr);
  SpTwoIndexAntiSymmetricRepresentation::LatticeField asU(GridPtr);
    //LatticeGaugeField U(GridPtr);

  RealD Fundmass = HMCParams.fund_mass;
  RealD ASmass = HMCParams.AS_mass;
  std::vector<Complex> boundary = {-1,+1,+1,+1};

  FundFermionAction::ImplParams bc(boundary);
  TwoIndexFermionAction::ImplParams bbc(boundary);

  FundFermionAction FundFermOp(fundU, *GridPtr, *GridRBPtr, Fundmass, bbc);
  TwoIndexFermionAction TwoIndexFermOp(asU, *GridPtr, *GridRBPtr, ASmass, bbc);
  ConjugateGradient<FundFermionField> fCG(1.0e-8, 2000);
  ConjugateGradient<TwoIndexFermionField> asCG(1.0e-8, 2000);
  OneFlavourRationalParams Params(1.0e-6, 64.0, 2000, 1.0e-6, 16);


  TwoFlavourEvenOddPseudoFermionAction<FundFermionImplPolicy> fundNf2(FundFermOp, fCG, fCG);
  TwoFlavourEvenOddPseudoFermionAction<TwoIndexFermionImplPolicy> asNf2(TwoIndexFermOp, asCG, asCG);
  OneFlavourEvenOddRationalPseudoFermionAction<TwoIndexFermionImplPolicy> asNf1(TwoIndexFermOp,Params);

  fundNf2.is_smeared = false;
  asNf2.is_smeared = false;
  asNf1.is_smeared = false;

  ActionLevel<HMCWrapper::Field, TheRepresentations > Level1(1);
  Level1.push_back(&fundNf2);
  Level1.push_back(&asNf2);
  Level1.push_back(&asNf1);

  ActionLevel<HMCWrapper::Field, TheRepresentations > Level2(4);
  Level2.push_back(&Waction);

  TheHMC.TheAction.push_back(Level1);
  TheHMC.TheAction.push_back(Level2);

  TheHMC.Parameters.MD.MDsteps = HMCParams.nsteps;
  TheHMC.Parameters.MD.trajL   = 1.0;
  std::cout << "Fund_mass: " << Fundmass << "\n";
  std::cout << "AS_mass: " << ASmass << "\n";
  std::cout << "nsteps: " << HMCParams.nsteps << "\n";
  TheHMC.ReadCommandLine(argc, argv);
  TheHMC.Run();

  Grid_finalize();
}

