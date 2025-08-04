#include <Grid/Grid.h>

struct hmc_params {
  int save_freq;
  int starttraj;
  double beta;
  double m;
  double tlen;
  int nsteps;
  std::string serial_seed = "1 2 3 4 5";
  std::string parallel_seed = "6 7 8 9 10";
  int Ls = 8;
  double M5 = 1.8; // domain wall height
  double b = 1.5; // controls exactly what action is being used
  double c = 0.5; // b-c=1.0 (fixed?); alpha=b+c=1 is Shamir (see e.g. 1411.7017)
  std::string cnfg_dir = ".";
};

hmc_params ReadCommandLineHMC(int argc, char** argv) {
  hmc_params HMCParams;
  if (Grid::GridCmdOptionExists(argv, argv + argc, "--savefreq")) {
    HMCParams.save_freq = std::stoi(Grid::GridCmdOptionPayload(argv, argv + argc, "--savefreq"));
  } else {
    std::cout << Grid::GridLogError << "--savefreq must be specified" << std::endl;
    exit(1);
  }
  if (Grid::GridCmdOptionExists(argv, argv + argc, "--beta")) {
    HMCParams.beta = std::stod(Grid::GridCmdOptionPayload(argv, argv + argc, "--beta"));
  } else {
    std::cout << Grid::GridLogError << "--beta must be specified" << std::endl;
    exit(1);
  }
    if (Grid::GridCmdOptionExists(argv, argv + argc, "--starttraj")) {
     HMCParams.starttraj = std::stod(Grid::GridCmdOptionPayload(argv, argv + argc, "--starttraj"));
  } else {
    std::cout << Grid::GridLogError << "--starttraj must be specified" << std::endl;
  }
  if (Grid::GridCmdOptionExists(argv, argv + argc, "--fermionmass")) {
    HMCParams.m = std::stod(Grid::GridCmdOptionPayload(argv, argv + argc, "--fermionmass"));
  } else {
    std::cout << Grid::GridLogError << "--fermionmass must be specified" << std::endl;
    exit(1);
  }
  if (Grid::GridCmdOptionExists(argv, argv + argc, "--tlen")) {
    HMCParams.tlen = std::stod(Grid::GridCmdOptionPayload(argv, argv + argc, "--tlen"));
  } else {
    std::cout << Grid::GridLogError << "--tlen must be specified" << std::endl;
    exit(1);
  }
  if (Grid::GridCmdOptionExists(argv, argv + argc, "--nsteps")) {
    HMCParams.nsteps = std::stoi(Grid::GridCmdOptionPayload(argv, argv + argc, "--nsteps"));
  } else {
    std::cout << Grid::GridLogError << "--nsteps must be specified" << std::endl;
    exit(1);
  }
  if (Grid::GridCmdOptionExists(argv, argv + argc, "--serialseed")) {
    HMCParams.serial_seed = Grid::GridCmdOptionPayload(argv, argv + argc, "--serialseed");
  }
  if (Grid::GridCmdOptionExists(argv, argv + argc, "--parallelseed")) {
    HMCParams.parallel_seed = Grid::GridCmdOptionPayload(argv, argv + argc, "--parallelseed");
  }
  if (Grid::GridCmdOptionExists(argv, argv + argc, "--cnfg_dir")) {
    HMCParams.cnfg_dir = Grid::GridCmdOptionPayload(argv, argv + argc, "--cnfg_dir");
  }

  return HMCParams;
}



int main(int argc, char **argv) {
  using namespace Grid;

  typedef Representations< SpFundamentalRepresentation > TheRepresentations;

  Grid_init(&argc, &argv);

  typedef GenericSpHMCRunnerHirep<TheRepresentations, MinimumNorm2> HMCWrapper; // ok
  typedef SpWilsonImplR FermionImplPolicy;                                    // ok
  typedef SpWilsonFermionD FermionAction;                                     // ok
  typedef typename FermionAction::FermionField FermionField;                  // ok?

  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  HMCWrapper TheHMC;
  hmc_params HMCParams = ReadCommandLineHMC(argc, argv);
  TheHMC.Resources.AddFourDimGrid("gauge");
  TheHMC.Parameters.MD.MDsteps = HMCParams.nsteps;
  TheHMC.Parameters.MD.trajL = HMCParams.tlen;

  TheHMC.Parameters.StartTrajectory = HMCParams.starttraj;

  CheckpointerParameters CPparams;
  CPparams.config_prefix = HMCParams.cnfg_dir + "/ckpoint_lat";
  CPparams.rng_prefix    = HMCParams.cnfg_dir + "/ckpoint_rng";
  CPparams.saveInterval  = HMCParams.save_freq;
  CPparams.format        = "IEEE64BIG";
  TheHMC.Resources.LoadNerscCheckpointer(CPparams);

  RNGModuleParameters RNGpar;
  RNGpar.serial_seeds = HMCParams.serial_seed;
  RNGpar.parallel_seeds = HMCParams.parallel_seed;
  TheHMC.Resources.SetRNGSeeds(RNGpar);

  // Construct observables
  typedef PlaquetteMod<HMCWrapper::ImplPolicy> PlaqObs;
  TheHMC.Resources.AddObservable<PlaqObs>();

  RealD beta = HMCParams.beta ;

  SpWilsonGaugeActionR Waction(beta);

  auto GridPtr = TheHMC.Resources.GetCartesian();
  auto GridRBPtr = TheHMC.Resources.GetRBCartesian();

  SpFundamentalRepresentation::LatticeField U(GridPtr);

  RealD mass = HMCParams.m;

  FermionAction FermOp(U, *GridPtr, *GridRBPtr, mass);

  ConjugateGradient<FermionField> CG(1.0e-8, 2000, false);

  TwoFlavourPseudoFermionAction<FermionImplPolicy> Nf2(FermOp, CG, CG);

  Nf2.is_smeared = false;

  ActionLevel<HMCWrapper::Field, TheRepresentations > Level1(1);
  Level1.push_back(&Nf2);

  ActionLevel<HMCWrapper::Field, TheRepresentations > Level2(4);
  Level2.push_back(&Waction);

  TheHMC.TheAction.push_back(Level1);
  TheHMC.TheAction.push_back(Level2);

  TheHMC.ReadCommandLine(argc, argv);
  TheHMC.Run();

  Grid_finalize();
}
