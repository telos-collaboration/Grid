/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/Test_hmc_EODWFRatio.cc

Copyright (C) 2015-2016

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
Author: Guido Cossu <guido.cossu@ed.ac.uk>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */
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
  if (Grid::GridCmdOptionExists(argv, argv + argc, "--dw_mass")) {
    HMCParams.M5 = std::stod(Grid::GridCmdOptionPayload(argv, argv + argc, "--dw_mass"));
  }
  if (Grid::GridCmdOptionExists(argv, argv + argc, "--mobius_b")) {
    HMCParams.b = std::stod(Grid::GridCmdOptionPayload(argv, argv + argc, "--mobius_b"));
  }
  if (Grid::GridCmdOptionExists(argv, argv + argc, "--mobius_c")) {
    HMCParams.c = std::stod(Grid::GridCmdOptionPayload(argv, argv + argc, "--mobius_c"));
  }
  if (Grid::GridCmdOptionExists(argv, argv + argc, "--Ls")) {
    HMCParams.Ls = std::stoi(Grid::GridCmdOptionPayload(argv, argv + argc, "--Ls"));
  }
  if (Grid::GridCmdOptionExists(argv, argv + argc, "--cnfg_dir")) {
    HMCParams.cnfg_dir = Grid::GridCmdOptionPayload(argv, argv + argc, "--cnfg_dir");
  }

  return HMCParams;
}



int main(int argc, char **argv) {
  using namespace Grid;
  hmc_params HMCParams = ReadCommandLineHMC(argc, argv);

  Grid_init(&argc, &argv);
  int threads = GridThread::GetThreads();
  // here make a routine to print all the relevant information on the run
  std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;

   // Typedefs to simplify notation
  typedef SpWilsonImplR FermionImplPolicy;
  //typedef WilsonImplR FermionImplPolicy;
  typedef MobiusFermion<SpWilsonImplD> FermionAction;
  //typedef MobiusFermion<WilsonImplD> FermionAction;
  //typedef MobiusFermionD FermionAction;
  typedef typename FermionAction::FermionField FermionField;
  typedef Representations<SpFundamentalRepresentation> TheRepresentations;
  //typedef Representations<FundamentalRepresentation> TheRepresentations;

  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  //  typedef GenericHMCRunner<LeapFrog> HMCWrapper;
  //  MD.name    = std::string("Leap Frog");
  //  typedef GenericHMCRunner<ForceGradient> HMCWrapper;
  //  MD.name    = std::string("Force Gradient");
  typedef GenericHMCRunnerHirep<TheRepresentations, MinimumNorm2> HMCWrapper;

  HMCWrapper TheHMC;
  TheHMC.Parameters.MD.MDsteps = HMCParams.nsteps;
  TheHMC.Parameters.MD.trajL = HMCParams.tlen;

  TheHMC.Parameters.StartTrajectory = HMCParams.starttraj;

  // Grid from the command line arguments --grid and --mpi
  TheHMC.Resources.AddFourDimGrid("gauge"); // use default simd lanes decomposition

  CheckpointerParameters CPparams;
  CPparams.config_prefix = HMCParams.cnfg_dir + "/ckpoint_EODWF_lat";
  CPparams.rng_prefix    = HMCParams.cnfg_dir + "/ckpoint_EODWF_rng";
  CPparams.saveInterval  = HMCParams.save_freq;
  CPparams.format        = "IEEE64BIG";
  TheHMC.Resources.LoadNerscCheckpointer(CPparams);

  RNGModuleParameters RNGpar;
  RNGpar.serial_seeds = HMCParams.serial_seed;
  RNGpar.parallel_seeds = HMCParams.parallel_seed;
  TheHMC.Resources.SetRNGSeeds(RNGpar);

  // Construct observables
  // here there is too much indirection
  typedef PlaquetteMod<HMCWrapper::ImplPolicy> PlaqObs;
  TheHMC.Resources.AddObservable<PlaqObs>();
  //////////////////////////////////////////////

  const int Ls      = HMCParams.Ls;
  Real beta         = HMCParams.beta;
  Real light_mass = HMCParams.m;
  RealD M5  = HMCParams.M5;
  RealD b   = HMCParams.b;
  RealD c   = HMCParams.c;


  std::cout << "M5: " << M5 << std::endl;
  std::cout << "Mass: " << light_mass << std::endl;
  std::cout << "Ls: " << Ls << std::endl;
  std::cout << "b: " << b << std::endl;
  std::cout << "c: " << c << std::endl;
  std::cout << "starttraj: " << HMCParams.starttraj << std::endl;
  std::cout << "beta: " << HMCParams.beta << std::endl;
  auto GridPtr   = TheHMC.Resources.GetCartesian();
  auto GridRBPtr = TheHMC.Resources.GetRBCartesian();
  auto FGrid     = SpaceTimeGrid::makeFiveDimGrid(Ls,GridPtr);
  auto FrbGrid   = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,GridPtr);

  SpWilsonGaugeActionR GaugeAction(beta);
  //WilsonGaugeActionR GaugeAction(beta);

  // temporarily need a gauge field
  SpFundamentalRepresentation::LatticeField U(GridPtr);
  //FundamentalRepresentation::LatticeField U(GridPtr);

  // These lines are unecessary if BC are all periodic
  std::vector<Complex> boundary = {1,1,1,-1};
  FermionAction::ImplParams Params(boundary);

  double StoppingCondition = 1e-10;
  double MaxCGIterations = 30000;
  ConjugateGradient<FermionField>  CG(StoppingCondition,MaxCGIterations);

  ////////////////////////////////////
  // Collect actions
  ////////////////////////////////////
  ActionLevel<HMCWrapper::Field, TheRepresentations> Level1(1);
  ActionLevel<HMCWrapper::Field, TheRepresentations> Level2(4);

  FermionAction FermOp(U, *FGrid, *FrbGrid, *GridPtr, *GridRBPtr, light_mass, M5, b, c, Params);
  TwoFlavourEvenOddPseudoFermionAction<FermionImplPolicy> Nf2(FermOp, CG, CG);
  Nf2.is_smeared = false;

  Level1.push_back(&Nf2);

  /////////////////////////////////////////////////////////////
  // Gauge action
  /////////////////////////////////////////////////////////////
  Level2.push_back(&GaugeAction);
  TheHMC.TheAction.push_back(Level1);
  TheHMC.TheAction.push_back(Level2);
  std::cout << GridLogMessage << " Action complete "<< std::endl;

  /////////////////////////////////////////////////////////////
  // HMC parameters are serialisable

  std::cout << GridLogMessage << " Running the HMC "<< std::endl;
  TheHMC.ReadCommandLine(argc, argv); // these can be parameters from file
  TheHMC.Run();  // no smearing

  Grid_finalize();
} // main

