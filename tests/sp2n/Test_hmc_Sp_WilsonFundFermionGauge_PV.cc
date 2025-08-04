/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./tests/Test_hmc_WilsonRatio.cc

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>

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

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
#include <Grid/Grid.h>

struct hmc_params {
  std::string save_path;
  double beta;
  double mass;
  double massPV;
};



hmc_params ReadCommandLineHMC(int argc, char** argv) {
  hmc_params HMCParams;
  if (Grid::GridCmdOptionExists(argv, argv + argc, "--savepath")) {
     HMCParams.save_path = Grid::GridCmdOptionPayload(argv, argv + argc, "--savepath");
  } else {
    std::cout << Grid::GridLogError << "--savepath must be specified" << std::endl;
  }

  if (Grid::GridCmdOptionExists(argv, argv + argc, "--beta")) {
     HMCParams.beta = std::stod(Grid::GridCmdOptionPayload(argv, argv + argc, "--beta"));
  } else {
    std::cout << Grid::GridLogError << "--beta must be specified" << std::endl;
  }

  if (Grid::GridCmdOptionExists(argv, argv + argc, "--mass")) {
     HMCParams.mass = std::stod(Grid::GridCmdOptionPayload(argv, argv + argc, "--mass"));
  } else {
    std::cout << Grid::GridLogError << "--mass must be specified" << std::endl;
  }

  if (Grid::GridCmdOptionExists(argv, argv + argc, "--massPV")) {
     HMCParams.massPV = std::stod(Grid::GridCmdOptionPayload(argv, argv + argc, "--massPV"));
  } else {
    std::cout << Grid::GridLogError << "--massPV must be specified" << std::endl;
  }


  return HMCParams;
}



int main(int argc, char **argv) {
  using namespace Grid;
   ;

  hmc_params HMCParams = ReadCommandLineHMC(argc, argv);
   Grid_init(&argc, &argv);
  int threads = GridThread::GetThreads();
  // here make a routine to print all the relevant information on the run
  std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;

  typedef Representations< SpFundamentalRepresentation > TheRepresentations;

   // Typedefs to simplify notation
  typedef GenericSpHMCRunnerHirep<TheRepresentations, MinimumNorm2> HMCWrapper;  // Uses the default minimum norm
  typedef SpWilsonImplR FermionImplPolicy;                                    // ok
  typedef SpWilsonFermionD FermionAction;
  typedef typename FermionAction::FermionField FermionField;


  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  HMCWrapper TheHMC;


  // Grid from the command line
  TheHMC.Resources.AddFourDimGrid("gauge");
  // Possibile to create the module by hand
  // hardcoding parameters or using a Reader


  CheckpointerParameters CPparams;
  CPparams.config_prefix = "./" + HMCParams.save_path +"ckpoint_lat";
  CPparams.rng_prefix = "./" + HMCParams.save_path + "ckpoint_rng";
  CPparams.saveInterval = 10000;
  CPparams.format = "IEEE64BIG";


  TheHMC.Resources.LoadNerscCheckpointer(CPparams);

  RNGModuleParameters RNGpar;
  RNGpar.serial_seeds = "1 2 3 4 5";
  RNGpar.parallel_seeds = "6 7 8 9 10";
  TheHMC.Resources.SetRNGSeeds(RNGpar);

  // Construct observables
  typedef PlaquetteMod<HMCWrapper::ImplPolicy> PlaqObs;
  TheHMC.Resources.AddObservable<PlaqObs>();
  //////////////////////////////////////////////

  /////////////////////////////////////////////////////////////
  // Collect actions, here use more encapsulation
  // need wrappers of the fermionic classes
  // that have a complex construction
  // standard
  RealD beta = HMCParams.beta ;
  SpWilsonGaugeActionR Waction(beta);
  auto GridPtr = TheHMC.Resources.GetCartesian();
  auto GridRBPtr = TheHMC.Resources.GetRBCartesian();

  // temporarily need a gauge field
  SpFundamentalRepresentation::LatticeField U(GridPtr);

  RealD mass = HMCParams.mass;
  RealD pv   = HMCParams.massPV;

  FermionAction DenOp(U, *GridPtr, *GridRBPtr, mass);
  FermionAction NumOp(U, *GridPtr, *GridRBPtr, pv);

  ConjugateGradient<FermionField> CG(1.0e-8, 2000);

  TwoFlavourRatioPseudoFermionAction<FermionImplPolicy> Nf2(NumOp, DenOp,CG,CG);

  // Set smearing (true/false), default: false
  Nf2.is_smeared = false;


    // Collect actions
  ActionLevel<HMCWrapper::Field, TheRepresentations> Level1(1);
  Level1.push_back(&Nf2);

  ActionLevel<HMCWrapper::Field, TheRepresentations> Level2(4);
  Level2.push_back(&Waction);

  TheHMC.TheAction.push_back(Level1);
  TheHMC.TheAction.push_back(Level2);
  /////////////////////////////////////////////////////////////
  /*
  double rho = 0.1;  // smearing parameter
  int Nsmear = 3;    // number of smearing levels
  Smear_Stout<HMCWrapper::ImplPolicy> Stout(rho);
  SmearedConfiguration<HMCWrapper::ImplPolicy> SmearingPolicy(GridPtr, Nsmear, Stout);
  */
  // HMC parameters are serialisable
  TheHMC.Parameters.MD.MDsteps = 6;
  TheHMC.Parameters.MD.trajL   = 3.0;

  TheHMC.ReadCommandLine(argc, argv); // these can be parameters from file
  TheHMC.Run();  // no smearing
  //TheHMC.Run(SmearingPolicy); // for smearing

  Grid_finalize();

} // main
