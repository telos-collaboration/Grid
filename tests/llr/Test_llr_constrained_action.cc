/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./HMC/LLR_HMC_SpWilsonGauge.cc

Copyright (C) 2015

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
Author: Frederic Bonnet  <frederic.bonnet@swansea.ac.uk>

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
// System
#include <iostream>
#include <fstream>
//#include <iomanip>
// Application includes
#include <Grid/Grid.h>

#include <Grid/qcd/llr_hmc/llr_hmc.h>

// Colors for test result printing
#define C_RED     ANSI_COLOR_BRIGHT_RED
#define C_GREEN   ANSI_COLOR_BRIGHT_GREEN
#define C_RESET   ANSI_COLOR_RESET //!<Short hand notation for color reset

/////////////////////////////////////////////////////////////
/// Main
/////////////////////////////////////////////////////////////
int main(int argc, char **argv) {
  
  std::cout << Grid::GridLogLLR <<"<---- Start Test_llr_contrained_action.cc ---->"<<std::endl;
  // Initializing Grid library environment.
  Grid::Grid_init(&argc, &argv);
  Grid::GridLogLayout();
  
  // Instantiating the inout parameter structure as a pointer.
  // IC TODO: do these structs have to be pointers?
  namespace_LLR::llrparams* s_llrparams_in =
    (struct namespace_LLR::llrparams*)malloc(sizeof(struct namespace_LLR::llrparams));
  // Initialising the structure
  s_llrparams_in->nrm = 700;
  s_llrparams_in->nth = 300;
  s_llrparams_in->it = 1;
  s_llrparams_in->umb_RM_freq = 1;
  s_llrparams_in->umb_meas_freq = 1;
  s_llrparams_in->umb_therm_freq = 1;
  s_llrparams_in->cfactor = 1;
  s_llrparams_in->starta = 5.66;
  s_llrparams_in->a = 5.66; // TODO: ?-- try the wrong value for interval --?
  s_llrparams_in->S0 = 13281.000;
  s_llrparams_in->dS = 3.0;
  // initialising the action and plaquette in the struc.
  s_llrparams_in->S = 1234.5678;
  s_llrparams_in->plaq = 4321.9876;

  // Create and initialise to some arbitrary values the hmc_parameters structure
  namespace_LLR::hmc_params_llr* s_hmc_params_llr_in =
    (struct namespace_LLR::hmc_params_llr*) malloc(sizeof (struct namespace_LLR::hmc_params_llr));
  // Initialising the hmc_params_llr structure
  s_hmc_params_llr_in->saveInterval = 1;
  s_hmc_params_llr_in->StartTrajectory = 1;
  s_hmc_params_llr_in->beta = 2.4;
  s_hmc_params_llr_in->trajL = 1.0;
  s_hmc_params_llr_in->MDsteps = 44;
  s_hmc_params_llr_in->Trajectories = 40;
  s_hmc_params_llr_in->Thermalizations = 20;

  // Create the llr_hmc object
  namespace_LLR::llr_hmc* p_llr_hmc_main_o = new namespace_LLR::llr_hmc();
  // Get the command line reader parameters and update data structure
  s_hmc_params_llr_in = p_llr_hmc_main_o->ReadCommandLine(argc, argv,s_hmc_params_llr_in);

  // Initialize robbins_monro
  p_llr_hmc_main_o->init_robbins_monro(s_llrparams_in);
  // printing the structures to see how they are constructed.
  p_llr_hmc_main_o->print_s_llrparams(s_llrparams_in);
  p_llr_hmc_main_o->print_s_hmc_params_llr(s_hmc_params_llr_in);


  // Distribution test
  int na = 1;
  Grid::RealD sigma = 3.0;
  Grid::RealD a_del = 0.0;
  Grid::RealD a_int = s_llrparams_in->starta - sigma;
  Grid::RealD a_fin = s_llrparams_in->starta + sigma;
  std::vector<Grid::RealD> a_array(na);
  if (na>1) {a_del = (a_fin - a_int) / (na - 1);}
  
  for (int i = 0; i < na; i++) {a_array[i] = a_int + i * a_del;}
  
  int mid = na / 2;
  a_array[mid] = s_llrparams_in->starta;
  
  // Printing out the vector list
  for (int i = 0; i < na; i++) {
    std::cout << Grid::GridLogLLR
	      << "Test_llr_contrained value     a ----->: ["
	      << std::setw(4) << i << "]: --->: " << a_array[i] << std::endl;
  }
  
  // The epsilon for the assertion
  Grid::RealD epsilon_a = 0;
  Grid::RealD epsilon_plaquette = 0.001;
  Grid::RealD epsilon_S = 20.0;
  if (na == 2) {
    epsilon_a = std::abs(a_fin - a_int) / 2;
    epsilon_plaquette = 0.1;  // For multiple a values use 0.1
    epsilon_S = 1000.0;       // For multiple a values use 1000.0
  } else if (na > 2){
    epsilon_a = std::abs(a_fin - a_int);
    epsilon_plaquette = 0.1;  // For multiple a values use 0.1
    epsilon_S = 1000.0;       // For multiple a values use 1000.0
  }
  std::cout << Grid::GridLogLLR
	    << "Epsilon S --->: MDsteps[" << s_hmc_params_llr_in->MDsteps << "] --->: "
	    << epsilon_S << std::endl;
  std::cout << Grid::GridLogLLR
	    << "Epsilon P --->: MDsteps[" << s_hmc_params_llr_in->MDsteps << "] --->: "
	    << epsilon_plaquette << std::endl;
  std::cout << Grid::GridLogLLR
	    << "Epsilon a --->: MDsteps[" << s_hmc_params_llr_in->MDsteps << "] --->: "
	    << epsilon_a << std::endl;

  // The expected values
  Grid::RealD expected_plaquette = 0.459;
  Grid::RealD expected_action = s_llrparams_in->S0;

  // Main loop for multiple a tests The runs will take about 30 minutes on 1 node
  for (int i = 0; i < na; i++) {
    
    std::cout << Grid::GridLogLLR << "Loop over a starts              ----->: [...]" <<std::endl;
    
    // Mapping the screening value from a_array[i] into the struct
    s_llrparams_in->a = a_array[i];    
    std::cout << Grid::GridLogLLR << "Test_llr_constrained value     a ----->: ["
	      << std::setw(4) << i << "]: --->: " << s_llrparams_in->a << std::endl;

    // Constructing the output files using the command line input parameters
    // Creating the output file
    std::cout << Grid::GridLogLLR << "<---- Creating/Opening output log file ... --->" << std::endl;
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(4) << s_llrparams_in->a;

    std::string filename_log = "test_llr_constrained_action_MDsteps-" + std::to_string(s_hmc_params_llr_in->MDsteps) +"-a" + oss.str() + ".log";
    std::string filename_csv = "test_llr_constrained_action_MDsteps-" + std::to_string(s_hmc_params_llr_in->MDsteps) +"-a" + oss.str() + ".csv";
    std::ofstream run_LLR_HMC_logfile(filename_log);
    std::ofstream run_LLR_HMC_csvfile(filename_csv);
    
    // Writing the header to csv file
    run_LLR_HMC_csvfile
      << "Grid::GridLogMessage"
      << "," << "traj"
      << "," << "beta"
      << "," << "MDsteps"
      << "," << "plaq"
      << "," << "action"
      << "," << "a"
      << "," << "S0"
      << "," << "vol"
      << std::endl;
    
    run_LLR_HMC_logfile << "<---- Start Test_llr_constrained_action.cc ---->" << std::endl;

    // Start of the main commands.
    bool with_llr = false;
    std::cout << Grid::GridLogLLR << "<---- with_llr (initialized)     ---->: " << with_llr << std::endl;
    std::cout << Grid::GridLogLLR << "<---- llr_config                 ---->: " << llr_config << std::endl;
#if defined(llr_config)
    with_llr = true;
#endif
    std::cout << Grid::GridLogLLR << "<---- with_llr (llr_config)      ---->: " << with_llr << std::endl;
    std::cout << Grid::GridLogLLR << "<---- Sp2n_config                ---->: " << Sp2n_config << std::endl;
    
    // if --enable--LLR switch is activated
    // IC TODO: this should be controlled by the cmake, not here. This test should be
    //          built only if llr is enabled.
    if (with_llr) {
      std::cout << Grid::GridLogLLR << "Start of if block ...  with_llr  ---->" <<std::endl;
      // Sp(2n) representation
      typedef Grid::GenericHMCRunnerSpLLR<Grid::MinimumNorm2> HMCWrapperSpLLR;
      HMCWrapperSpLLR TheHMC;
      
      // Grid from the command line
      TheHMC.ReadCommandLine(argc, argv);
      TheHMC.Resources.AddFourDimGrid("gauge");
      
      // Getting the check point parameters setup
      Grid::CheckpointerParameters CPparams;
      CPparams.config_prefix = "ckpoint_lat";
      CPparams.rng_prefix = "ckpoint_rng";
      CPparams.saveInterval = s_hmc_params_llr_in->saveInterval; //100;
      CPparams.format = "IEEE64BIG";
      TheHMC.Resources.LoadNerscCheckpointer(CPparams);
      // Getting the random numbers setup
      Grid::RNGModuleParameters RNGpar;
      RNGpar.serial_seeds = "1 2 3 4 5";
      RNGpar.parallel_seeds = "6 7 8 9 10";
      TheHMC.Resources.SetRNGSeeds(RNGpar);
      
      // Construct observables
      typedef Grid::PlaquetteMod<HMCWrapperSpLLR::ImplPolicy> PlaqObs;
      TheHMC.Resources.AddObservable<PlaqObs>();
      
      typedef Grid::TopologicalChargeMod<HMCWrapperSpLLR::ImplPolicy> QObs;
      Grid::TopologyObsParameters TopParams;
      TopParams.interval = 1;
      TopParams.do_smearing = false;
      TopParams.Smearing.init_step_size = 0.01;
      TopParams.Smearing.tolerance = 1e-5;
      //TopParams.Smearing.steps = 200;
      //TopParams.Smearing.step_size = 0.01;
      TopParams.Smearing.meas_interval = 50;
      TopParams.Smearing.maxTau = 2.0;
      TheHMC.Resources.AddObservable<QObs>(TopParams);

      typedef LLRPlaquetteMod<HMCWrapperSpLLR::ImplPolicy> LLRSpPlaqObs;
      TheHMC.Resources.AddObservable<LLRSpPlaqObs>();
      
      typedef LLRActionMod<HMCWrapperSpLLR::ImplPolicy> LLRSpActObs;
      ActionLoggerObsParameters ActParams;
      ActParams.a_ = s_llrparams_in->a;
      ActParams.beta_ = s_hmc_params_llr_in->beta;
      ActParams.MDsteps_ = s_hmc_params_llr_in->MDsteps;
      ActParams.s_llrparams_in_ = s_llrparams_in;
      ActParams.logFile_ = &run_LLR_HMC_logfile;
      ActParams.csvFile_ = &run_LLR_HMC_csvfile;
      TheHMC.Resources.AddObservable<LLRSpActObs>(ActParams);
      
      // Collect actions, here use more encapsulation
      // standard LLRGaugeAction with Sp(2n) gauge action
      typedef Grid::LLRGaugeAction<Grid::SpWilsonGaugeActionR, Grid::PeriodicGimplR> LLRGaugeActionR;
      LLRGaugeActionR LLRaction(s_llrparams_in, s_hmc_params_llr_in->beta);
      
      Grid::ActionLevel<HMCWrapperSpLLR::Field> Level1(1);
      Level1.push_back(&LLRaction);
      TheHMC.TheAction.push_back(Level1);
      
      // HMC parameters
      TheHMC.Parameters.MD.MDsteps = s_hmc_params_llr_in->MDsteps;      // 40;
      TheHMC.Parameters.MD.trajL   = float(s_hmc_params_llr_in->trajL); // 1.0;
      
      TheHMC.Run();  // no smearing
      
      // Creating the test PASS/FAIL
      std::cout << Grid::GridLogLLR << "--------------------------------------------------"<<std::endl;
      std::cout << Grid::GridLogLLR << "Final action and Plaquette:" << std::endl;
      
      // Compute diffs
      Grid::RealD diff_plaquette = std::abs(s_llrparams_in->plaq - expected_plaquette);
      Grid::RealD diff_action    = std::abs(s_llrparams_in->S    - expected_action);
      Grid::RealD diff_a         = std::abs(s_llrparams_in->a    - s_llrparams_in->starta);
      
      // Determine pass/fail
      bool pass_plaquette = (diff_plaquette < epsilon_plaquette);
      bool pass_action    = (diff_action    < epsilon_S);
      bool pass_a         = (diff_a         <= epsilon_a);

      // Print results
      std::cout << Grid::GridLogLLR
		<< "starta[" << s_llrparams_in->starta << "] --->: "
		<< "current a[" << s_llrparams_in->a << "]" << std::endl;  
      std::cout << Grid::GridLogLLR
		<< "a values  --->: MDsteps[" << s_hmc_params_llr_in->MDsteps << "] --->: "
		<< std::setw(12) << s_llrparams_in->a
		<< "   (starta    " << s_llrparams_in->starta
		<< ", diff = " << std::setw(14) << diff_a
		<< ")  ==> " << (pass_a ? C_GREEN "PASS" : C_RED "FAIL") << C_RESET
                << std::endl;
      std::cout << Grid::GridLogLLR
		<< "Action    --->: MDsteps[" << s_hmc_params_llr_in->MDsteps << "] --->: "
		<< std::setw(12) << s_llrparams_in->S
		<< "   (expected " << expected_action
		<< ", diff = " << std::setw(14) << diff_action
		<< ")  ==> " << (pass_action ? C_GREEN "PASS" : C_RED "FAIL") << C_RESET
                << std::endl;
      std::cout << Grid::GridLogLLR
		<< "Plaquette --->: MDsteps[" << s_hmc_params_llr_in->MDsteps << "] --->: "
		<< std::setw(12) <<s_llrparams_in->plaq
		<< "   (expected " << expected_plaquette
		<< ", diff = " << std::setw(14) << diff_plaquette
		<< ")  ==> " << (pass_plaquette ? C_GREEN "PASS" : C_RED "FAIL") << C_RESET
		<< std::endl;

      SOFT_ASSERT(diff_plaquette < epsilon_plaquette,
		  "Plaquette within tolerance: value = " +
		  std::to_string(s_llrparams_in->plaq) +
		  "     expected = " + std::to_string(expected_plaquette) +
		  "     diff = " + std::to_string(diff_plaquette));
      
      SOFT_ASSERT(diff_action < epsilon_S,
		  "Action within tolerance   : value = " +
		  std::to_string(s_llrparams_in->S) +
		  " expected = " + std::to_string(expected_action) +
		  " diff = " + std::to_string(diff_action));
      
      // Assert plaquette
      assert( diff_plaquette < epsilon_plaquette &&
	      "Plaquette assertion FAILED: |plaq - expected_plaq| >= epsilon" );
      
      // Assert action
      assert( diff_action < epsilon_S &&
	      "Action assertion FAILED: |S - expected_S| >= epsilon" );
      
      assert((diff_plaquette < epsilon_plaquette) &&
	     ("Plaquette FAIL: value=" +
	      std::to_string(s_llrparams_in->plaq) +
	      " expected=" + std::to_string(expected_plaquette) +
	      " diff=" + std::to_string(diff_plaquette)).c_str());
      
      assert((diff_action < epsilon_S) &&
	     ("Action FAIL: value=" +
	      std::to_string(s_llrparams_in->S) +
	      " expected=" + std::to_string(expected_action) +
	      " diff=" + std::to_string(diff_action)).c_str());
      
      std::cout << Grid::GridLogLLR << "End of if block ...    with_llr ----->: "<< with_llr << std::endl;
    } /* [end-if] with_llr */
    
    
    // End statement
    std::cout << Grid::GridLogLLR << "<---- End Test_llr_contrained_action.cc ---->" << std::endl;
    run_LLR_HMC_logfile << "<---- End Test_llr_contrained_action.cc ---->" << std::endl;
    // Closing log file
    run_LLR_HMC_logfile.close();
    run_LLR_HMC_csvfile.close();

} /* [end-For-Loop] Main loop (int i = 0; i < na; i++) */

  // Finalising Grid environment.
  Grid::Grid_finalize();
}
