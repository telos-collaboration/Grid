/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/hmc/HMC.h

Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: neo <cossu@post.kek.jp>
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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
			   /*  END LEGAL */
			   //--------------------------------------------------------------------
			   /*! @file HMC.h
			    * @brief Classes for Hybrid Monte Carlo update
			    *
			    * @author Guido Cossu
			    * @author Peter Boyle
			    */
			   //--------------------------------------------------------------------
#pragma once

#include <string>
#include <list>

#include <Grid/qcd/hmc/integrators/Integrator.h>
#include <Grid/qcd/hmc/integrators/Integrator_algorithm.h>

// Temporary include for now may be removed later
#include <Grid/qcd/llr_hmc/llr_hmc.h>

NAMESPACE_BEGIN(Grid);

struct HMCparameters: Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(HMCparameters,
                                  Integer, StartTrajectory,
                                  Integer, Trajectories, /* @brief Number of sweeps in this run */
                                  bool, MetropolisTest,
                                  Integer, NoMetropolisUntil,
				  bool, PerformRandomShift, /* @brief Randomly shift the gauge configuration at the start of a trajectory */
                                  std::string, StartingType,
                                  IntegratorParameters, MD)

  HMCparameters() {
    ////////////////////////////// Default values
    MetropolisTest    = true;
    NoMetropolisUntil = 10;
    StartTrajectory   = 0;
    Trajectories      = 10;
    StartingType      = "HotStart";
    PerformRandomShift = true;
    /////////////////////////////////
  }

  template <class ReaderClass >
  HMCparameters(Reader<ReaderClass> & TheReader){
    initialize(TheReader);
  }

  template < class ReaderClass > 
  void initialize(Reader<ReaderClass> &TheReader){
    std::cout << GridLogMessage << "Reading HMC\n";
    read(TheReader, "HMC", *this);
  }


  void print_parameters() const {
    std::cout << GridLogMessage << "[HMC parameters] Trajectories            : " << Trajectories << "\n";
    std::cout << GridLogMessage << "[HMC parameters] Start trajectory        : " << StartTrajectory << "\n";
    std::cout << GridLogMessage << "[HMC parameters] Metropolis test (on/off): " << std::boolalpha << MetropolisTest << "\n";
    std::cout << GridLogMessage << "[HMC parameters] Thermalization trajs    : " << NoMetropolisUntil << "\n";
    std::cout << GridLogMessage << "[HMC parameters] Doing random shift      : " << std::boolalpha << PerformRandomShift << "\n";
    std::cout << GridLogMessage << "[HMC parameters] Starting type           : " << StartingType << "\n";
    MD.print_parameters();
  }
  
};
	
template <class IntegratorType>
class HybridMonteCarlo {
private:
  const HMCparameters Params;

  // llr switch passed in the constructor of
  bool with_llr = false;

  typedef typename IntegratorType::Field Field;
  typedef typename IntegratorType::FieldImplementation FieldImplementation;
  typedef std::vector< HmcObservable<Field> * > ObsListType;

  //pass these from the resource manager
  GridSerialRNG &sRNG;   
  GridParallelRNG &pRNG; 

  Field &Ucur;
  
  IntegratorType &TheIntegrator;
  ObsListType Observables;

  /////////////////////////////////////////////////////////
  // Metropolis step
  /////////////////////////////////////////////////////////
  bool metropolis_test(const RealD DeltaH) {
    RealD rn_test;

    RealD prob = std::exp(-DeltaH);

    random(sRNG, rn_test);

    std::cout << GridLogHMC << "--------------------------------------------------\n";
    std::cout << GridLogHMC << "exp(-dH) = " << prob << "  Random = " << rn_test << "\n";
    std::cout << GridLogHMC << "Acc. Probability = " << ((prob < 1.0) ? prob : 1.0) << "\n";

    if ((prob > 1.0) || (rn_test <= prob)) {  // accepted
      std::cout << GridLogHMC << "Metropolis_test -- ACCEPTED\n";
      std::cout << GridLogHMC << "--------------------------------------------------\n";
      return true;
    } else {  // rejected
      std::cout << GridLogHMC << "Metropolis_test -- REJECTED\n";
      std::cout << GridLogHMC << "--------------------------------------------------\n";
      return false;
    }
  }

  /////////////////////////////////////////////////////////
  // Evolution
  /////////////////////////////////////////////////////////
  RealD evolve_hmc_step(Field &U) {

    GridBase *Grid = U.Grid();

    if(Params.PerformRandomShift){
#if 0
      //////////////////////////////////////////////////////////////////////////////////////////////////////
      // Mainly for DDHMC perform a random translation of U modulo volume
      //////////////////////////////////////////////////////////////////////////////////////////////////////
      std::cout << GridLogMessage << "--------------------------------------------------\n";
      std::cout << GridLogMessage << "Random shifting gauge field by [";

      std::vector<typename FieldImplementation::GaugeLinkField> Umu(Grid->Nd(), U.Grid());
      for(int mu=0;mu<Grid->Nd();mu++) Umu[mu] = PeekIndex<LorentzIndex>(U, mu);

      for(int d=0;d<Grid->Nd();d++) {

	int L = Grid->GlobalDimensions()[d];

	RealD rn_uniform;  random(sRNG, rn_uniform);

	int shift = (int) (rn_uniform*L);

	std::cout << shift;
	if(d<Grid->Nd()-1) std::cout <<",";
	else               std::cout <<"]\n";
      
	//shift all fields together in a way that respects the gauge BCs
	for(int mu=0; mu < Grid->Nd(); mu++)
	  Umu[mu] = FieldImplementation::CshiftLink(Umu[mu],d,shift);

	for(int mu=0;mu<Grid->Nd();mu++) PokeIndex<LorentzIndex>(U,Umu[mu],mu);
      }
      std::cout << GridLogMessage << "--------------------------------------------------\n";
#endif	
    }

    TheIntegrator.reset_timer();

      std::cout << GridLogHMC << C_RED << "with_llr at Integrator.refresh ---->: "<< with_llr <<C_RESET <<std::endl;
      int set_llr_if = 0;
      if (with_llr == true) { set_llr_if = 1234; }
      std::cout << GridLogHMC << B_RED << "set_llr_if at Int.refresh      ---->: "<< set_llr_if <<C_RESET <<std::endl;

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // set U and initialize P and phi's
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    std::cout << GridLogMessage << "--------------------------------------------------\n";
    std::cout << GridLogMessage << "Refresh momenta and pseudofermions"<< std::endl;
    TheIntegrator.refresh(U, sRNG, pRNG);
    std::cout << GridLogMessage << "--------------------------------------------------\n";
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // initial state action
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    std::cout << GridLogMessage << "--------------------------------------------------\n";
    std::cout << GridLogMessage << "Compute initial action"<< std::endl;

      RealD S0 = 0.0;
      RealD H0 = 0.0;
      RealD S0_llr = 0.0;
      //RealD H0_llr = 0.0;
      if (set_llr_if == 1234) {
          //////////////////////////////////////////////////////////////////////////////////////////////////////
          // Initial state action with LLR
          //////////////////////////////////////////////////////////////////////////////////////////////////////
          S0_llr = TheIntegrator.Sinitial_llr(U);

          H0 = S0_llr;  // Mapping back to the original variable name to minimize the footprint
          std::cout << GridLogHMC << C_MAGENTA << "H0_llr from Sinitial_llr(U)     ----->: " << H0 << C_RESET << std::endl;

      } else if (set_llr_if == 0) {
          //////////////////////////////////////////////////////////////////////////////////////////////////////
          // Initial state action with no LLR
          //////////////////////////////////////////////////////////////////////////////////////////////////////
          H0 = TheIntegrator.Sinitial(U);
          std::cout << GridLogHMC << C_MAGENTA << "H0 from Sinitial(U)             ----->: " << H0 << C_RESET<< std::endl;
      }
      //////////////////////////////////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////////////////////////////////

    std::cout << GridLogMessage << "--------------------------------------------------\n";

    std::streamsize current_precision = std::cout.precision();
    std::cout.precision(15);
    std::cout << GridLogHMC << "Total H before trajectory = " << H0 << "\n";
    std::cout.precision(current_precision);

    std::cout << GridLogMessage << "--------------------------------------------------\n";
    std::cout << GridLogMessage << " Molecular Dynamics evolution "<< std::endl;
    TheIntegrator.integrate(U);
    std::cout << GridLogMessage << "--------------------------------------------------\n";

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // updated state action
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    std::cout << GridLogMessage << "--------------------------------------------------\n";
    std::cout << GridLogMessage << "Compute final action"<< std::endl;

      RealD H1 = 0.0;
      RealD S1_llr = 0.0;
      //RealD H1_llr = 0.0;
      if (set_llr_if == 1234) {
          ////////////////////////////////////////////////////////////////////////////////////////////////
          // final state action with LLR
          ////////////////////////////////////////////////////////////////////////////////////////////////
          S1_llr = TheIntegrator.S_llr(U);
          H1 = S1_llr;
          std::cout << GridLogHMC << C_MAGENTA << "H1_llr from S_llr(U)            ----->: " << H1 << C_RESET << std::endl;
          std::cout << GridLogMessage << "--------------------------------------------------\n";
      } else if (set_llr_if == 0) {
          ////////////////////////////////////////////////////////////////////////////////////////////////
          // final state action with no LLR
          ////////////////////////////////////////////////////////////////////////////////////////////////
          H1 = TheIntegrator.S(U);
          std::cout << GridLogHMC << C_MAGENTA << "H1 from S(U)                    ----->: " << H1 << C_RESET << std::endl;
          std::cout << GridLogMessage << "--------------------------------------------------\n";
      }
    ///////////////////////////////////////////////////////////
    if(0){
      std::cout << "------------------------- Reversibility test" << std::endl;

      TheIntegrator.reverse_momenta();
      TheIntegrator.integrate(U);

        if (set_llr_if == 1234) {
            H1 = TheIntegrator.S_llr(U);
            std::cout << GridLogHMC << B_YELLOW  << "H1(llr) Reversibility test      ----->: "<< H1 << C_RESET <<std::endl;
        } else if (set_llr_if == 0) {
            H1 = TheIntegrator.S(U);  // updated state action
            std::cout << GridLogHMC << B_YELLOW  << "H1 Reversibility test           ----->: "<< H1 << C_RESET <<std::endl;
        }
      std::cout << "--------------------------------------------" << std::endl;
    }
    ///////////////////////////////////////////////////////////
    std::cout.precision(15);
      ////////////////////////////////////////////////////////////////////////////////////////////////////
      // Delta in the final and initial state action with LLR
      ////////////////////////////////////////////////////////////////////////////////////////////////////
      RealD dH = 0.0;
      //RealD dH_llr = 0.0;
      if (set_llr_if == 1234) {
          //dH_llr = H1_llr - H0_llr;
          dH = H1 - H0;
          std::cout << GridLogHMC << B_MAGENTA << "H1_llr HMC                      ----->: " << H1 << C_RESET << std::endl;
          std::cout << GridLogHMC << B_MAGENTA << "dH = H1 - H0                    ----->: " << dH << C_RESET << std::endl;
          std::cout << GridLogHMC << "--------------------------------------------------\n";
          std::cout << GridLogHMC << "Total H1_llr after trajectory  = " << H1 << "  the_deltaH = " << dH << "\n";
          std::cout << GridLogHMC << "--------------------------------------------------\n";
      } else if (set_llr_if == 0) {
          dH = H1 - H0;
          std::cout << GridLogHMC << B_YELLOW  << "H1 HMC                          ----->: "<< H1 << C_RESET <<std::endl;
          std::cout << GridLogHMC << B_YELLOW  << "dH = H1 - H0                    ----->: "<< dH << C_RESET <<std::endl;
          std::cout << GridLogHMC << "--------------------------------------------------\n";
          std::cout << GridLogHMC << "Total H after trajectory  = " << H1 << "  dH = " << dH << "\n";
          std::cout << GridLogHMC << "--------------------------------------------------\n";
      }
    std::cout.precision(current_precision);

    return (H1 - H0);
  }

public:
  /////////////////////////////////////////
  // public variables
  /////////////////////////////////////////

  /////////////////////////////////////////
  // Constructor
  /////////////////////////////////////////
  HybridMonteCarlo(HMCparameters _Pams, IntegratorType &_Int,
                   GridSerialRNG &_sRNG, GridParallelRNG &_pRNG, 
                   ObsListType _Obs, Field &_U,
                   bool _with_llr)
    : Params(_Pams), TheIntegrator(_Int), sRNG(_sRNG), pRNG(_pRNG), Observables(_Obs), Ucur(_U), with_llr(_with_llr) {}
  ~HybridMonteCarlo(){};

  void evolve(void) {
    Real DeltaH;

      std::cout << GridLogHMC << "with_llr in       ---------->: "<< with_llr     <<std::endl;
      std::cout << GridLogHMC << "HybridMonteCarlo  ---------->: "<< __FUNCTION__ <<std::endl;

    Field Ucopy(Ucur.Grid());

    Params.print_parameters();
    TheIntegrator.print_actions();

    // Actual updates (evolve a copy Ucopy then copy back eventually)
    unsigned int FinalTrajectory = Params.Trajectories + Params.NoMetropolisUntil + Params.StartTrajectory;

    for (int traj = Params.StartTrajectory; traj < FinalTrajectory; ++traj) {

      std::cout << GridLogHMC << "-- # Trajectory = " << traj << "\n";

      if (traj < Params.StartTrajectory + Params.NoMetropolisUntil) {
      	std::cout << GridLogHMC << "-- Thermalization" << std::endl;
      }
      
      double t0=usecond();
      Ucopy = Ucur;

        std::cout << GridLogHMC << "\x1b[33m"<<"DeltaH b4 evolve_hmc_step      --------->: "<< DeltaH <<"\x1b[0m"<<std::endl;
      DeltaH = evolve_hmc_step(Ucopy);
        std::cout << GridLogHMC << "\x1b[33m"<<"DeltaH                         --------->: "<< DeltaH <<"\x1b[0m"<<std::endl;
      // Metropolis-Hastings test
      bool accept = true;
      if (Params.MetropolisTest && traj >= Params.StartTrajectory + Params.NoMetropolisUntil) {
        accept = metropolis_test(DeltaH);
          std::cout <<GridLogMessage << "DeltaH      (accept) ---------->: "<< DeltaH <<std::endl;
      } else {
      	std::cout << GridLogHMC << "Skipping Metropolis test" << std::endl;
      }

      if (accept)
        Ucur = Ucopy; 
      
      double t1=usecond();
      std::cout << GridLogHMC << "Total time for trajectory (s): " << (t1-t0)/1e6 << std::endl;

      TheIntegrator.print_timer();
      
      TheIntegrator.Smearer.set_Field(Ucur);
      for (int obs = 0; obs < Observables.size(); obs++) {
      	std::cout << GridLogDebug << "Observables # " << obs << std::endl;
      	std::cout << GridLogDebug << "Observables total " << Observables.size() << std::endl;
      	std::cout << GridLogDebug << "Observables pointer " << Observables[obs] << std::endl;
        Observables[obs]->TrajectoryComplete(traj + 1, TheIntegrator.Smearer, sRNG, pRNG);
      }
      std::cout << GridLogHMC << ":::::::::::::::::::::::::::::::::::::::::::" << std::endl;
    }
  }

};

NAMESPACE_END(Grid);
// april 11 2017 merge, Guido, commenting out
//#include <Grid/parallelIO/NerscIO.h>
//#include <Grid/qcd/hmc/NerscCheckpointer.h>
//#include <Grid/qcd/hmc/HmcRunner.h>

