/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./HMC/HMC_WilsonGauge.cc

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
#include <Grid/Grid.h>

// Application includes

int main(int argc, char **argv) {
    std::cout<<"Start HMC_SUn_WilsonGauge.cc" <<std::endl;
    // Initializing Grid library

    Grid::Grid_init(&argc, &argv);
    Grid::GridLogLayout();

    // SU(N) representation
    typedef Grid::GenericHMCRunner <Grid::MinimumNorm2> HMCWrapper;
    HMCWrapper TheHMC;

    // Grid from the command line
    TheHMC.ReadCommandLine(argc, argv);
    TheHMC.Resources.AddFourDimGrid("gauge");
    // Possible to create the module by hand
    // hardcoding parameters or using a Reader

    // Checking the parameters in the used, we will use the same as the
    // Standard wilson ones.
    Grid::CheckpointerParameters CPparams;
    CPparams.config_prefix = "ckpoint_lat";
    CPparams.rng_prefix = "ckpoint_rng";
    CPparams.saveInterval = 100;
    CPparams.format = "IEEE64BIG";

    TheHMC.Resources.LoadNerscCheckpointer(CPparams);

    Grid::RNGModuleParameters RNGpar;
    RNGpar.serial_seeds = "1 2 3 4 5";
    RNGpar.parallel_seeds = "6 7 8 9 10";
    TheHMC.Resources.SetRNGSeeds(RNGpar);

    // Construct observables
    // here there is too much indirection
    typedef Grid::PlaquetteMod<HMCWrapper::ImplPolicy > PlaqObs;
    typedef Grid::TopologicalChargeMod<HMCWrapper::ImplPolicy> QObs;

    TheHMC.Resources.AddObservable<PlaqObs>();
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
    //////////////////////////////////////////////
    /////////////////////////////////////////////////////////////
    // Collect actions, here use more encapsulation
    // need wrappers of the fermionic classes
    // that have a complex construction
    // standard
    Grid::RealD beta = 6.0;

    typedef Grid::WilsonGaugeAction<Grid::PeriodicGimplR> GaugeActionR;

    GaugeActionR action(beta);

    Grid::ActionLevel<HMCWrapper::Field> Level1(1);
    Level1.push_back(&action);
    //Level1.push_back(WGMod.getPtr());
    TheHMC.TheAction.push_back(Level1);
    /////////////////////////////////////////////////////////////

    // HMC parameters are serialisable
    TheHMC.Parameters.MD.MDsteps = 20;
    TheHMC.Parameters.MD.trajL   = 1.0;

    TheHMC.Run();  // no smearing


    // Finalizing the Grid library
    Grid::Grid_finalize();
    // End statement
    std::cout<<"<---- End HMC_SUn_WilsonGauge.cc ---->" <<std::endl;

}