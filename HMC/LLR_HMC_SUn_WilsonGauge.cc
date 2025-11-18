/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./HMC/LLR_HMC_WilsonGauge.cc

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

#include <Grid/qcd/llr_hmc/llr_hmc.h>

int main(int argc, char **argv) {
    std::cout<<"Start LLR_HMC_SUn_WilsonGauge.cc" <<std::endl;
    // Initializing Grid library
    Grid::Grid_init(&argc, &argv);
    Grid::GridLogLayout();

    // Instantiating the inout parameter structure as a pointer.
    namespace_LLR::llrparams* s_llrparams_in = (struct namespace_LLR::llrparams*)malloc(sizeof(struct namespace_LLR::llrparams));
    // Initialising the structure
    s_llrparams_in->nrm = 700;
    s_llrparams_in->nth = 300;
    s_llrparams_in->it = 1;
    s_llrparams_in->umb_RM_freq = 1;
    s_llrparams_in->umb_meas_freq = 1;
    s_llrparams_in->umb_therm_freq = 1;
    s_llrparams_in->cfactor = 1;
    s_llrparams_in->starta = 5.6;
    s_llrparams_in->a = 5.66;
    s_llrparams_in->S0 = 13281.000;
    s_llrparams_in->dS = 3.0;

    namespace_LLR::llr_hmc* p_llr_hmc_main_o = new namespace_LLR::llr_hmc(s_llrparams_in);

    // printing the structure to see how it is constructed.
    p_llr_hmc_main_o->init_robbins_monro(s_llrparams_in);
    p_llr_hmc_main_o->print_s_llrparams(s_llrparams_in);

    // Start of the main commands.
    bool with_llr = false;
    std::cout<<"<---- with_llr (initialized)     ---->: " << with_llr << std::endl;
    std::cout<<"<---- llr_config                 ---->: " << llr_config << std::endl;
#if defined(llr_config)
    with_llr = true;
#endif
    std::cout<<"<---- with_llr (llr_config)      ---->: " << with_llr << std::endl;
    std::cout<<"<---- Sp2n_config                ---->: " << Sp2n_config << std::endl;

    // if --enable--LLR switch is activated
    if (with_llr) {
        // SU(N) representation
        typedef Grid::GenericHMCRunnerLLR <Grid::MinimumNorm2> HMCWrapperLLR;
        HMCWrapperLLR TheHMC;
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
        typedef Grid::PlaquetteMod<HMCWrapperLLR::ImplPolicy> PlaqObs;
        typedef Grid::TopologicalChargeMod<HMCWrapperLLR::ImplPolicy> QObs;
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

        //Grid::SpWilsonGaugeActionR Waction(beta);

        //typedef Grid::LLRGaugeAction<Grid::SpWilsonGaugeActionR , Grid::PeriodicGimplR> LLRGaugeActionR;
        typedef Grid::LLRGaugeAction<Grid::WilsonGaugeActionR , Grid::PeriodicGimplR> LLRGaugeActionR;

        LLRGaugeActionR LLRaction(s_llrparams_in, beta);

        Grid::ActionLevel<HMCWrapperLLR::Field> Level1(1);
        Level1.push_back(&LLRaction);
        //Level1.push_back(WGMod.getPtr());
        TheHMC.TheAction.push_back(Level1);
        /////////////////////////////////////////////////////////////

        // HMC parameters are serialisable
        TheHMC.Parameters.MD.MDsteps = 20;
        TheHMC.Parameters.MD.trajL   = 1.0;

        TheHMC.Run();  // no smearing


    } /* [end-if] with_llr */

    // Finalizing the Grid library
    Grid::Grid_finalize();
    // End statement
    std::cout<<"<---- End LLR_HMC_SUn_WilsonGauge.cc ---->" <<std::endl;

} /* end of main LLR-HMC-WilsonGauge.cc */
