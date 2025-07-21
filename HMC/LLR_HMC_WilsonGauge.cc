/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./HMC/LLR-HMC-WilsonGauge.cc

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

int main(int argc, char **argv) {
    std::cout<<"Start LLR-HMC-WilsonGauge.cc" <<std::endl;
    // Initializing Grid library
    Grid::Grid_init(&argc, &argv);
    Grid::GridLogLayout();

    typedef Grid::GenericHMCRunnerLLR<Grid::MinimumNorm2> HMCWrapperLLR;  // Uses the default minimum norm
    HMCWrapperLLR TheHMC;

    // Grid from the command line
    TheHMC.Resources.AddFourDimGrid("gauge");
    // Possible to create the module by hand
    // hardcoding parameters or using a Reader

    // Checking the parameters in the used, we will use the same as the
    // Standard wilson ones.
    Grid::CheckpointerParameters CPparams;
    CPparams.config_prefix = "ckpoint_lat";
    CPparams.rng_prefix = "ckpoint_rng";
    CPparams.saveInterval = 1;
    CPparams.format = "IEEE64BIG";

    TheHMC.Resources.LoadNerscCheckpointer(CPparams);

    Grid::RNGModuleParameters RNGpar;
    RNGpar.serial_seeds = "1 2 3 4 5";
    RNGpar.parallel_seeds = "6 7 8 9 10";
    TheHMC.Resources.SetRNGSeeds(RNGpar);



    // Finalizing the Grid library
    Grid::Grid_finalize();
    // End statement
    std::cout<<"End LLR-HMC-WilsonGauge.cc" <<std::endl;

} /* end of main LLR-HMC-WilsonGauge.cc */
