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

/////////////////////////////////////////////////////////////
/// Main
/////////////////////////////////////////////////////////////
int main(int argc, char **argv) {
  
  std::cout << Grid::GridLogLLR <<"<---- Start Test_llr_contrained_action.cc ---->"<<std::endl;
  // Initializing Grid library environment.
  Grid::Grid_init(&argc, &argv);
  Grid::GridLogLayout();
  
  // Instantiating the inout parameter structure as a pointer.
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
  
  // Finalising Grid environment.
  Grid::Grid_finalize();
}
