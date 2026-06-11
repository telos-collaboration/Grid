/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/observables/hmc_observable.h

Copyright (C) 2017

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

#pragma once

NAMESPACE_BEGIN(Grid);

template <class Field>
class HmcObservable 
/*! @brief HMC Observable class,
 *  HMC observables inherit the method 
 *  TrajectoryComplete from here, which
 *  is used to run the computation of
 *  on-the-fly observables.
 *
 *  Expects:
 *   an integer traj : trajectory number
 *   a reference to ConfigurationBase : the gauge field configuration
 *   a reference to a serial RNG
 *   a reference to a parallel RNG
 */
{
 public:
  virtual void TrajectoryComplete(int traj,
                                  ConfigurationBase<Field> &SmartConfig,
                                  GridSerialRNG &sRNG,
                                  GridParallelRNG &pRNG)
  {
    TrajectoryComplete(traj,SmartConfig.get_U(false),sRNG,pRNG); // Unsmeared observable
  };
  virtual void TrajectoryComplete(int traj,
                                  Field &U,
                                  GridSerialRNG &sRNG,
                                  GridParallelRNG &pRNG) = 0;

  // allow backward compatibility with current observables that do not have 
  // acceptance argument
  virtual void TrajectoryComplete(int traj,
                                  ConfigurationBase<Field> &SmartConfig,
                                  GridSerialRNG &sRNG,
                                  GridParallelRNG &pRNG, 
                                  bool accept)
  {
    TrajectoryComplete(traj,SmartConfig.get_U(false),sRNG,pRNG,accept); // Unsmeared observable
  };
  virtual void TrajectoryComplete(int traj, 
                                  Field &U,
                                  GridSerialRNG &sRNG,
                                  GridParallelRNG &pRNG,
                                  bool accept) 
  {
    TrajectoryComplete(traj, U, sRNG, pRNG);
  };
};

NAMESPACE_END(Grid);

NAMESPACE_CHECK(HmcObs);
#include "plaquette.h"
NAMESPACE_CHECK(Plaq);
#include "topological_charge.h"
NAMESPACE_CHECK(Topo);
#include "polyakov_loop.h"
NAMESPACE_CHECK(Polyakov);


