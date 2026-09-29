/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: RobbinsMonroSolver.h

Copyright (C) 2026

Author: Ryan Hill <Ryan.Hill@ed.ac.uk>

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

#include <cassert>

#include <Grid/qcd/action/gauge/ConstrainedAction.h>

NAMESPACE_BEGIN(Grid);

enum class RobbinsMonroPhase
{
  Thermalising,
  Accumulating
};

struct RobbinsMonroParameters
{
  int initial_thermalisation_trajectories;
  int rethermalisation_trajectories;
  int trajectories_per_update;
  RealD gain;

  RobbinsMonroParameters(int initial_thermalisation_trajectories = 0,
                         int rethermalisation_trajectories_ = 0,
                         int trajectories_per_update_ = 1,
                         RealD gain_ = 1.0)
      : initial_thermalisation_trajectories(initial_thermalisation_trajectories),
        rethermalisation_trajectories(rethermalisation_trajectories_),
        trajectories_per_update(trajectories_per_update_),
        gain(gain_)
  {}
};

struct RobbinsMonroUpdate
{
  int trajectory;
  int iteration;
  RealD mean_action;
  RealD residual;
  RealD previous_a;
  RealD updated_a;
};

struct RobbinsMonroStatus
{
  RobbinsMonroPhase phase;
  int iteration;
  int trajectories_remaining_in_phase;
  RealD accumulated_action;
  RobbinsMonroUpdate last_update;
};

template <class ConstrainedActionType>
class RobbinsMonroSolver {
public:
  typedef ConstrainedActionType ActionType;
  typedef typename ActionType::GaugeField Field;

  RobbinsMonroSolver(ActionType &action, RobbinsMonroParameters parameters)
      : action_(action)
      , parameters_(parameters)
      , status_{
          .phase=RobbinsMonroPhase::Accumulating,
          .iteration=1,
          .trajectories_remaining_in_phase=parameters.trajectories_per_update,
          .accumulated_action=0
        }
  {
    if (parameters_.initial_thermalisation_trajectories != 0)
    {
      status_.phase = RobbinsMonroPhase::Thermalising;
      status_.trajectories_remaining_in_phase = parameters_.initial_thermalisation_trajectories;
    }
  }

  void record_configuration(int trajectory, Field &U)
  {
    status_.trajectories_remaining_in_phase--;
    if (status_.phase == RobbinsMonroPhase::Accumulating)
    {
      std::cout << GridLogMessage << "RM Phase: Accumulating (remaining trajectories: " << status_.trajectories_remaining_in_phase << ")" << std::endl;
      accumulate(trajectory, U);
    }
    else // Thermalising
    {
      std::cout << GridLogMessage << "RM Phase: Thermalising (remaining trajectories: " << status_.trajectories_remaining_in_phase << ")" << std::endl;
      thermalise();
    }
  }

  void record_configuration(int trajectory, ConfigurationBase<Field> &configuration)
  {
    Field &U = configuration.get_U(action_.is_smeared);
    record_configuration(trajectory, U);
  }

  // Accessors for replica-swapping
  const RobbinsMonroParameters &parameters() const { return parameters_; }
  void restore_state(const RobbinsMonroStatus &state) { status_ = state; }
  RobbinsMonroStatus status() const { return status_; }
public:
  ActionType &action_;
private:
  RobbinsMonroParameters parameters_;
  RobbinsMonroStatus status_;

  void accumulate(int trajectory, Field &U)
  {
    // Accumulation step
    status_.accumulated_action += action_.Sunconstrained(U);
    if (status_.trajectories_remaining_in_phase > 0) // Ready to update? If not, return to continue sampling.
    {
      return;
    }

    // After we've gathered enough samples, perform the update of a.
    RobbinsMonroUpdate update;
    update.trajectory = trajectory;
    update.iteration = status_.iteration;
    update.mean_action = status_.accumulated_action / parameters_.trajectories_per_update;
    update.residual = update.mean_action - action_.parameters().S0;
    update.previous_a = action_.parameters().a;
    update.updated_a = update.previous_a + parameters_.gain * update.residual / (action_.parameters().sigma * action_.parameters().sigma * status_.iteration);

    action_.set_a(update.updated_a);
    status_.last_update = update;
    ++status_.iteration;
    status_.accumulated_action = 0.0;

    // Swap over the thermalisation steps to adjust for new value of a.
    status_.phase = RobbinsMonroPhase::Thermalising;
    status_.trajectories_remaining_in_phase = parameters_.rethermalisation_trajectories; 
  }

  void thermalise()
  {
    // Only role here is to switch state to 'Accumulating' when # thermalisation steps have passed.
    if (status_.trajectories_remaining_in_phase <= 0)
    {
      status_.phase = RobbinsMonroPhase::Accumulating;
      status_.trajectories_remaining_in_phase = parameters_.trajectories_per_update; 
    }
  }
};

NAMESPACE_END(Grid);

