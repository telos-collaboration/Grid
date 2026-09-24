/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ConstrainedAction.h

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
#include <sstream>

#include <Grid/qcd/action/ActionBase.h>

NAMESPACE_BEGIN(Grid);

/* Gaussian-constrained action for the LLR algorithm.
 *
 * Implements the constrained action
 *   (beta) S[U] -> aS[U] + (S[U]-S_0)^2/(2 * sigma^2).
 * given a gauge action S.
 * Assumes that beta can be freely set to 1 in order to perform
 * the direct replacement of S[U].
 */
struct ConstrainedActionParameters
{
  RealD a;     // LLR parameter - tuned to d ln[rho(S)] / dS at S0
  RealD S0;    // Constraint centre
  RealD sigma; // Gaussian width
};

template <class WrappedAction>
class ConstrainedAction : public Action<typename WrappedAction::GaugeField>
{
public:
  using GaugeField = typename WrappedAction::GaugeField;

  using Action<GaugeField>::S;
  using Action<GaugeField>::Sinitial;
  using Action<GaugeField>::deriv;
  using Action<GaugeField>::refresh;

  ConstrainedAction(WrappedAction &wrapped, ConstrainedActionParameters parameters)
      : wrapped_(wrapped), parameters_(parameters)
  {
    this->is_smeared = wrapped_.is_smeared;
  }

  virtual void refresh(const GaugeField &U, GridSerialRNG &sRNG, GridParallelRNG &pRNG)
  {
    wrapped_.refresh(U, sRNG, pRNG);
  }

  virtual RealD S(const GaugeField &U)
  {
    return constrained_value(wrapped_.S(U));
  }

  virtual RealD Sinitial(const GaugeField &U)
  {
    return constrained_value(wrapped_.Sinitial(U));
  }

  RealD Sunconstrained(const GaugeField &U)
  {
    return wrapped_.S(U);
  }

  virtual void deriv(const GaugeField &U, GaugeField &force)
  {
    RealD base_action = wrapped_.S(U);
    wrapped_.deriv(U, force);

    RealD scale = parameters_.a + (base_action - parameters_.S0) / (parameters_.sigma * parameters_.sigma);
    force *= scale;
  }

  virtual std::string action_name()
  {
    return "ConstrainedAction<" + wrapped_.action_name() + ">";
  }

  virtual std::string LogParameters()
  {
    std::stringstream sstream;
    sstream << GridLogMessage << "[" << action_name() << "] a:     " << parameters_.a << std::endl;
    sstream << GridLogMessage << "[" << action_name() << "] S0:    " << parameters_.S0 << std::endl;
    sstream << GridLogMessage << "[" << action_name() << "] sigma: " << parameters_.sigma << std::endl;
    sstream << wrapped_.LogParameters();
    return sstream.str();
  }

  const ConstrainedActionParameters &parameters() const
  {
    return parameters_;
  }

  void set_parameters(ConstrainedActionParameters parameters)
  {
    parameters_ = parameters;
  }

  void set_a(RealD a)
  {
    parameters_.a = a;
  }

  void set_S0(RealD S0)
  {
    parameters_.S0 = S0;
  }

  void set_sigma(RealD sigma)
  {
    parameters_.sigma = sigma;
  }

private:
  WrappedAction &wrapped_;
  ConstrainedActionParameters parameters_;

  RealD constrained_value(RealD base_action) const
  {
    RealD displacement = base_action - parameters_.S0;
    return (parameters_.a * base_action) + displacement * displacement / (2.0 * parameters_.sigma * parameters_.sigma);
  }
};

NAMESPACE_END(Grid);

