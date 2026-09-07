/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./tests/core/Test_ConstrainedAction.cc

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

    See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */
#include <cassert>
#include <cmath>

#include <Grid/Grid.h>

using namespace Grid;

namespace {

bool isclose(RealD actual, RealD expected, RealD relative_tolerance)
{
  return std::abs(actual - expected) <= relative_tolerance * (1.0 + std::abs(expected));
}

void assert_force_matches(const LatticeGaugeField &actual,
                          const LatticeGaugeField &expected,
                          RealD relative_tolerance)
{
  LatticeGaugeField difference(actual.Grid());
  difference = actual - expected;
  assert(norm2(difference) <= relative_tolerance * (1.0 + norm2(expected)));
}

RealD constrained_value(RealD base_action,
                        const ConstrainedActionParameters &parameters)
{
  RealD displacement = base_action - parameters.S0;
  return parameters.a * base_action + displacement * displacement / (2.0 * parameters.sigma * parameters.sigma);
}

void assert_action_and_force(WilsonGaugeActionR &wrapped,
                             ConstrainedAction<WilsonGaugeActionR> &constrained,
                             const LatticeGaugeField &U,
                             const ConstrainedActionParameters &parameters)
{
  RealD base_action = wrapped.S(U);

  // Check the constrained action value is the stated LLR potential.
  std::cout << GridLogMessage << "  Constrained action value: " << constrained.S(U) << std::endl;
  std::cout << GridLogMessage << "  Expected:                 " << constrained_value(base_action, parameters) << std::endl;
  assert(isclose(constrained.S(U),
                 constrained_value(base_action, parameters),
                 1.0e-12));

  // Check the initial action is transformed by the same potential.
  std::cout << GridLogMessage << "  Constrained initial value: " << constrained.Sinitial(U) << std::endl;
  std::cout << GridLogMessage << "  Expected:                  " << constrained_value(wrapped.Sinitial(U), parameters) << std::endl;
  assert(isclose(constrained.Sinitial(U),
                 constrained_value(wrapped.Sinitial(U), parameters),
                 1.0e-12));

  LatticeGaugeField base_force(U.Grid());
  LatticeGaugeField constrained_force(U.Grid());
  LatticeGaugeField expected_force(U.Grid());

  wrapped.deriv(U, base_force);
  constrained.deriv(U, constrained_force);

  RealD scale = parameters.a + (base_action - parameters.S0) / (parameters.sigma * parameters.sigma);
  expected_force = scale * base_force;

  // Check the constrained force obeys the LLR chain rule.
  assert_force_matches(constrained_force, expected_force, 1.0e-12);
}

}  // namespace

int main(int argc, char **argv)
{
  Grid_init(&argc, &argv);

  Coordinate lattice_size({4, 4, 4, 4});
  Coordinate simd_layout = GridDefaultSimd(4, vComplex::Nsimd());
  Coordinate mpi_layout = GridDefaultMpi();
  GridCartesian grid(lattice_size, simd_layout, mpi_layout);

  GridParallelRNG parallel_rng(&grid);
  parallel_rng.SeedFixedIntegers({1, 2, 3, 4});

  LatticeGaugeField U(&grid);
  SU<Nc>::HotConfiguration(parallel_rng, U);

  WilsonGaugeActionR wrapped(1.0);
  ConstrainedActionParameters parameters{1.7, wrapped.S(U) + 0.25, 2.5};
  ConstrainedAction<WilsonGaugeActionR> constrained(wrapped, parameters);

  std::cout << GridLogMessage << "--- Test constrained values from constructor ---" << std::endl;
  assert(!constrained.is_smeared);
  assert_action_and_force(wrapped, constrained, U, parameters);

  std::cout << GridLogMessage << "--- Test set individual parameters ---" << std::endl;
  parameters.a = 2.2;
  constrained.set_a(parameters.a);
  std::cout << GridLogMessage << "  Expected a: " << parameters.a << std::endl;
  std::cout << GridLogMessage << "  Actual a:   " << constrained.parameters().a << std::endl;
  assert_action_and_force(wrapped, constrained, U, parameters);

  parameters.S0 -= 0.4;
  constrained.set_S0(parameters.S0);
  std::cout << GridLogMessage << "  Expected S0: " << parameters.S0 << std::endl;
  std::cout << GridLogMessage << "  Actual S0:   " << constrained.parameters().S0 << std::endl;
  assert_action_and_force(wrapped, constrained, U, parameters);

  parameters.sigma = 1.75;
  constrained.set_sigma(parameters.sigma);
  std::cout << GridLogMessage << "  Expected Sigma: " << parameters.sigma << std::endl;
  std::cout << GridLogMessage << "  Actual Sigma:   " << constrained.parameters().sigma << std::endl;
  assert_action_and_force(wrapped, constrained, U, parameters);

  std::cout << GridLogMessage << "--- Test set all parameters ---" << std::endl;
  parameters = ConstrainedActionParameters{0.8, wrapped.S(U) - 0.5, 1.5};
  constrained.set_parameters(parameters);
  std::cout << GridLogMessage << "  Expected a: " << parameters.a << std::endl;
  std::cout << GridLogMessage << "  Actual a:   " << constrained.parameters().a << std::endl;
  std::cout << GridLogMessage << "  Expected S0: " << parameters.S0 << std::endl;
  std::cout << GridLogMessage << "  Actual S0:   " << constrained.parameters().S0 << std::endl;
  std::cout << GridLogMessage << "  Expected Sigma: " << parameters.sigma << std::endl;
  std::cout << GridLogMessage << "  Actual Sigma:   " << constrained.parameters().sigma << std::endl;
  assert(constrained.parameters().a == parameters.a);
  assert(constrained.parameters().S0 == parameters.S0);
  assert(constrained.parameters().sigma == parameters.sigma);
  assert_action_and_force(wrapped, constrained, U, parameters);

  // Check 6: the wrapper propagates smearing and refresh.
  wrapped.is_smeared = true;
  ConstrainedAction<WilsonGaugeActionR> smeared_constrained(wrapped, parameters);
  std::cout << GridLogMessage << "Check setting smear works:" << std::endl;
  std::cout << GridLogMessage << "  " << smeared_constrained.is_smeared << std::endl;
  assert(smeared_constrained.is_smeared);
  wrapped.is_smeared = false;

  GridSerialRNG serial_rng;
  serial_rng.SeedFixedIntegers({5, 6, 7, 8});
  constrained.refresh(U, serial_rng, parallel_rng); // no-op for pure gauge action
  assert_action_and_force(wrapped, constrained, U, parameters);

  Grid_finalize();
  return 0;
}
