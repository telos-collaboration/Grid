    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/Actions.h

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Peter Boyle <peterboyle@Peters-MacBook-Pro-2.local>
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

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
#ifndef GRID_QCD_ACTION_H
#define GRID_QCD_ACTION_H

/*! \dir
 * This directory holds implementations of lattice quantum fields that Grid can work with.
 *
 * To construct a particular action of interest,
 * one should use instances of the relevant ``Action`` types for each term contributing to the action.
 * For example,
 * for 2+1 flavour QCD,
 * one would need one gauge action instance (e.g. ``WilsonGaugeAction``),
 * one a two-flavour pseudofermion action (e.g. ``TwoFlavourPseudoFermionAction``)
 * for a given fermion action (e.g. ``MobiusFermionD``),
 * and a one-flavour pseudofermion action (e.g. ``OneFlavourRationalPseudoFermionAction``)
 * for the same fermion action.
 *
 * - The gauge/ directory contains actions for gauge fields
 * - The fermion/ directory contains fermionic actions
 * - The pseudofermion/ directory contains actions that wrap fermion fields for use in HMC;
 *   for example,
 *   as a single flavour (RHMC),
 *   two flavours (HMC),
 *   or a ratio (Hasenbusch acceleration, Pauli-Villars fields).
 * - The scalar/ directory contains actions for scalar fields and interactions
 * - The filter/ directory contains tools to manipulate
 *   the distribution of conjugate momenta as part of the update/refresh
 *
 * All actions implement the interface defined in ActionBase.h.
 */

////////////////////////////////////////////
// Abstract base interface
////////////////////////////////////////////
#include <Grid/qcd/action/ActionCore.h>
/// \cond DO_NOT_DOCUMENT
NAMESPACE_CHECK(ActionCore);
/// \endcond
////////////////////////////////////////////////////////////////////////
// Fermion actions; prevent coupling fermion.cc files to other headers
////////////////////////////////////////////////////////////////////////
#include <Grid/qcd/action/fermion/FermionCore.h>
/// \cond DO_NOT_DOCUMENT
NAMESPACE_CHECK(FermionCore);
/// \endcond
#include <Grid/qcd/action/fermion/Fermion.h>
/// \cond DO_NOT_DOCUMENT
NAMESPACE_CHECK(Fermion);
/// \endcond
////////////////////////////////////////
// Pseudo fermion combinations for HMC
////////////////////////////////////////
#include <Grid/qcd/action/pseudofermion/PseudoFermion.h>
/// \cond DO_NOT_DOCUMENT
NAMESPACE_CHECK(PseudoFermion);
/// \endcond

#endif
