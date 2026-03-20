/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/gauge/Gauge_aggregate.h

Copyright (C) 2015

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
#ifndef GRID_QCD_GAUGE_H
#define GRID_QCD_GAUGE_H

/*! \dir
 * This directory contains implementations of a variety of discretisations of a non-Abelian gauge field.
 * To get started,
 * see Gauge.h.
 */

/*! \file
 * # Lattice gauge fields
 *
 * This header `#include`s all relevant headers for gauge actions.
 */

#include <Grid/qcd/action/gauge/GaugeImplementations.h>
#include <Grid/qcd/utils/WilsonLoops.h>
#include <Grid/qcd/action/gauge/WilsonGaugeAction.h>
#include <Grid/qcd/action/gauge/PlaqPlusRectangleAction.h>

/// \cond DO_NOT_DOCUMENT
NAMESPACE_BEGIN(Grid);
/// \endcond

/*! \name Wilson gauge action
* See WilsonGaugeAction.h for more details. */
/*! @{ */
/*! Wilson gauge action for \f$SU(N)\f$ or \f$Sp(2n)\f$, in default precision */
typedef WilsonGaugeAction<PeriodicGimplR>          WilsonGaugeActionR;
/*! Single-precision Wilson gauge action for \f$SU(N)\f$ or \f$Sp(2n)\f$ */
typedef WilsonGaugeAction<PeriodicGimplF>          WilsonGaugeActionF;
/*! Double-precision Wilson gauge action for \f$SU(N)\f$ or \f$Sp(2n)\f$ */
typedef WilsonGaugeAction<PeriodicGimplD>          WilsonGaugeActionD;
/*! Wilson gauge action for \f$Sp(2n)\f$, in default precision */
typedef WilsonGaugeAction<SpPeriodicGimplR>        SpWilsonGaugeActionR;
/*! Single-precision Wilson gauge action for \f$Sp(2n)\f$ */
typedef WilsonGaugeAction<SpPeriodicGimplF>        SpWilsonGaugeActionF;
/*! Double-precision Wilson gauge action for \f$Sp(2n)\f$ */
typedef WilsonGaugeAction<SpPeriodicGimplD>        SpWilsonGaugeActionD;
/*! @} */

/*! \name Improved gauge actions
 * See PlaqPlusRectangleAction.h for more details. */
/*! @{ */
/*! General plaquette plus rectangle gauge action for \f$SU(N)\f$ or \f$Sp(2n)\f$, in default precision */
typedef PlaqPlusRectangleAction<PeriodicGimplR>    PlaqPlusRectangleActionR;
/*! Single-precision, general plaquette plus rectangle gauge action for \f$SU(N)\f$ or \f$Sp(2n)\f$ */
typedef PlaqPlusRectangleAction<PeriodicGimplF>    PlaqPlusRectangleActionF;
/*! Double-precision, general plaquette plus rectangle gauge action for \f$SU(N)\f$ or \f$Sp(2n)\f$ */
typedef PlaqPlusRectangleAction<PeriodicGimplD>    PlaqPlusRectangleActionD;
/*! Iwasaki gauge action for \f$SU(N)\f$ or \f$Sp(2n)\f$, in default precision */
typedef IwasakiGaugeAction<PeriodicGimplR>         IwasakiGaugeActionR;
/*! Single-precision Iwasaki gauge action for \f$SU(N)\f$ or \f$Sp(2n)\f$ */
typedef IwasakiGaugeAction<PeriodicGimplF>         IwasakiGaugeActionF;
/*! Double-precision Iwasaki gauge action for \f$SU(N)\f$ or \f$Sp(2n)\f$ */
typedef IwasakiGaugeAction<PeriodicGimplD>         IwasakiGaugeActionD;
/*! Symanzik gauge action for \f$SU(N)\f$ or \f$Sp(2n)\f$, in default precision */
typedef SymanzikGaugeAction<PeriodicGimplR>        SymanzikGaugeActionR;
/*! Single-precision Symanzik gauge action for \f$SU(N)\f$ or \f$Sp(2n)\f$ */
typedef SymanzikGaugeAction<PeriodicGimplF>        SymanzikGaugeActionF;
/*! Double-precision Symanzik gauge action for \f$SU(N)\f$ or \f$Sp(2n)\f$ */
typedef SymanzikGaugeAction<PeriodicGimplD>        SymanzikGaugeActionD;
/*! @} */

/*! \name Gauge actions with conjugated boundary conditions */
/*! @{ */
/*! Wilson gauge action for \f$SU(N)\f$ or \f$Sp(2n)\f$, in default precision, with conjugated boundary conditions */
typedef WilsonGaugeAction<ConjugateGimplR>          ConjugateWilsonGaugeActionR;
/*! Single-precision Wilson gauge action for \f$SU(N)\f$ or \f$Sp(2n)\f$, with conjugated boundary conditions */
typedef WilsonGaugeAction<ConjugateGimplF>          ConjugateWilsonGaugeActionF;
/*! Double-precision Wilson gauge action for \f$SU(N)\f$ or \f$Sp(2n)\f$, with conjugated boundary conditions */
typedef WilsonGaugeAction<ConjugateGimplD>          ConjugateWilsonGaugeActionD;

/*!  General plaquette plus rectangle gauge action for \f$SU(N)\f$ or \f$Sp(2n)\f$, in default precision, with conjugated boundary conditions */
typedef PlaqPlusRectangleAction<ConjugateGimplR>    ConjugatePlaqPlusRectangleActionR;
/*!  Single-precision, general plaquette plus rectangle gauge action for \f$SU(N)\f$ or \f$Sp(2n)\f$, with conjugated boundary conditions */
typedef PlaqPlusRectangleAction<ConjugateGimplF>    ConjugatePlaqPlusRectangleActionF;
/*!  Double-precision, general plaquette plus rectangle gauge action for \f$SU(N)\f$ or \f$Sp(2n)\f$, with conjugated boundary conditions */
typedef PlaqPlusRectangleAction<ConjugateGimplD>    ConjugatePlaqPlusRectangleActionD;
/*!  Iwasaki gauge action for \f$SU(N)\f$ or \f$Sp(2n)\f$, in default precision, with conjugated boundary conditions */
typedef IwasakiGaugeAction<ConjugateGimplR>         ConjugateIwasakiGaugeActionR;
/*!  Single-precision Iwasaki gauge action for \f$SU(N)\f$ or \f$Sp(2n)\f$, with conjugated boundary conditions */
typedef IwasakiGaugeAction<ConjugateGimplF>         ConjugateIwasakiGaugeActionF;
/*!  Double-precision Iwasaki gauge action for \f$SU(N)\f$ or \f$Sp(2n)\f$, with conjugated boundary conditions */
typedef IwasakiGaugeAction<ConjugateGimplD>         ConjugateIwasakiGaugeActionD;
/*!  Symanzik gauge action for \f$SU(N)\f$ or \f$Sp(2n)\f$, in default precision, with conjugated boundary conditions */
typedef SymanzikGaugeAction<ConjugateGimplR>        ConjugateSymanzikGaugeActionR;
/*!  Single-precision Symanzik gauge action for \f$SU(N)\f$ or \f$Sp(2n)\f$, with conjugated boundary conditions */
typedef SymanzikGaugeAction<ConjugateGimplF>        ConjugateSymanzikGaugeActionF;
/*!  Double-precision Symanzik gauge action for \f$SU(N)\f$ or \f$Sp(2n)\f$, with conjugated boundary conditions */
typedef SymanzikGaugeAction<ConjugateGimplD>        ConjugateSymanzikGaugeActionD;
/*! @} */

/// \cond DO_NOT_DOCUMENT
NAMESPACE_END(Grid);
/// \endcond

#endif
