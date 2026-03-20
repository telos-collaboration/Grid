/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/gauge/PlaqPlusRectangleAction.h

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
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
#ifndef QCD_PLAQ_PLUS_RECTANGLE_ACTION_H
#define QCD_PLAQ_PLUS_RECTANGLE_ACTION_H

/// \cond DO_NOT_DOCUMENT
NAMESPACE_BEGIN(Grid);
/// \endcond

/*! @file
 * This file contains improved gauge actions,
 * expected to give more rapid convergence to the continuum limit than the WilsonGaugeAction.
 */
    
template<class Gimpl>
class PlaqPlusRectangleAction : public Action<typename Gimpl::GaugeField>
/*! @brief Plaquette plus rectangle actions
 *
 * The simplest way to add additional contributions to a gauge action
 * is in addition to the elementary plaquette used in the WilsonGaugeAction,
 * to additionally incorporate the next-smallest element,
 * the product of two adjacent plaquettes forming a rectangle in a common plane.
 *
 * This may be added with different coefficients for the two contributions,
 * called \f$b\f$ and \f$c\f$ here
 * (but frequently called \f$c_0\f$ and \f$c_1\f$ in other sources);
 * specific choices of these
 * give rise to the different named actions inheriting from this base class.
 */
{
public:
  /*! @brief All attributes of the gauge implementation are available from the action.
   * See GaugeImplTypes.h
   */
  INHERIT_GIMPL_TYPES(Gimpl);

  using Action<GaugeField>::S;
  using Action<GaugeField>::Sinitial;
  using Action<GaugeField>::deriv;
  using Action<GaugeField>::refresh;

private:
  RealD c_plaq;  ///< @brief The coefficient \f$b=c_0\f$ of the plaquette contribution
  RealD c_rect;  ///< @brief The coefficient \f$c=c_1\f$ of the rectangle contribution
  typename WilsonLoops<Gimpl>::StapleAndRectStapleAllWorkspace workspace;
public:
  /*! @brief Initialise the action with the relevant coefficients.
   *
   * @param b: The coefficient of the elementary plaquette contribution
   * @param c: The coefficient of the \f$1\times2\f$ rectangle contribution
   */
  PlaqPlusRectangleAction(RealD b,RealD c): c_plaq(b),c_rect(c){};

  virtual std::string action_name(){return "PlaqPlusRectangleAction";}
      
  /*! @brief Gauge fields do not have pseudofermions, so this is a no-op */
  virtual void refresh(const GaugeField &U, GridSerialRNG &sRNG, GridParallelRNG& pRNG) {};
      
  virtual std::string LogParameters(){
    std::stringstream sstream;
    sstream << GridLogMessage << "["<<action_name() <<"] c_plaq: " << c_plaq << std::endl;
    sstream << GridLogMessage << "["<<action_name() <<"] c_rect: " << c_rect << std::endl;
    return sstream.str();
  }

  /*! @brief The plaquette plus rectangle gauge action itself.
   *
   * \f[S_g = - \frac{\beta}{3} \left[b \sum_{x;\mu<\nu} P[U]_{x,\mu\nu} + c \sum_{x;\mu\ne\nu} R[U]_{x,\mu\nu}\right]\f]
   * where \f$P[U]\f$ is the mean of the real part of the trace of
   * the path ordered product of links around the \f$1\times1\f$ plaquette,
   * and \f$R[U]\f$ is the equivalent quantity for the \f$1\times2\f$ rectangle.
   *
   * @param U: The gauge field on which to compute the action. */
  virtual RealD S(const GaugeField &U) {
    RealD vol = U.Grid()->gSites();

    RealD plaq = WilsonLoops<Gimpl>::avgPlaquette(U);
    RealD rect = WilsonLoops<Gimpl>::avgRectangle(U);

    RealD action=c_plaq*(1.0 -plaq)*(Nd*(Nd-1.0))*vol*0.5
      +c_rect*(1.0 -rect)*(Nd*(Nd-1.0))*vol;

    return action;
  };

  /*! @brief The derivative of any PlaqPlusRectangleGaugeAction subclass.
   *
   * @param U: The gauge field on which to compute the derivative
   * @param dSdU: Output field into which to write the derivative
   */
  virtual void deriv(const GaugeField &U, GaugeField &dSdU) {
    //extend Ta to include Lorentz indexes
    RealD factor_p = c_plaq/RealD(Nc)*0.5;
    RealD factor_r = c_rect/RealD(Nc)*0.5;

    GridBase *grid = U.Grid();

    std::vector<GaugeLinkField> Umu (Nd,grid);
    for(int mu=0;mu<Nd;mu++){
      Umu[mu] = PeekIndex<LorentzIndex>(U,mu);
    }
    std::vector<GaugeLinkField> RectStaple(Nd,grid), Staple(Nd,grid);
    WilsonLoops<Gimpl>::StapleAndRectStapleAll(Staple, RectStaple, Umu, workspace);

    GaugeLinkField dSdU_mu(grid);
    GaugeLinkField staple(grid);

    for (int mu=0; mu < Nd; mu++){
      dSdU_mu = Ta(Umu[mu]*Staple[mu])*factor_p;
      dSdU_mu = dSdU_mu + Ta(Umu[mu]*RectStaple[mu])*factor_r;

      PokeIndex<LorentzIndex>(dSdU, dSdU_mu, mu);
    }

  };

};

template<class Gimpl>
class RBCGaugeAction : public PlaqPlusRectangleAction<Gimpl>
/*! @brief Convenience class for common physically defined cases.
 *
 * This defines a special case of PlaqPlusRectangleAction
 * for which the coefficient of the plaquette is \f$\beta(1-8c_1)\f$,
 * and the coefficient of the rectangle is \f$\beta c_1\f$.
 *
 * RBC c1 parameterisation was not really invented by RBC,
 * but they are not aware of a good reference.
 * This name may be changed if prior use of the
 * parameterisation in terms of plaquette coefficients is made known;
 * pull requests with references are welcome.
 */
{
public:
  INHERIT_GIMPL_TYPES(Gimpl);
  RBCGaugeAction(RealD beta,RealD c1) : PlaqPlusRectangleAction<Gimpl>(beta*(1.0-8.0*c1), beta*c1) {};
  virtual std::string action_name(){return "RBCGaugeAction";}
};

template<class Gimpl>
class IwasakiGaugeAction : public RBCGaugeAction<Gimpl>
  /*! @brief The Iwasaki gauge action, model IM11 in arXiv:1111.7054.
   *
   * Iwasaki introduced this action as one of a number of candidates
   * in a 1983 report,
   * later uploaded to arXiv as https://arxiv.org/abs/1111.7054
   */
{
public:
  INHERIT_GIMPL_TYPES(Gimpl);
  IwasakiGaugeAction(RealD beta) : RBCGaugeAction<Gimpl>(beta,-0.331) {};
  virtual std::string action_name(){return "IwasakiGaugeAction";}
};

template<class Gimpl>
class SymanzikGaugeAction : public RBCGaugeAction<Gimpl>
/*! @brief The tree-level Symanzik improved gauge action,
 * introduced in Nucl.Phys.B 236 (1984) 397
 *
 * See DOI https://doi.org/10.1016/0550-3213(84)90543-1; Eq. (1.3)
 */
{
public:
  INHERIT_GIMPL_TYPES(Gimpl);
  SymanzikGaugeAction(RealD beta) : RBCGaugeAction<Gimpl>(beta,-1.0/12.0) {};
  virtual std::string action_name(){return "SymanzikGaugeAction";}
};

template<class Gimpl>
class DBW2GaugeAction : public RBCGaugeAction<Gimpl>
/*! @brief The "doubly blocked from Wilson action in two coupling space"
 * (DBW2 for short),
 * introduced by Takaishi in Phys.Rev.D 54(1996) 1050.
 *
 * See the original paper at https://doi.org/10.1103/PhysRevD.54.1050,
 * the explanation of the DBW2 acronym at https://doi.org/10.1016/S0550-3213(00)00145-0.
 * The coefficients in the form used here
 * were not explicitly written down in either of these;
 * see Iwasaki's 1983 report https://arxiv.org/abs/1111.7054
 * (underneath Eq. (3))
 * for this.
 */
{
public:
  INHERIT_GIMPL_TYPES(Gimpl);
  DBW2GaugeAction(RealD beta) : RBCGaugeAction<Gimpl>(beta,-1.4067) {};
  virtual std::string action_name(){return "DBW2GaugeAction";}
};

/// \cond DO_NOT_DOCUMENT
NAMESPACE_END(Grid);
/// \endcond

#endif
