/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/ActionBase.h

Copyright (C) 2015-2016

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: neo <cossu@post.kek.jp>
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

#ifndef ACTION_BASE_H
#define ACTION_BASE_H

/// \cond DO_NOT_DOCUMENT
NAMESPACE_BEGIN(Grid);
/// \endcond

/*! \file
 * Classes used for all actions
 *
 * - ConfigurationBase defines a "smart configuration", which permits optional smearing.
 * - Action defines the interface that all actions must follow.
 * - EmptyAction defines a trivial action.
 */

template< class Field >
class ConfigurationBase
/*! @brief Base class for "smart configuration",
 * a container for a gauge configuration holding unsmeared and smeared versions */
{
public:
  ConfigurationBase() {}
  virtual ~ConfigurationBase() {}
  virtual void set_Field(Field& U) =0;
  virtual void smeared_force(Field&) = 0;
  virtual Field& get_SmearedU() =0;
  virtual Field &get_U(bool smeared = false) = 0;
};

template <class GaugeField >
class Action
/*! @brief Base class for all actions. */
{
public:
  bool is_smeared = false;  ///< @brief Is the action computed on a smeared gauge field
  RealD deriv_norm_sum;  ///< @brief Accumulator for absolute value of derivative of action
  RealD deriv_max_sum;  ///< @brief Accumulator for maximum value of derivative of action
  RealD Fdt_norm_sum;  ///< @brief Accumulator for absolute value of impulse (force integrated over time step)
  RealD Fdt_max_sum;  ///< @brief Accumulator for maximum impulse (force integrated over time step)
  int   deriv_num;  ///< @brief Number of times deriv_log() has been called
  RealD deriv_us;  ///< @brief Time spent computing the derivative in microseconds
  RealD S_us;  ///< @brief Time spent computing the action in microseconds
  RealD refresh_us;  ///< @brief Time spent refreshing pseudofermion fields in microseconds

  /*! \name Timing methods */
  /*! @{ */
  /*! @brief %Zero timers associated with this action. */
  void  reset_timer(void)        {
    deriv_us = S_us = refresh_us = 0.0;
    deriv_norm_sum = deriv_max_sum=0.0;
    Fdt_max_sum =  Fdt_norm_sum = 0.0;
    deriv_num=0;
  }
  /*! @brief Record statistics on the action and derivative.
   * @param nrm: Absolute value of the action.
   * @param max: Maximum value of the action.
   * @param Fdt_nrm: Absolute value of the force integrated across this time step.
   * @param Fdt_max: Maximum value of the force integrated across this time step.
   */
  void  deriv_log(RealD nrm, RealD max,RealD Fdt_nrm,RealD Fdt_max) {
    if ( max > deriv_max_sum ) {
      deriv_max_sum=max;
    }
    deriv_norm_sum+=nrm;
    if ( Fdt_max > Fdt_max_sum ) {
      Fdt_max_sum=Fdt_max;
    }
    Fdt_norm_sum+=Fdt_nrm; deriv_num++;
  }
  /*! @brief Mean maximum derivative for a trajectory. */
  RealD deriv_max_average(void)       { return deriv_max_sum; };
  /*! @brief Mean (absolute value of) derivative for a trajectory. */
  RealD deriv_norm_average(void)      { return deriv_norm_sum/deriv_num; };
  /*! @brief Mean maximum value of impulse (force integrated over a time step) for a trajectory. */
  RealD Fdt_max_average(void)         { return Fdt_max_sum; };
  /*! @brief Mean (absolute value of) impulse (force integrated over a time step) for a trajectory.*/
  RealD Fdt_norm_average(void)        { return Fdt_norm_sum/deriv_num; };
  /*! @brief Time spent computing the derivative in microseconds. */
  RealD deriv_timer(void)        { return deriv_us; };
  /*! @brief Time spent computing the action in microseconds. */
  RealD S_timer(void)            { return S_us; };
  /*! @brief Time spent refreshing pseudofermion fields in microseconds. */
  RealD refresh_timer(void)      { return refresh_us; };
  /*! @brief Start timing derivative computation. */
  void deriv_timer_start(void)   { deriv_us-=usecond(); }
  /*! @brief Stop timing derivative computation. */
  void deriv_timer_stop(void)    { deriv_us+=usecond(); }
  /*! @brief Start timing pseudofermion field refresh. */
  void refresh_timer_start(void) { refresh_us-=usecond(); }
  /*! @brief Stop timing pseudofermion field refresh. */
  void refresh_timer_stop(void)  { refresh_us+=usecond(); }
  /*! @brief Start timing action computation. */
  void S_timer_start(void)       { S_us-=usecond(); }
  /*! @brief Stop timing action computation. */
  void S_timer_stop(void)        { S_us+=usecond(); }
  /*! @} */

  /*! \name Methods working on plain gauge fields */
  /*! @{ */
  /////////////////////////////
  // Heatbath?
  /////////////////////////////
  /*! @brief Refresh pseudofermion fields
   * @param U: The gauge field
   * @param sRNG: The serial random number generator to use
   * @param pRNG: The parallel random number generator to use
   */
  virtual void refresh(const GaugeField& U, GridSerialRNG &sRNG, GridParallelRNG& pRNG) = 0;

  /*! @brief Evaluate this action with the given gauge field
   * @param U: The gauge field to evaluate the action of. */
  virtual RealD S(const GaugeField& U) = 0;

  /*! @brief Get the action at the start of the trajectory.
   *
   * If the refresh computes the action,
   * it can be cached so that
   * the value at the start of a Monte Carlo trajectory can be used mid-trajectory.
   * @param U: The gauge field.
   */
  /* An alternative approach would be to have a refreshAndAction() function,
   * with the caller responsible for caching. */
  virtual RealD Sinitial(const GaugeField& U) { return this->S(U); } ;

  /*! Evaluate the derivative of the action,
   * to be integrated in the molecular dynamics part of the HMC.
   * @param U: The gauge field to compute the action's derivative on.
   * @param dSdU: Array into which to output the resulting derivative.
   */
  virtual void deriv(const GaugeField& U, GaugeField& dSdU) = 0;
  /*! @} */

  /*! \name Methods working on smart configuration containers
   *
   * If ``is_smeared``, then these methods work on the smeared field;
   * otherwise, they work on the unsmeared field.
   */
  /*! @{ */
  /*! @brief Refresh pseudofermion fields
   * @param U: The gauge field
   * @param sRNG: The serial random number generator to use
   * @param pRNG: The parallel random number generator to use
   */
  virtual void refresh(ConfigurationBase<GaugeField> & U, GridSerialRNG &sRNG, GridParallelRNG& pRNG)
  {
    refresh(U.get_U(is_smeared),sRNG,pRNG);
  }

  /*! @brief Evaluate this action with the given gauge field
   * @param U: The gauge field to evaluate the action of. */
  virtual RealD S(ConfigurationBase<GaugeField>& U)
  {
    return S(U.get_U(is_smeared));
  }
  /*! @brief Get the action at the start of the trajectory.
   *
   * If the refresh computes the action,
   * it can be cached so that
   * the value at the start of a Monte Carlo trajectory can be used mid-trajectory.
   * @param U: The gauge field.
   */
  /* An alternative approach would be to have a refreshAndAction() function,
   * with the caller responsible for caching. */
  virtual RealD Sinitial(ConfigurationBase<GaugeField>& U) 
  {
    return Sinitial(U.get_U(is_smeared));
  }
  /*! Evaluate the derivative of the action,
   * to be integrated in the molecular dynamics part of the HMC.
   * @param U: The gauge field to compute the action's derivative on.
   * @param dSdU: Array into which to output the resulting derivative.
   */
  virtual void deriv(ConfigurationBase<GaugeField>& U, GaugeField& dSdU)
  {
    deriv(U.get_U(is_smeared),dSdU); 
    if ( is_smeared ) {
      U.smeared_force(dSdU);
    }
  }
  /*! @} */
  /*! \name Logging */
  /*! @{*/
  virtual std::string action_name()    = 0;                             ///< Report the name of the action
  virtual std::string LogParameters()  = 0;                             ///< Print the parameters of the action
  /*! @} */
  virtual ~Action(){}
};

template <class GaugeField >
class EmptyAction : public Action <GaugeField>
/*! @brief A trivial action, which may be used as a placeholder. */
{
  using Action<GaugeField>::refresh;
  using Action<GaugeField>::Sinitial;
  using Action<GaugeField>::deriv;

  virtual void refresh(const GaugeField& U, GridSerialRNG &sRNG, GridParallelRNG& pRNG) { assert(0);}; // refresh pseudofermions
  virtual RealD S(const GaugeField& U) { return 0.0;};                             // evaluate the action
  virtual void deriv(const GaugeField& U, GaugeField& dSdU) { assert(0); };        // evaluate the action derivative

  ///////////////////////////////
  // Logging
  ///////////////////////////////
  virtual std::string action_name()    { return std::string("Level Force Log"); };
  virtual std::string LogParameters()  { return std::string("No parameters");};
};


/// \cond DO_NOT_DOCUMENT
NAMESPACE_END(Grid);
/// \endcond

#endif // ACTION_BASE_H
