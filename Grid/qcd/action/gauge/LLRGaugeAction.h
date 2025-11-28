/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/gauge/LLRGaugeAction.h

Copyright (C) 2025

Author: Frederic Bonnet <frederic.bonnet@swansea.ac.uk>

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
#include <Grid/qcd/llr_hmc/llr_hmc.h>

#ifndef QCD_LLR_GAUGE_ACTION_H
#define QCD_LLR_GAUGE_ACTION_H

NAMESPACE_BEGIN(Grid);

template<typename GaugeAction, class Gimpl>
class LLRGaugeAction  : public Action<typename Gimpl::GaugeField> {
public:

    INHERIT_GIMPL_TYPES(Gimpl);

    using Action<GaugeField>::S;
    using Action<GaugeField>::Sinitial;
    using Action<GaugeField>::deriv;
    using Action<GaugeField>::refresh;
    using Action<GaugeField>::set_struct_llrparams;
    using Action<GaugeField>::get_struct_llrparams;

    // instantiating the llr_hmc helper class via the namespace LLR
    namespace_LLR::llrparams* p_llrparams_s;
    namespace_LLR::llr_hmc* p_llr_hmc_o = new namespace_LLR::llr_hmc();
    /////////////////////////////////////////////////////////////
    /// [Constructors]
    /////////////////////////////////////////////////////////////
    explicit LLRGaugeAction(RealD beta_):
                beta(beta_){};

    explicit LLRGaugeAction(namespace_LLR::llrparams* llrparams_):
                p_llrparams_s(llrparams_){};

    explicit LLRGaugeAction(namespace_LLR::llrparams* llrparams_, RealD beta_):
                beta(beta_),
                p_llrparams_s(llrparams_){};

    virtual std::string action_name() {return "LLR_WilsonGaugeAction";}
    /////////////////////////////////////////////////////////////
    /// [Logger]
    /////////////////////////////////////////////////////////////
    virtual std::string LogParameters(){
        std::stringstream sstream;
        sstream << GridLogLLR << "[LLR_WilsonGaugeAction] Beta: " << beta << std::endl;
        return sstream.str();
    }
    /////////////////////////////////////////////////////////////
    /// [Actions]
    /////////////////////////////////////////////////////////////
    virtual void refresh(const GaugeField &U, GridSerialRNG &sRNG, GridParallelRNG &pRNG){};

    virtual RealD S(const GaugeField &U) {
        RealD plaq = WilsonLoops<Gimpl>::avgPlaquette(U);
        RealD vol = U.Grid()->gSites();
        RealD action = (1.0 - plaq) * (Nd * (Nd - 1.0)) * vol * 0.5;
        std::cout << GridLogLLR << B_GREEN << "plaq(nb) LLRGaugeAction         ----->: "<< plaq   << C_RESET <<std::endl;
        std::cout << GridLogLLR << B_GREEN << "Vol(nb)  LLRGaugeAction         ----->: "<< vol    << C_RESET <<std::endl;
        std::cout << GridLogLLR << B_GREEN << "S(nb)    LLRGaugeAction         ----->: "<< action << C_RESET <<std::endl;
        return action;
    };

    virtual RealD Sa(const GaugeField &U) {
        RealD plaq = WilsonLoops<Gimpl>::avgPlaquette(U);
        RealD vol = U.Grid()->gSites();
        RealD action = p_llrparams_s->a * (1.0 - plaq) * (Nd * (Nd - 1.0)) * vol * 0.5;
        std::cout << GridLogLLR << B_MAGENTA << "plaq(a) LLRGaugeAction          ----->: "<< plaq   << C_RESET <<std::endl;
        std::cout << GridLogLLR << B_MAGENTA << "Vol(a)  LLRGaugeAction          ----->: "<< vol    << C_RESET <<std::endl;
        std::cout << GridLogLLR << B_MAGENTA << "S(a)    LLRGaugeAction          ----->: "<< action << C_RESET <<std::endl;
        return action;
    };

    virtual RealD Sbeta(const GaugeField &U) {
        RealD plaq = WilsonLoops<Gimpl>::avgPlaquette(U);
        RealD vol = U.Grid()->gSites();
        RealD action = beta * (1.0 - plaq) * (Nd * (Nd - 1.0)) * vol * 0.5;
        std::cout << GridLogLLR << B_CYAN << "plaq(beta) LLRGaugeAction       ----->: "<< plaq   << C_RESET <<std::endl;
        std::cout << GridLogLLR << B_CYAN << "Vol(beta)  LLRGaugeAction       ----->: "<< vol    << C_RESET <<std::endl;
        std::cout << GridLogLLR << B_CYAN << "S(beta)    LLRGaugeAction       ----->: "<< action << C_RESET <<std::endl;
        return action;
    };
    /////////////////////////////////////////////////////////////
    /// [Deriv]
    /////////////////////////////////////////////////////////////
    virtual void deriv(const GaugeField &U, GaugeField &dSdU) {
        // extend Ta to include Lorentz indexes
        std::cout << GridLogLLR << B_GREEN   << "(deriv)     LLR_WilsonGaugeAction"       << C_RESET   << std::endl;
        std::cout << GridLogLLR << B_MAGENTA << "Nc(deriv) LLRGaugeAction        ----->: "<< RealD(Nc) << C_RESET <<std::endl;

        // Removing the beta dependence on 0.5 * beta / RealD(Nc);
        RealD factor = 0.5  / RealD(Nc);
        GridBase *grid = U.Grid();

        std::cout << GridLogLLR << B_YELLOW  << "factor(deriv) LLRGaugeAction    ----->: "<< factor << C_RESET <<std::endl;

        GaugeLinkField dSdU_mu(grid);
        std::vector<GaugeLinkField> Umu(Nd, grid);

        for (int mu = 0; mu < Nd; mu++) {
            Umu[mu] = PeekIndex<LorentzIndex>(U, mu);
        }

        // Get the average Plaquette
        RealD plaq = WilsonLoops<Gimpl>::avgPlaquette(U);
        RealD vol = U.Grid()->gSites();

        std::cout << GridLogLLR << B_RED     << "Vol(deriv)   LLRGaugeAction     ----->: "<< vol               << C_RESET <<std::endl;
        std::cout << GridLogLLR << B_CYAN    << "plaq(deriv)  LLRGaugeAction     ----->: "<< plaq              << C_RESET <<std::endl;
        std::cout << GridLogLLR << B_MAGENTA << "p->S0(deriv) LLRGaugeAction     ----->: "<< p_llrparams_s->S0 << C_RESET <<std::endl;
        std::cout << GridLogLLR << B_YELLOW  << "p->a(deriv)  LLRGaugeAction     ----->: "<< p_llrparams_s->a  << C_RESET <<std::endl;
        std::cout << GridLogLLR << B_YELLOW  << "p->dS(deriv) LLRGaugeAction     ----->: "<< p_llrparams_s->dS << C_RESET <<std::endl;

        //GaugeLinkField new_dSdU_mu(grid);
        RealD action_u = 0.0;
        action_u = S(U);
        for (int mu = 0; mu < Nd; mu++) {
            // Staple in direction mu
            WilsonLoops<Gimpl>::Staple(dSdU_mu, Umu, mu);

            // Implementation of Eq.(3.3) https://doi.org/10.22323/1.256.0276
            dSdU_mu *= ( p_llrparams_s->a + (action_u - p_llrparams_s->S0) / (p_llrparams_s->dS * p_llrparams_s->dS) );
            dSdU_mu = Ta(Umu[mu] * dSdU_mu) * factor;

            PokeIndex<LorentzIndex>(dSdU, dSdU_mu, mu);
        }
    }
    /////////////////////////////////////////////////////////////
    /// [Setters]
    /////////////////////////////////////////////////////////////
    virtual void set_struct_llrparams(namespace_LLR::llrparams* llrparams_s_in) {
        p_llrparams_s = llrparams_s_in;
    }
    /////////////////////////////////////////////////////////////
    /// [Getters]
    /////////////////////////////////////////////////////////////
    virtual namespace_LLR::llrparams* get_struct_llrparams() {return p_llrparams_s;}
    /////////////////////////////////////////////////////////////
    /// [Helpers]
    /////////////////////////////////////////////////////////////
private:
    RealD beta;
};


NAMESPACE_END(Grid);


#endif //QCD_LLR_GAUGE_ACTION_H
