//
// Created by dc-bonn2 on 9/17/25.
//

#include <Grid/qcd/llr_hmc/llr_hmc.h>

#ifndef GRID_TELOS_DEVEL_LLRGAUGEACTION_H
#define GRID_TELOS_DEVEL_LLRGAUGEACTION_H

NAMESPACE_BEGIN(Grid);

//template<class Gimpl>
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

    explicit LLRGaugeAction(namespace_LLR::llrparams* llrparams_, RealD beta_):
                beta(beta_),
                p_llrparams_s(llrparams_){};

    virtual std::string action_name() {return "LLRGaugeAction";}
    /////////////////////////////////////////////////////////////
    /// [Logger]
    /////////////////////////////////////////////////////////////
    virtual std::string LogParameters(){
        std::stringstream sstream;
        sstream << GridLogMessage << "[LLRGaugeAction] Beta: " << beta << std::endl;
        return sstream.str();
    }
    /////////////////////////////////////////////////////////////
    /// [Actions]
    /////////////////////////////////////////////////////////////
    virtual RealD GActionS(const GaugeAction &p_GaugeAction, const GaugeField &U) {
        RealD gauge_action = 0.0;
        //gauge_action =  p_GaugeAction.S(U);
        std::cout <<GridLogIntegrator << B_GREEN << "gauge_action LLRGaugeAction      ----->: "<< gauge_action << C_RESET <<std::endl;
        return gauge_action;
    };

    virtual void refresh(const GaugeField &U, GridSerialRNG &sRNG, GridParallelRNG &pRNG){};

    virtual RealD S(const GaugeField &U) {
        RealD plaq = WilsonLoops<Gimpl>::avgPlaquette(U);
        RealD vol = U.Grid()->gSites();
        //RealD action = beta * (1.0 - plaq) * (Nd * (Nd - 1.0)) * vol * 0.5;
        RealD action = (1.0 - plaq) * (Nd * (Nd - 1.0)) * vol * 0.5;
        std::cout << GridLogIntegrator << B_GREEN << "plaq(nb) LLRGaugeAction         ----->: "<< plaq   << C_RESET <<std::endl;
        std::cout << GridLogIntegrator << B_GREEN << "Vol(nb)  LLRGaugeAction         ----->: "<< vol    << C_RESET <<std::endl;
        std::cout << GridLogIntegrator << B_GREEN << "S(nb)    LLRGaugeAction         ----->: "<< action << C_RESET <<std::endl;
        return action;
    };

    virtual RealD Sa(const GaugeField &U) {
        RealD plaq = WilsonLoops<Gimpl>::avgPlaquette(U);
        RealD vol = U.Grid()->gSites();
        RealD action = p_llrparams_s->a * (1.0 - plaq) * (Nd * (Nd - 1.0)) * vol * 0.5;
        std::cout << GridLogIntegrator << B_MAGENTA << "plaq(a) LLRGaugeAction          ----->: "<< plaq   << C_RESET <<std::endl;
        std::cout << GridLogIntegrator << B_MAGENTA << "Vol(a)  LLRGaugeAction          ----->: "<< vol    << C_RESET <<std::endl;
        std::cout << GridLogIntegrator << B_MAGENTA << "S(a)    LLRGaugeAction          ----->: "<< action << C_RESET <<std::endl;
        return action;
    };

    virtual RealD Sbeta(const GaugeField &U) {
        RealD plaq = WilsonLoops<Gimpl>::avgPlaquette(U);
        RealD vol = U.Grid()->gSites();
        RealD action = beta * (1.0 - plaq) * (Nd * (Nd - 1.0)) * vol * 0.5;
        std::cout << GridLogIntegrator << B_CYAN << "plaq(beta) LLRGaugeAction       ----->: "<< plaq   << C_RESET <<std::endl;
        std::cout << GridLogIntegrator << B_CYAN << "Vol(beta)  LLRGaugeAction       ----->: "<< vol    << C_RESET <<std::endl;
        std::cout << GridLogIntegrator << B_CYAN << "S(beta)    LLRGaugeAction       ----->: "<< action << C_RESET <<std::endl;
        return action;
    };
    /////////////////////////////////////////////////////////////
    /// [Deriv]
    /////////////////////////////////////////////////////////////
    virtual void deriv(const GaugeField &U, GaugeField &dSdU) {
        // not optimal implementation FIXME
        // extend Ta to include Lorentz indexes
        std::cout << GridLogIntegrator << B_GREEN   << "(deriv)     LLRGaugeAction"              << C_RESET <<std::endl;
        std::cout << GridLogIntegrator << B_MAGENTA << "beta(deriv) LLRGaugeAction      ----->: "<< beta             << C_RESET <<std::endl;
        std::cout << GridLogIntegrator << B_CYAN    << "p->a(deriv) LLRGaugeAction b4   ----->: "<< p_llrparams_s->a << C_RESET <<std::endl;
        std::cout << GridLogIntegrator << B_MAGENTA << "Nc(deriv)   LLRGaugeAction      ----->: "<< RealD(Nc)        << C_RESET <<std::endl;

        RealD factor = 0.5  / RealD(Nc);//0.5 * beta / RealD(Nc); //0.5 *
        GridBase *grid = U.Grid();

        std::cout <<GridLogIntegrator << B_YELLOW   << "factor(deriv) LLRGaugeAction    ----->: "<< factor <<"\x1b[0m"<<std::endl;

        GaugeLinkField dSdU_mu(grid);
        std::vector<GaugeLinkField> Umu(Nd, grid);

        for (int mu = 0; mu < Nd; mu++) {
            Umu[mu] = PeekIndex<LorentzIndex>(U, mu);
        }

        // Get the average Plaquette to compute S_0
        RealD plaq = WilsonLoops<Gimpl>::avgPlaquette(U);
        RealD vol = U.Grid()->gSites();
        // now construct the S_0
        RealD S0_plaq_ = 6 * vol * plaq;
        // hard coded value for now TODO: remove this value when properly setup
        RealD S0_ = p_llrparams_s->S0; //173088.00000;

        std::cout << GridLogLLR        << B_RED     << "Vol(deriv)      LLRGaugeAction  ----->: "<< vol               << C_RESET <<std::endl;
        std::cout << GridLogIntegrator << B_RED     << "Vol(deriv)      LLRGaugeAction  ----->: "<< vol               << C_RESET <<std::endl;
        std::cout << GridLogIntegrator << B_CYAN    << "plaq(deriv)     LLRGaugeAction  ----->: "<< plaq              << C_RESET <<std::endl;
        std::cout << GridLogIntegrator << B_YELLOW  << "S0_plaq_(deriv) LLRGaugeAction  ----->: "<< S0_plaq_          << C_RESET <<std::endl;
        std::cout << GridLogIntegrator << B_MAGENTA << "p->S0(deriv)    LLRGaugeAction  ----->: "<< p_llrparams_s->S0 << C_RESET <<std::endl;

        Real a_ = p_llrparams_s->a;  //7.60982;
        RealD dS_ = p_llrparams_s->dS; //37.16129;

        std::cout << GridLogIntegrator << B_YELLOW  << "p->a(deriv)  LLRGaugeAction aft ----->: "<< p_llrparams_s->a  << C_RESET <<std::endl;
        std::cout << GridLogIntegrator << B_YELLOW  << "p->dS(deriv) LLRGaugeAction     ----->: "<< p_llrparams_s->dS << C_RESET <<std::endl;

        //GaugeLinkField new_dSdU_mu(grid);
        RealD action_a = 0.0;
        RealD action_beta = 0.0;
        RealD action_u = 0.0;
        // TODO: to fix the bug one need to initialize back dSdU_mu as *= makes it a divergent function.
        // TODO:  here we need an additional parameter that when we finish the loop over mu reset_coef ---> 0.0

        action_a = Sa(U);
        action_beta = Sbeta(U);

        action_u = S(U);

        for (int mu = 0; mu < Nd; mu++) {
            // Getting the action for when action
            // Staple in direction mu
            WilsonLoops<Gimpl>::Staple(dSdU_mu, Umu, mu);

            //dSdU_mu *= 1.0 * beta;
            //dSdU_mu *= ( p_llrparams_s->a + (S(U) - p_llrparams_s->S0) / (p_llrparams_s->dS * p_llrparams_s->dS) ) ;

            // TODO: do not forget to put that back
            //dSdU_mu *= p_llrparams_s->a;
            dSdU_mu *= ( p_llrparams_s->a + (action_u - p_llrparams_s->S0) / (p_llrparams_s->dS * p_llrparams_s->dS) );
            // TODO: do not forget to put that back

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


#endif //GRID_TELOS_DEVEL_LLRGAUGEACTION_H
