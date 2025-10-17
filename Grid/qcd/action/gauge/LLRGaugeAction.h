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

    // instantiating the llr_hmc helper class via the namespace LLR
    namespace_LLR::llr_hmc* p_llr_hmc_o = new namespace_LLR::llr_hmc();

    explicit LLRGaugeAction(RealD beta_):
                beta(beta_){};

    explicit LLRGaugeAction(namespace_LLR::llrparams* llrparams_, RealD beta_):
                beta(beta_),
                p_llrparams_s(llrparams_){};

    virtual std::string action_name() {return "LLRGaugeAction";}

    virtual std::string LogParameters(){
        std::stringstream sstream;
        sstream << GridLogMessage << "[LLRGaugeAction] Beta: " << beta << std::endl;
        return sstream.str();
    }

    virtual void refresh(const GaugeField &U, GridSerialRNG &sRNG, GridParallelRNG &pRNG){};  // noop as no pseudoferms

    virtual RealD S(const GaugeField &U) {
        RealD plaq = WilsonLoops<Gimpl>::avgPlaquette(U);
        RealD vol = U.Grid()->gSites();
        RealD action = beta * (1.0 - plaq) * (Nd * (Nd - 1.0)) * vol * 0.5;
        std::cout <<GridLogIntegrator << "\x1b[35m"<<"Vol LLRGaugeAction              ----->: "<< vol <<"\x1b[0m"<<std::endl;
        std::cout <<GridLogIntegrator << "\x1b[33m"<<"Action LLRGaugeAction           ----->: "<< action <<"\x1b[0m"<<std::endl;
        return action;
    };

    virtual RealD GActionS(const GaugeAction &p_GaugeAction, const GaugeField &U) {
        RealD gauge_action = 0.0;

        //gauge_action =  p_GaugeAction.S(U);
        std::cout <<GridLogIntegrator << "\x1b[33m"<<"gauge_action LLRGaugeAction     ----->: "<< gauge_action <<"\x1b[0m"<<std::endl;

        return gauge_action;
    };

    virtual void deriv(const GaugeField &U, GaugeField &dSdU) {
        // not optimal implementation FIXME
        // extend Ta to include Lorentz indexes
        std::cout <<GridLogIntegrator << "\x1b[32m"<<"deriv LLRGaugeAction" <<"\x1b[0m"<<std::endl;
        std::cout <<GridLogIntegrator << "\x1b[35m"<<"deriv LLRGaugeAction       beta ----->: "<< beta <<"\x1b[0m"<<std::endl;
        std::cout <<GridLogIntegrator << "\x1b[35m"<<"deriv LLRGaugeAction p_llr_s->a ----->: "<< p_llrparams_s->a <<"\x1b[0m"<<std::endl;
        std::cout <<GridLogIntegrator << "\x1b[35m"<<"deriv LLRGaugeAction         Nc ----->: "<< RealD(Nc) <<"\x1b[0m"<<std::endl;

        RealD factor = 0.5 * p_llrparams_s->a / RealD(Nc);
        GridBase *grid = U.Grid();

        std::cout <<GridLogIntegrator << "\x1b[32m"<<"factor LLRGaugeAction (deriv)   ----->: "<< factor <<"\x1b[0m"<<std::endl;

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

        std::cout <<GridLogIntegrator << "\x1b[35m"<<"vol LLRGaugeAction (deriv)      ----->: "<< vol <<"\x1b[0m"<<std::endl;
        std::cout <<GridLogIntegrator << "\x1b[32m"<<"plaq LLRGaugeAction (deriv)     ----->: "<< plaq <<"\x1b[0m"<<std::endl;
        std::cout <<GridLogIntegrator << "\x1b[33m"<<"S0_plaq_ LLRGaugeAction (deriv) ----->: "<< S0_plaq_ <<"\x1b[0m"<<std::endl;
        std::cout <<GridLogIntegrator << "\x1b[34m"<<"S0_ LLRGaugeAction (deriv)      ----->: "<< S0_ <<"\x1b[0m"<<std::endl;

        Real a_ = p_llrparams_s->a;  //7.60982;
        RealD dS_ = p_llrparams_s->dS; //37.16129;

        std::cout <<GridLogIntegrator << "\x1b[35m"<<"a_ LLRGaugeAction (deriv)       ----->: "<< a_ <<"\x1b[0m"<<std::endl;
        std::cout <<GridLogIntegrator << "\x1b[36m"<<"dS_ LLRGaugeAction (deriv)      ----->: "<< dS_ <<"\x1b[0m"<<std::endl;


        GaugeLinkField new_dSdU_mu(grid);

        for (int mu = 0; mu < Nd; mu++) {
            // Staple in direction mu
            WilsonLoops<Gimpl>::Staple(dSdU_mu, Umu, mu);

            dSdU_mu = Ta(Umu[mu] * dSdU_mu) * factor;

            dSdU_mu *= ( p_llrparams_s->a + (S(U) - p_llrparams_s->S0) / (p_llrparams_s->dS * p_llrparams_s->dS) ) ;

            PokeIndex<LorentzIndex>(dSdU, dSdU_mu, mu);

        }
    }


private:
    RealD beta;
    namespace_LLR::llrparams* p_llrparams_s;
};


NAMESPACE_END(Grid);


#endif //GRID_TELOS_DEVEL_LLRGAUGEACTION_H
