/*
 * Created by dc-bonn2 on 8/10/25.
 */

//#include <Grid/Grid.h>

//#ifndef GRID_TELOS_DEVEL_LLR_HMC_H
//#define GRID_TELOS_DEVEL_LLR_HMC_H

namespace namespace_LLR {
  struct llrparams {
    int nrm, nth;
    int umb_RM_freq, umb_meas_freq, umb_therm_freq;
    int cfactor;
    int it;
    double starta;
    double a;
    double S0;
    double dS;
    double S;
    double plaq;
  };
  
  struct hmc_params_llr {
    int saveInterval;
    int StartTrajectory;
    double beta;
    double trajL;
    int MDsteps;
    int Thermalizations;
    int Trajectories;
    //std::string serial_seeds; // = "1 2 3 4 5";
    //std::string parallel_seeds; // = "6 7 8 9 10";
    //std::string cnfg_dir; // = ".";
  };
  
  class llr_hmc {
  private:
  public:
    llrparams* s_llrparams;
    hmc_params_llr* s_hmc_params_llr;

    /// [Constructors]
    llr_hmc();
    llr_hmc(namespace_LLR::llrparams *llpr_in);
    llr_hmc(namespace_LLR::hmc_params_llr *hmcParamsLlr_in);
    llr_hmc(namespace_LLR::llrparams *llpr_in,
	    namespace_LLR::hmc_params_llr *hmcParamsLlr_in);
    /// [Initializers]
    int _initialize();
    void init_llr_params ();
    void init_hmc_params_llr(namespace_LLR::hmc_params_llr *hmcParLLR_in);
    void init_robbins_monro(llrparams *llrp_in);
    /// [Helpers]
    void restart_robbins_monro(int start_it_in);
    /// [Drivers]
    int llr_robbins_monro();
    int llr_therm_robbins_monro();
    /// [Printers]
    int print_s_llrparams(llrparams* s_llr_parameters);
    int print_s_hmc_params_llr(hmc_params_llr* s_hmc_llr_parameters);
    /// [Readers]
    hmc_params_llr* ReadCommandLine(int argc, char **argv,
				    namespace_LLR::hmc_params_llr* HMCParams);
    /// [Getters]
    double get_llr_a(void);
    double get_llr_S0(void);
    double get_llr_dS(void);
    /// Destructors
    ~llr_hmc();
  };
} /* end of namespace_LLR. */
//#endif
