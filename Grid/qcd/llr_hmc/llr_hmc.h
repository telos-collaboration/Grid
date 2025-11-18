/*
 * Created by dc-bonn2 on 8/10/25.
 */

//#include <Grid/Grid.h>

#ifndef GRID_TELOS_DEVEL_LLR_HMC_H
#define GRID_TELOS_DEVEL_LLR_HMC_H

///////////////////////////////////////////////////////////////////////////////////
// Application includes
///////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
/// Definitions
//////////////////////////////////////////////////////////////////////////////
/* Return variables */
#define RC_SUCCESS           0 //!<Overall return code (rc) for SUCCESS
#define RC_FAIL             -1 //!<Overall return code (rc) for FAIL
#define RC_STOP             -3 //!<Overall return code (rc) for STOP
#define RC_WARNING          -4 //!<Overall return code (rc) for WARNING
/* Color definitions for color output*/
#define ANSI_COLOR_RESET          "\x1b[0m"
#define ANSI_COLOR_BRIGHT_BLACK   "\x1b[0;30m"
#define ANSI_COLOR_BRIGHT_RED     "\x1b[0;31m"
#define ANSI_COLOR_BRIGHT_GREEN   "\x1b[0;32m"
#define ANSI_COLOR_BRIGHT_YELLOW  "\x1b[0;33m"
#define ANSI_COLOR_BRIGHT_BLUE    "\x1b[0;34m"
#define ANSI_COLOR_BRIGHT_MAGENTA "\x1b[0;35m"
#define ANSI_COLOR_BRIGHT_CYAN    "\x1b[0;36m"
#define ANSI_COLOR_BRIGHT_WHITE   "\x1b[0;37m"
/* Bold faced coloring */
#define B_ANSI_COLOR_BRIGHT_BLACK   "\x1b[0;90m"
#define B_ANSI_COLOR_BRIGHT_RED     "\x1b[0;91m"
#define B_ANSI_COLOR_BRIGHT_GREEN   "\x1b[0;92m"
#define B_ANSI_COLOR_BRIGHT_YELLOW  "\x1b[0;93m"
#define B_ANSI_COLOR_BRIGHT_BLUE    "\x1b[0;94m"
#define B_ANSI_COLOR_BRIGHT_MAGENTA "\x1b[0;95m"
#define B_ANSI_COLOR_BRIGHT_CYAN    "\x1b[0;96m"
#define B_ANSI_COLOR_BRIGHT_WHITE   "\x1b[0;97m"
/*! Shorter notation */
#define C_BLACK   ANSI_COLOR_BRIGHT_BLACK
#define C_RED     ANSI_COLOR_BRIGHT_RED
#define C_GREEN   ANSI_COLOR_BRIGHT_GREEN
#define C_YELLOW  ANSI_COLOR_BRIGHT_YELLOW
#define C_BLUE    ANSI_COLOR_BRIGHT_BLUE
#define C_MAGENTA ANSI_COLOR_BRIGHT_MAGENTA
#define C_CYAN    ANSI_COLOR_BRIGHT_CYAN
#define C_WHITE   ANSI_COLOR_BRIGHT_WHITE
/* bold */
#define B_BLACK   B_ANSI_COLOR_BRIGHT_BLACK
#define B_RED     B_ANSI_COLOR_BRIGHT_RED
#define B_GREEN   B_ANSI_COLOR_BRIGHT_GREEN
#define B_YELLOW  B_ANSI_COLOR_BRIGHT_YELLOW
#define B_BLUE    B_ANSI_COLOR_BRIGHT_BLUE
#define B_MAGENTA B_ANSI_COLOR_BRIGHT_MAGENTA
#define B_CYAN    B_ANSI_COLOR_BRIGHT_CYAN
#define B_WHITE   B_ANSI_COLOR_BRIGHT_WHITE
/* Resets all the coloring */
#define C_RESET        ANSI_COLOR_RESET //!<Short hand notation for color reset
///////////////////////////////////////////////////////////////////////////////////
// Namespace declaration
///////////////////////////////////////////////////////////////////////////////////
namespace namespace_LLR {
    ///////////////////////////////////////////////////////////////////////////////
    /// Data structure
    ///////////////////////////////////////////////////////////////////////////////
    typedef struct llrparams {
        int nrm, nth;
        int umb_RM_freq, umb_meas_freq, umb_therm_freq;
        int cfactor;
        int it;
        double starta;
        double a;
        double S0;
        double dS;
    } llrparams_t;

    //static llrparams llrp;

    typedef struct reppar {
        double S_llr;
        double dS;
        double S0;
        double a;
        double starta;
        double deltaS;
        int rep;
        int repnext;
    } reppar_t;

/////////////////////////////////////////////////////////////
/// Struct
/////////////////////////////////////////////////////////////
    typedef struct hmc_params_llr {
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
    } hmc_params_llr_t;
    ///////////////////////////////////////////////////////////////////////////////
    /// Class
    ///////////////////////////////////////////////////////////////////////////////
    class llr_hmc {
    private:
    public:
        llrparams* s_llrparams;
        hmc_params_llr* s_hmc_params_llr;
        //////////////////////////////////////////////////////////////////////////
        /// [Constructors]
        //////////////////////////////////////////////////////////////////////////
        llr_hmc();
        llr_hmc(namespace_LLR::llrparams *llpr_in);
        llr_hmc(namespace_LLR::llrparams *llpr_in,
                namespace_LLR::hmc_params_llr *hmcParamsLlr_in);
        //////////////////////////////////////////////////////////////////////////
        /// [Initializors]
        //////////////////////////////////////////////////////////////////////////
        int _initialize();
        void init_llr_params ();
        void init_hmc_params_llr(namespace_LLR::hmc_params_llr *hmcParLLR_in);
        void init_robbins_monro(llrparams *llrp_in);
        //////////////////////////////////////////////////////////////////////////
        /// [Helpers]
        //////////////////////////////////////////////////////////////////////////
        void restart_robbins_monro(int start_it_in);
        //////////////////////////////////////////////////////////////////////////
        /// [Drivers]
        //////////////////////////////////////////////////////////////////////////
        int llr_robbins_monro();
        int llr_therm_robbins_monro();
        //////////////////////////////////////////////////////////////////////////
        /// [Printers]
        //////////////////////////////////////////////////////////////////////////
        int print_s_llrparams(llrparams* s_llr_parameters);
        int print_s_hmc_params_llr(hmc_params_llr* s_hmc_llr_parameters);
        //////////////////////////////////////////////////////////////////////////
        /// [Readers]
        //////////////////////////////////////////////////////////////////////////
        hmc_params_llr* ReadCommandLine(int argc, char **argv,
                                        namespace_LLR::hmc_params_llr* HMCParams);
        //////////////////////////////////////////////////////////////////////////
        /// [Setters]
        //////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////
        /// [Getters]
        //////////////////////////////////////////////////////////////////////////
        double get_llr_a(void);
        double getS0(void);
        double getdS(void);
        //////////////////////////////////////////////////////////////////////////
        /// Destructors
        //////////////////////////////////////////////////////////////////////////
        ~llr_hmc();
    };
} /* end of namespace_LLR. */
//extern namespace_LLR::llrparams* s_llrparams;
#endif
