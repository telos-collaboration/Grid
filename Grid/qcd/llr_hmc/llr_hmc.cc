/*
 * Created by dc-bonn2 on 8/10/25.
*/
//////////////////////////////////////////////////////////////////////////////
// System imports
//////////////////////////////////////////////////////////////////////////////
#include <iostream>
//////////////////////////////////////////////////////////////////////////////
// Application imports
//////////////////////////////////////////////////////////////////////////////
#include "llr_hmc.h"
//////////////////////////////////////////////////////////////////////////////
/// Definitions
//////////////////////////////////////////////////////////////////////////////
/* Return variables */
#define RC_SUCCESS           0 //!<Overall return code (rc) for SUCCESS
#define RC_FAIL             -1 //!<Overall return code (rc) for FAIL
#define RC_STOP             -3 //!<Overall return code (rc) for STOP
#define RC_WARNING          -4 //!<Overall return code (rc) for WARNING
/* Color definitions for color output*/
#define ANSI_COLOR_BRIGHT_BLACK   "\x1b[0;90m"
#define ANSI_COLOR_BRIGHT_RED     "\x1b[0;91m"
#define ANSI_COLOR_BRIGHT_GREEN   "\x1b[0;92m"
#define ANSI_COLOR_BRIGHT_YELLOW  "\x1b[0;93m"
#define ANSI_COLOR_BRIGHT_BLUE    "\x1b[0;94m"
#define ANSI_COLOR_BRIGHT_MAGENTA "\x1b[0;95m"
#define ANSI_COLOR_BRIGHT_CYAN    "\x1b[0;96m"
#define ANSI_COLOR_BRIGHT_WHITE   "\x1b[0;97m"
#define ANSI_COLOR_RESET          "\x1b[0m"
/*! Shorter notation */
#define B_BLACK   ANSI_COLOR_BRIGHT_BLACK
#define B_RED     ANSI_COLOR_BRIGHT_RED
#define B_GREEN   ANSI_COLOR_BRIGHT_GREEN
#define B_YELLOW  ANSI_COLOR_BRIGHT_YELLOW
#define B_BLUE    ANSI_COLOR_BRIGHT_BLUE
#define B_MAGENTA ANSI_COLOR_BRIGHT_MAGENTA
#define B_CYAN    ANSI_COLOR_BRIGHT_CYAN
#define B_WHITE   ANSI_COLOR_BRIGHT_WHITE
/* Resets all the coloring */
#define C_RESET        ANSI_COLOR_RESET //!<Short hand notation for color reset
//////////////////////////////////////////////////////////////////////////////
// Namespace declaration
//////////////////////////////////////////////////////////////////////////////
namespace namespace_LLR {
    //////////////////////////////////////////////////////////////////////////
    /// Class constructors
    //////////////////////////////////////////////////////////////////////////
    llr_hmc::llr_hmc() {
        int rc = RC_SUCCESS;
        std::cout<<B_BLUE<<
                 "Class namespace_LLR::llr_hmc() has been instantiated, return code: "
                 <<B_GREEN<<rc<<C_RESET<<std::endl;
        // initialize the struct data structure
        rc = _initialize();
    };
    llr_hmc::llr_hmc(namespace_LLR::llrparams *llpr_in) {
        int rc = RC_SUCCESS;
        std::cout<<B_BLUE<<
                 "Class namespace_LLR::llr_hmc() has been instantiated, return code: "
                 <<B_GREEN<<rc<<C_RESET<<std::endl;
        // initialize the struct data structure
        init_robbinsmonro(llpr_in);
    };
    //////////////////////////////////////////////////////////////////////////
    /// Helpers
    //////////////////////////////////////////////////////////////////////////
    //reset it to 0 and rhoa to initial value
    void llr_hmc::restart_robbinsmonro(int startit) {
        s_llrparams->it = startit;
        s_llrparams->a = s_llrparams->starta;
    }
    //////////////////////////////////////////////////////////////////////////
    /// Initialisors
    //////////////////////////////////////////////////////////////////////////
    int llr_hmc::_initialize() {
        int rc = RC_SUCCESS;
        // initialized all values to 0
        s_llrparams = (struct llrparams*)malloc(sizeof(struct llrparams));
        init_llr_params();
        // print the reseted values.
        rc = print_s_llrparams(s_llrparams);
        return rc;
    } /* end of _initialize method */

    void llr_hmc::init_robbinsmonro(namespace_LLR::llrparams *llrp_in) {
        s_llrparams->nrm = llrp_in->nrm;
        s_llrparams->nth = llrp_in->nth;
        s_llrparams->it = llrp_in->it;
        s_llrparams->cfactor = llrp_in->cfactor;
        s_llrparams->umb_RM_freq = llrp_in->umb_RM_freq;
        s_llrparams->umb_meas_freq = llrp_in->umb_meas_freq;
        s_llrparams->umb_therm_freq = llrp_in->umb_therm_freq;
        s_llrparams->starta = llrp_in->starta;
        s_llrparams->a = llrp_in->starta;
        s_llrparams->dS = llrp_in->dS;
        s_llrparams->S0 = llrp_in->S0;
    } /* end of init_robbinsmonro method */

    void llr_hmc::init_llr_params() {
        s_llrparams->nrm = 0;
        s_llrparams->nth = 0;
        s_llrparams->it = 0;
        s_llrparams->cfactor = 0;
        s_llrparams->umb_RM_freq = 0;
        s_llrparams->umb_meas_freq = 0;
        s_llrparams->umb_therm_freq = 0;
        s_llrparams->starta = 0.0;
        s_llrparams->a = 0.0;
        s_llrparams->dS = 0.0;
        s_llrparams->S0 = 0.0;
    } /* end of init_llr_params method */
    //////////////////////////////////////////////////////////////////////////
    /// Printers
    //////////////////////////////////////////////////////////////////////////
    int llr_hmc::print_s_llrparams(namespace_LLR::llrparams *s_llr_parameters) {
        int rc = RC_SUCCESS;
        std::cout<<B_BLUE<<"*--------- llrparams struct -----------*"<<std::endl;
        std::cout<<B_BLUE<<"s_llr_parameters->nrm              ---> "<<B_YELLOW<<s_llr_parameters->nrm<<std::endl;
        std::cout<<B_BLUE<<"s_llr_parameters->nth              ---> "<<B_YELLOW<<s_llr_parameters->nth<<std::endl;
        std::cout<<B_BLUE<<"s_llr_parameters->umb_RM_freq      ---> "<<B_YELLOW<<s_llr_parameters->umb_RM_freq<<std::endl;
        std::cout<<B_BLUE<<"s_llr_parameters->umb_meas_freq    ---> "<<B_YELLOW<<s_llr_parameters->umb_meas_freq<<std::endl;
        std::cout<<B_BLUE<<"s_llr_parameters->umb_therm_freq   ---> "<<B_YELLOW<<s_llr_parameters->umb_therm_freq<<std::endl;
        std::cout<<B_BLUE<<"s_llr_parameters->cfactor          ---> "<<B_YELLOW<<s_llr_parameters->cfactor<<std::endl;
        std::cout<<B_BLUE<<"s_llr_parameters->it               ---> "<<B_YELLOW<<s_llr_parameters->it<<std::endl;
        std::cout<<B_BLUE<<"s_llr_parameters->starta           ---> "<<B_YELLOW<<s_llr_parameters->starta<<std::endl;
        std::cout<<B_BLUE<<"s_llr_parameters->a                ---> "<<B_YELLOW<<s_llr_parameters->a<<std::endl;
        std::cout<<B_BLUE<<"s_llr_parameters->S0               ---> "<<B_YELLOW<<s_llr_parameters->S0<<std::endl;
        std::cout<<B_BLUE<<"s_llr_parameters->dS               ---> "<<B_YELLOW<<s_llr_parameters->dS<<std::endl;
        std::cout<<B_BLUE<<"*--------------------------------------*"<<std::endl;
        std::cout<<C_RESET;
        return rc;
    } /* end of print_IPAddresses_data_structure_t(IPAddresses_struct *s_IPAdds) constructor */
    //////////////////////////////////////////////////////////////////////////
    /// Getters
    //////////////////////////////////////////////////////////////////////////
    double llr_hmc::get_llr_a(void) { return s_llrparams->a; }
    double llr_hmc::getS0(void) { return s_llrparams->S0; }
    double llr_hmc::getdS(void) { return s_llrparams->dS; }

    /*
    I would suggest the following priorities:

    Understand how the code in:

    force_llr_0.c
    force_llr_hmc.c
    mon_llr_hmc.c
    mon_llr_gauge.c

    relate to the equations in the paper that Davide has been showing us.
    Port these into Grid into the actions directory.

    Identify how these actions may connect up with the Metropolis test,
    possibly via adjusting the integrator. Adjust the Metropolis test.

    */

} /* end of namespace_LLR */
