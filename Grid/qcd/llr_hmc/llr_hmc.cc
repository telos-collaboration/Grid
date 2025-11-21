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
#include <Grid/Grid.h>

#include "llr_hmc.h"
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
        init_robbins_monro(llpr_in);
    };

    llr_hmc::llr_hmc(namespace_LLR::llrparams *llpr_in,
                     namespace_LLR::hmc_params_llr *hmcParamsLlr_in) {
        int rc = RC_SUCCESS;
        std::cout<<B_BLUE<<
                 "Class namespace_LLR::llr_hmc() has been instantiated, return code: "
                 <<B_GREEN<<rc<<C_RESET<<std::endl;
        // initialize the struct data structure
        init_robbins_monro(llpr_in);
        //init_hmc_params_llr(hmcParamsLlr_in);
        //print_s_hmc_params_llr(hmcParamsLlr_in);
    };

    //////////////////////////////////////////////////////////////////////////
    /// [Helpers]
    //////////////////////////////////////////////////////////////////////////
    //reset it to 0 and rhoa to initial value
    void llr_hmc::restart_robbins_monro(int start_it_in) {
        s_llrparams->it = start_it_in;
        s_llrparams->a = s_llrparams->starta;
    }
    //////////////////////////////////////////////////////////////////////////
    /// [Initializors]
    //////////////////////////////////////////////////////////////////////////
    int llr_hmc::_initialize() {
        int rc = RC_SUCCESS;
        // initialized all values to 0
        s_llrparams = (struct llrparams*)malloc(sizeof(struct llrparams));
        init_llr_params();
        // print the initialized values fro struct s_llrparams.
        rc = print_s_llrparams(s_llrparams);
        // initialized values in struct hmc_params_llr
        s_hmc_params_llr = (struct hmc_params_llr*) malloc(sizeof(struct hmc_params_llr));
        init_hmc_params_llr(s_hmc_params_llr);
        rc = print_s_hmc_params_llr(s_hmc_params_llr);
        if (rc != RC_SUCCESS) {rc = RC_WARNING;}

        return rc;
    } /* end of _initialize method */

    void llr_hmc::init_robbins_monro(namespace_LLR::llrparams *llrp_in) {
        s_llrparams->nrm            = llrp_in->nrm;
        s_llrparams->nth            = llrp_in->nth;
        s_llrparams->it             = llrp_in->it;
        s_llrparams->cfactor        = llrp_in->cfactor;
        s_llrparams->umb_RM_freq    = llrp_in->umb_RM_freq;
        s_llrparams->umb_meas_freq  = llrp_in->umb_meas_freq;
        s_llrparams->umb_therm_freq = llrp_in->umb_therm_freq;
        s_llrparams->starta         = llrp_in->starta;
        s_llrparams->a              = llrp_in->starta;
        s_llrparams->dS             = llrp_in->dS;
        s_llrparams->S0             = llrp_in->S0;
    } /* end of init_robbins_monro method */

    void llr_hmc::init_llr_params() {
        s_llrparams->nrm            = 0;
        s_llrparams->nth            = 0;
        s_llrparams->it             = 0;
        s_llrparams->cfactor        = 0;
        s_llrparams->umb_RM_freq    = 0;
        s_llrparams->umb_meas_freq  = 0;
        s_llrparams->umb_therm_freq = 0;
        s_llrparams->starta         = 0.0;
        s_llrparams->a              = 0.0;
        s_llrparams->dS             = 0.0;
        s_llrparams->S0             = 0.0;
    } /* end of init_llr_params method */

    void namespace_LLR::llr_hmc::init_hmc_params_llr(namespace_LLR::hmc_params_llr *hmcParLLR_in) {
        hmcParLLR_in->saveInterval = 1;
        hmcParLLR_in->StartTrajectory = 1;
        hmcParLLR_in->beta = 2.4;
        hmcParLLR_in->trajL = 1.0;
        hmcParLLR_in->MDsteps = 1;
        hmcParLLR_in->Thermalizations = 21;
        hmcParLLR_in->Trajectories = 101;
        /*
        hmcParLLR_in->serial_seeds = "1 2 3 4 5";
        hmcParLLR_in->parallel_seeds = "6 7 8 9 10";
        hmcParLLR_in->cnfg_dir = ".";
         */
    } /* end of hmc_params_llr method */
    //////////////////////////////////////////////////////////////////////////
    /// [Readers]
    //////////////////////////////////////////////////////////////////////////
    hmc_params_llr* llr_hmc::ReadCommandLine(int argc, char **argv,
                                             namespace_LLR::hmc_params_llr* HMCParams ) {
        std::cout<<B_BLUE<<"*------- ReadCommandLine       --------*"<<std::endl;
        std::cout<<B_BLUE<<"argc                               ---> "<<B_YELLOW<<argc<<std::endl;
        for (int i = 0; i < argc; ++i) {
            std::cout
            <<B_BLUE<<"**argv                             ---> "
            <<C_CYAN<<"argv[" << C_GREEN<< i <<C_CYAN<< "] = " << B_YELLOW << argv[i] <<std::endl;
        }

        if (Grid::GridCmdOptionExists(argv, argv + argc, "--savefreq")) {
            HMCParams->saveInterval =
                    std::stoi(Grid::GridCmdOptionPayload(argv, argv + argc, "--savefreq"));
        } else {
            std::cout << Grid::GridLogError << "--savefreq must be specified" << std::endl;
            exit(1);
        }

        if (Grid::GridCmdOptionExists(argv, argv + argc, "--starttraj")) {
            HMCParams->StartTrajectory =
                    std::stod(Grid::GridCmdOptionPayload(argv, argv + argc, "--starttraj"));
        } else {
            std::cout << Grid::GridLogError << "--starttraj must be specified" << std::endl;
        }

        if (Grid::GridCmdOptionExists(argv, argv + argc, "--beta")) {
            HMCParams->beta =
                    std::stod(Grid::GridCmdOptionPayload(argv, argv + argc, "--beta"));
        } else {
            std::cout << Grid::GridLogError << "--beta must be specified" << std::endl;
            exit(1);
        }

        if (Grid::GridCmdOptionExists(argv, argv + argc, "--tlen")) {
            HMCParams->trajL =
                    std::stod(Grid::GridCmdOptionPayload(argv, argv + argc, "--tlen"));
        } else {
            std::cout << Grid::GridLogError << "--tlen must be specified" << std::endl;
            exit(1);
        }

        if (Grid::GridCmdOptionExists(argv, argv + argc, "--nsteps")) {
            HMCParams->MDsteps =
                    std::stoi(Grid::GridCmdOptionPayload(argv, argv + argc, "--nsteps"));
        } else {
            std::cout << Grid::GridLogError << "--nsteps must be specified" << std::endl;
            exit(1);
        }

        if (Grid::GridCmdOptionExists(argv, argv + argc, "--Trajectories")) {
            HMCParams->Trajectories =
                    std::stoi(Grid::GridCmdOptionPayload(argv, argv + argc, "--Trajectories"));
        } else {
            std::cout << Grid::GridLogError << "--Trajectories must be specified" << std::endl;
            exit(1);
        }

        if (Grid::GridCmdOptionExists(argv, argv + argc, "--Thermalizations")) {
            HMCParams->Thermalizations =
                    std::stoi(Grid::GridCmdOptionPayload(argv, argv + argc, "--Thermalizations"));
        } else {
            std::cout << Grid::GridLogError << "--Thermalizations must be specified" << std::endl;
            exit(1);
        }

        std::cout<<B_BLUE<<"*--------------------------------------*"<<std::endl;
        std::cout<<C_RESET;
        return HMCParams;
    }
    //////////////////////////////////////////////////////////////////////////
    /// [Drivers]
    //////////////////////////////////////////////////////////////////////////
    int llr_hmc::llr_robbins_monro() {
        int rc = RC_SUCCESS;
        // initialized all values to 0
        return rc;
    }

    int llr_hmc::llr_therm_robbins_monro() {
        int rc = RC_SUCCESS;
        // initialized all values to 0
        return rc;
    }
    //////////////////////////////////////////////////////////////////////////
    /// [Printers]
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
    } /* end of print_s_llrparams method */

    int llr_hmc::print_s_hmc_params_llr(namespace_LLR::hmc_params_llr *s_hmc_llr_parameters) {
        int rc = RC_SUCCESS;
        std::cout<<B_BLUE<<"*------- hmc_params_llr struct --------*"<<std::endl;
        std::cout<<B_BLUE<<"s_hmc_params_llr->beta             ---> "<<B_YELLOW<<s_hmc_llr_parameters->beta<<std::endl;
        std::cout<<B_BLUE<<"s_hmc_params_llr->saveInterval     ---> "<<B_YELLOW<<s_hmc_llr_parameters->saveInterval<<std::endl;
        std::cout<<B_BLUE<<"s_hmc_params_llr->StartTrajectory  ---> "<<B_YELLOW<<s_hmc_llr_parameters->StartTrajectory<<std::endl;
        std::cout<<B_BLUE<<"s_hmc_params_llr->trajL            ---> "<<B_YELLOW<<s_hmc_llr_parameters->trajL<<std::endl;
        std::cout<<B_BLUE<<"s_hmc_params_llr->MDsteps          ---> "<<B_YELLOW<<s_hmc_llr_parameters->MDsteps<<std::endl;
        std::cout<<B_BLUE<<"s_hmc_params_llr->Thermalizations  ---> "<<B_YELLOW<<s_hmc_llr_parameters->Thermalizations<<std::endl;
        std::cout<<B_BLUE<<"s_hmc_params_llr->Trajectories     ---> "<<B_YELLOW<<s_hmc_llr_parameters->MDsteps<<std::endl;
        std::cout<<B_BLUE<<"*--------------------------------------*"<<std::endl;
        std::cout<<C_RESET;
        return rc;
    } /* end of print_s_hmc_params_llr method */
    //////////////////////////////////////////////////////////////////////////
    /// [Setters]
    //////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////
    /// [Getters]
    //////////////////////////////////////////////////////////////////////////
    double llr_hmc::get_llr_a(void) { return s_llrparams->a;  }
    double llr_hmc::getS0(void)     { return s_llrparams->S0; }
    double llr_hmc::getdS(void)     { return s_llrparams->dS; }
    //////////////////////////////////////////////////////////////////////////
    /// [Comments]
    //////////////////////////////////////////////////////////////////////////
    /*
        I would suggest the following priorities:

        Understand how the code in:

        force_llr_0.c
        force_llr_hmc.c
        mon_llr_hmc.c
        mon_llr_gauge.c

        understand the updater
        update_llr_ghmc

        relate to the equations in the paper that Davide has been showing us.
        Port these into Grid into the actions directory.

        Identify how these actions may connect up with the Metropolis test,
        possibly via adjusting the integrator. Adjust the Metropolis test.

    */
    //////////////////////////////////////////////////////////////////////////
    /// [End]
    //////////////////////////////////////////////////////////////////////////
} /* end of namespace_LLR */
