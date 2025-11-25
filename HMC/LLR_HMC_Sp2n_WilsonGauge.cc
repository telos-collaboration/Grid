/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./HMC/LLR_HMC_SpWilsonGauge.cc

Copyright (C) 2015

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
Author: Frederic Bonnet  <frederic.bonnet@swansea.ac.uk>

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
// System
#include <iostream>
#include <fstream>
// Application includes
#include <Grid/Grid.h>

#include <Grid/qcd/llr_hmc/llr_hmc.h>
//////////////////////////////////////////////////////////////////////////////
// Namespace declaration
//////////////////////////////////////////////////////////////////////////////
namespace namespace_LLR {
//////////////////////////////////////////////////////////////////////////////
//// [LLR-Structure]
//////////////////////////////////////////////////////////////////////////////
////////////////////////
/// [ActionLoggerObsParameters]
////////////////////////
struct ActionLoggerObsParameters: Grid::Serializable {
    /*
        Grid::GRID_SERIALIZABLE_CLASS_MEMBERS(ActionLoggerObsParameters,
                Grid::RealD, beta,
                Grid::RealD, a,
                namespace_LLR::llrparams*, s_llrparams_in );
    */
    Grid::RealD beta_ = 0.0;
	int MDsteps_ = 1;
    Grid::RealD a_ = 0.0;
    std::ofstream* logFile_ = nullptr;
    std::ofstream* csvFile_ = nullptr;
    namespace_LLR::llrparams* s_llrparams_in_ = nullptr;
    ActionLoggerObsParameters(Grid::RealD beta = 2.4,
							  int MDsteps = 1,
                              Grid::RealD a = 1.0,
                              namespace_LLR::llrparams* s_llrparams_in = nullptr,
                              std::ofstream* logFile = nullptr,
                              std::ofstream* csvFile = nullptr)
            : beta_(beta), MDsteps_(MDsteps), a_(a), s_llrparams_in_(s_llrparams_in),
              logFile_(logFile), csvFile_(csvFile){}

    template <class ReaderClass >
    ActionLoggerObsParameters(Grid::Reader<ReaderClass>& Reader){
        read(Reader, "ActionLoggerMeasurement", *this);
    }
};
//////////////////////////////////////////////////////////////////////////////
/// Template classes helpers
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//// [LLR-Logger]
//////////////////////////////////////////////////////////////////////////////
template <class Impl>
class LLRActionLogger : public Grid::HmcObservable<typename Impl::Field> {
    ActionLoggerObsParameters Pars;
public:
    // necessary for HmcObservable compatibility
    typedef typename Impl::Field Field;
    Grid::RealD computed_action = 0.0;
    Grid::RealD computed_plaquette = 0.0;
    //////////////////////////////////////////////////////////////////////////
    /// [Constructor]
    //////////////////////////////////////////////////////////////////////////
    LLRActionLogger(Grid::RealD beta = 0.4, Grid::RealD a = 0.1): Pars(beta, a){}
    LLRActionLogger(ActionLoggerObsParameters Pars_): Pars(Pars_){}
    //////////////////////////////////////////////////////////////////////////
    /// [Printers]
    //////////////////////////////////////////////////////////////////////////
    int print(ActionLoggerObsParameters Pars_in, int traj, Grid::RealD action_in,
              int vol_in){
        int rc = RC_SUCCESS;
        //TODO: may implement this later.
        return rc;
    }
    //////////////////////////////////////////////////////////////////////////
    /// [Setters]
    //////////////////////////////////////////////////////////////////////////
    void set_plaquette(Grid::RealD plaq_in){computed_plaquette = plaq_in;}
    void set_action(Grid::RealD act_in){computed_action = act_in;}
    //////////////////////////////////////////////////////////////////////////
    /// [Getters]
    //////////////////////////////////////////////////////////////////////////
    Grid::RealD get_plaquette(){return computed_plaquette;}
    Grid::RealD get_action(){return computed_action;}
    //////////////////////////////////////////////////////////////////////////
    /// [Helpers]
    //////////////////////////////////////////////////////////////////////////
    virtual void TrajectoryComplete(int traj,
                                    Grid::ConfigurationBase<Field> &SmartConfig,
                                    Grid::GridSerialRNG &sRNG,
                                    Grid::GridParallelRNG &pRNG)
    {
        std::cout << Grid::GridLogMessage << "+++++++++++++++++++"<<std::endl;
        std::cout << Grid::GridLogMessage << "LLR Unsmeared plaquette"<<std::endl;
        // Un-smeared
        TrajectoryComplete(traj,SmartConfig.get_U(false),sRNG,pRNG);
        std::cout << Grid::GridLogMessage << "LLR Smeared plaquette"<<std::endl;
        // Smeared
        TrajectoryComplete(traj,SmartConfig.get_U(true),sRNG,pRNG);
        std::cout << Grid::GridLogMessage << "+++++++++++++++++++"<<std::endl;
    };

    void TrajectoryComplete(int traj,
                            Field &U,
                            Grid::GridSerialRNG &sRNG,
                            Grid::GridParallelRNG &pRNG) {

        Grid::RealD plaq = Grid::WilsonLoops<Impl>::avgPlaquette(U);
        Grid::RealD vol = U.Grid()->gSites();

        Grid::RealD action = (1.0 - plaq) * (Grid::Nd * (Grid::Nd - 1.0)) * vol * 0.5;

        // putting the results in the setter
        set_action(action);
        set_plaquette(plaq);
        // Try this way see if this works.
        Pars.s_llrparams_in_->S = action;
        Pars.s_llrparams_in_->plaq = plaq;

        // Printing to standard output
        int def_prec = std::cout.precision();
        std::cout
                << Grid::GridLogLLR
                << std::setprecision(std::numeric_limits<Grid::Real>::digits10 + 1)
                << "LLR beta: [ " << Pars.beta_<< " ] < -- > "
                << "MDsteps: [ " << Pars.MDsteps_<< " ] < -- > "
                << "LLR Plaquette: [ " << traj << " ] "<< plaq << " ] < -- > "
                << "LLR Action: [ " << action << " ] < -- > "
                << "LLR a: [ " << Pars.a_ << " ] < -- > "
                << "LLR S0: [ " << Pars.s_llrparams_in_->S0 << " ] < -- > "
                << "Volume: [ " << vol << " ] "
                << std::endl;

        // Writing to log file after checking the file being open.
        std::ostream &out_log = (Pars.logFile_ && Pars.logFile_->is_open()) ? *Pars.logFile_ : std::cout;
        out_log
                << Grid::GridLogLLR
                << std::setprecision(std::numeric_limits<Grid::Real>::digits10 + 1)
                << "LLR beta: [ " << Pars.beta_<< " ] < -- > "
				<< "MDsteps: [ " << Pars.MDsteps_<< " ] < -- > "
                << "LLR Plaquette: [ " << traj << " ] " << plaq << " ] < -- > "
                << "LLR Action: [ " << action << " ] < -- > "
                << "LLR a: [ " << Pars.a_ << " ] < -- > "
                << "LLR S0: [ " << Pars.s_llrparams_in_->S0 << " ] < -- > "
                << "Volume: [ " << vol << " ] "
                << std::endl;

        // Writing to csv file.
        std::ostream &out_csv = (Pars.csvFile_ && Pars.csvFile_->is_open()) ? *Pars.csvFile_ : std::cout;
        out_csv
                << Grid::GridLogLLR
                << std::setprecision(std::numeric_limits<Grid::Real>::digits10 + 1)
                << "," << traj
                << "," << Pars.beta_
				<< "," << Pars.MDsteps_
                << "," << plaq
                << "," << action
                << "," << Pars.a_
                << "," << Pars.s_llrparams_in_->S0
                << "," << vol
                << std::endl;

        std::cout.precision(def_prec);

    }
};
/////////////////////////////////////////////////////////////
//// [LLR-Module]
/////////////////////////////////////////////////////////////
////////////////////////
/// [Plaquette]
////////////////////////
template < class Impl >
class LLRPlaquetteMod: public Grid::ObservableModule<Grid::PlaquetteLogger<Impl>, Grid::NoParameters>{
typedef Grid::ObservableModule<Grid::PlaquetteLogger<Impl>, Grid::NoParameters> ObsBase;
using ObsBase::ObsBase; // for constructors

// acquire resource
virtual void initialize(){
    this->ObservablePtr.reset(new Grid::PlaquetteLogger<Impl>());
}
public:
LLRPlaquetteMod(): ObsBase(Grid::NoParameters()){}
};
////////////////////////
/// [Action]
////////////////////////
template < class Impl >
class LLRActionMod: public Grid::ObservableModule<LLRActionLogger<Impl>, ActionLoggerObsParameters>{
    typedef Grid::ObservableModule< LLRActionLogger<Impl>, ActionLoggerObsParameters > ObsBase;
    using ObsBase::ObsBase; // for constructors
// acquire resource
    virtual void initialize(){
        //this->ObservablePtr.reset(new LLRActionLogger<Impl>(this->Par_.beta_, this->Par_.a_)); //this->Par_
        this->ObservablePtr.reset(new LLRActionLogger<Impl>(this->Par_)); //this->Par_
    }
public:
    LLRActionMod(ActionLoggerObsParameters Par): ObsBase(Par){}
    LLRActionMod(): ObsBase(){}
    Grid::RealD log_action = LLRActionLogger<Impl>(this->Par_).get_action();
    Grid::RealD log_plaquette = LLRActionLogger<Impl>(this->Par_).get_plaquette();
};
} /* end of namespace namespace_LLR */
/////////////////////////////////////////////////////////////
/// Helpers
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/// Main
/////////////////////////////////////////////////////////////
int main(int argc, char **argv) {
    std::cout<<C_RED<<"Start LLR_HMC_Sp2n_WilsonGauge.cc"<<C_RESET<<std::endl;
    // Initializing Grid library environment.
    Grid::Grid_init(&argc, &argv);
    Grid::GridLogLayout();

    // Instantiating the inout parameter structure as a pointer.
    namespace_LLR::llrparams* s_llrparams_in = (struct namespace_LLR::llrparams*)malloc(sizeof(struct namespace_LLR::llrparams));
    // Initialising the structure
    s_llrparams_in->nrm = 700;
    s_llrparams_in->nth = 300;
    s_llrparams_in->it = 1;
    s_llrparams_in->umb_RM_freq = 1;
    s_llrparams_in->umb_meas_freq = 1;
    s_llrparams_in->umb_therm_freq = 1;
    s_llrparams_in->cfactor = 1;
    s_llrparams_in->starta = 5.6; // TODO: ?-- try the wrong value for interval --?
    s_llrparams_in->a = 5.66;
    s_llrparams_in->S0 = 13281.000;
    s_llrparams_in->dS = 3.0;
    // initialising the action and plaquette in the struc.
    s_llrparams_in->S = 1234.5678;
    s_llrparams_in->plaq = 4321.9876;

    // Create and initialise to some arbitrary values the hmc_parameters structure
    namespace_LLR::hmc_params_llr* s_hmc_params_llr_in =
            (struct namespace_LLR::hmc_params_llr*) malloc(sizeof (struct namespace_LLR::hmc_params_llr));
    // Initialising the hmc_params_llr structure
    s_hmc_params_llr_in->saveInterval = 1;
    s_hmc_params_llr_in->StartTrajectory = 1;
    s_hmc_params_llr_in->beta = 2.4;
    s_hmc_params_llr_in->trajL = 1.0;
    s_hmc_params_llr_in->MDsteps = 44;
    s_hmc_params_llr_in->Trajectories = 100;
    s_hmc_params_llr_in->Thermalizations = 20;

    // bringing the llr_hmc object class
    namespace_LLR::llr_hmc* p_llr_hmc_main_o = new namespace_LLR::llr_hmc(s_llrparams_in,
                                                                          s_hmc_params_llr_in);

    // Get the command line reader parameters and update data structure
    s_hmc_params_llr_in = p_llr_hmc_main_o->ReadCommandLine(argc, argv,s_hmc_params_llr_in);

    // Initialize robbins_monro
    p_llr_hmc_main_o->init_robbins_monro(s_llrparams_in);
    // printing the structure to see how it is constructed.
    p_llr_hmc_main_o->print_s_llrparams(s_llrparams_in);
    p_llr_hmc_main_o->print_s_hmc_params_llr(s_hmc_params_llr_in);

    // Constructing the output files using the command line input parameters
    // Creating the output file
    std::cout<<C_RED<<"<---- Creating/Opening output csv/log file ... --->"<<C_RESET<<std::endl;
    int md = s_hmc_params_llr_in->MDsteps;
    std::ofstream run_LLR_HMC_logfile("LLR_HMC_Sp2n_WilsonGauge_MDsteps-"+std::to_string(md)+".log");
    std::ofstream run_LLR_HMC_csvfile("LLR_HMC_Sp2n_WilsonGauge_MDsteps-"+std::to_string(md)+".csv");
    // Writing to csv file the header.
    run_LLR_HMC_csvfile
            << "Grid::GridLogMessage"
            << "," << "traj"
            << "," << "beta"
			<< "," << "MDsteps"
            << "," << "plaq"
            << "," << "action"
            << "," << "a"
            << "," << "S0"
            << "," << "vol"
            << std::endl;
    run_LLR_HMC_logfile<<"<---- Start LLR_HMC_Sp2n_WilsonGauge.cc ---->"<<C_RESET<<std::endl;

    // Start of the main commands.
    bool with_llr = false;
    std::cout<<"<---- with_llr (initialized)     ---->: " << with_llr << std::endl;
    std::cout<<"<---- llr_config                 ---->: " << llr_config << std::endl;
#if defined(llr_config)
    with_llr = true;
#endif
    std::cout<<"<---- with_llr (llr_config)      ---->: " << with_llr << std::endl;
    std::cout<<"<---- Sp2n_config                ---->: " << Sp2n_config << std::endl;

    // if --enable--LLR switch is activated
    if (with_llr) {
        std::cout<<C_CYAN<<"Start of if block ...  with_llr ----->: "<< with_llr << C_RESET <<std::endl;
        // Sp(2n) representation
        typedef Grid::GenericHMCRunnerSpLLR<Grid::MinimumNorm2> HMCWrapperSpLLR;
        HMCWrapperSpLLR TheHMC;

        // Grid from the command line
        TheHMC.ReadCommandLine(argc, argv);
        TheHMC.Resources.AddFourDimGrid("gauge");

        // Getting the check point parameters setup
        Grid::CheckpointerParameters CPparams;
        CPparams.config_prefix = "ckpoint_lat";
        CPparams.rng_prefix = "ckpoint_rng";
        CPparams.saveInterval = s_hmc_params_llr_in->saveInterval;
        CPparams.format = "IEEE64BIG";
        TheHMC.Resources.LoadNerscCheckpointer(CPparams);
        // Getting the random numbers setup
        Grid::RNGModuleParameters RNGpar;
        RNGpar.serial_seeds = "1 2 3 4 5";
        RNGpar.parallel_seeds = "6 7 8 9 10";
        TheHMC.Resources.SetRNGSeeds(RNGpar);

        // Construct observables
        typedef Grid::PlaquetteMod<HMCWrapperSpLLR::ImplPolicy> PlaqObs;
        typedef Grid::TopologicalChargeMod<HMCWrapperSpLLR::ImplPolicy> QObs;
        TheHMC.Resources.AddObservable<PlaqObs>();
        //////////////////////////////////////////////
        Grid::TopologyObsParameters TopParams;
        TopParams.interval = 1;
        TopParams.do_smearing = false;
        TopParams.Smearing.init_step_size = 0.01;
        TopParams.Smearing.tolerance = 1e-5;
        //TopParams.Smearing.steps = 200;
        //TopParams.Smearing.step_size = 0.01;
        TopParams.Smearing.meas_interval = 50;
        TopParams.Smearing.maxTau = 2.0;
        TheHMC.Resources.AddObservable<QObs>(TopParams);
        //////////////////////////////////////////////
        typedef namespace_LLR::LLRPlaquetteMod<HMCWrapperSpLLR::ImplPolicy> LLRSpPlaqObs;
        TheHMC.Resources.AddObservable<LLRSpPlaqObs>();

        namespace_LLR::ActionLoggerObsParameters LLRActParams;
        LLRActParams.a_ = s_llrparams_in->a;
        LLRActParams.beta_ = s_hmc_params_llr_in->beta;
        LLRActParams.MDsteps_ = s_hmc_params_llr_in->MDsteps;
        LLRActParams.s_llrparams_in_ = s_llrparams_in;
        LLRActParams.logFile_ = &run_LLR_HMC_logfile;
        LLRActParams.csvFile_ = &run_LLR_HMC_csvfile;

        typedef namespace_LLR::LLRActionMod<HMCWrapperSpLLR::ImplPolicy> LLRSpActObs;
        TheHMC.Resources.AddObservable<LLRSpActObs>(LLRActParams);

        /////////////////////////////////////////////////////////////
        // Collect actions, here use more encapsulation
        // standard LLRGaugeAction with Sp(2n) gauge action

        typedef Grid::LLRGaugeAction<Grid::SpWilsonGaugeActionR , Grid::PeriodicGimplR> LLRGaugeActionR;

        LLRGaugeActionR LLRaction(s_llrparams_in, s_hmc_params_llr_in->beta);

        Grid::ActionLevel<HMCWrapperSpLLR::Field> Level1(1);
        Level1.push_back(&LLRaction);
        TheHMC.TheAction.push_back(Level1);
        /////////////////////////////////////////////////////////////
        // HMC parameters are serialisable
        TheHMC.Parameters.MD.MDsteps = s_hmc_params_llr_in->MDsteps;      // 40;
        TheHMC.Parameters.MD.trajL   = float(s_hmc_params_llr_in->trajL); // 1.0;

        TheHMC.Run();  // no smearing

        // Creating the test PASS/FAIL
        std::cout << Grid::GridLogLLR << "--------------------------------------------------"<<std::endl;
        std::cout
                << Grid::GridLogLLR
                << "Final action and Plaquette:"
                << std::endl;

        std::cout
                << Grid::GridLogLLR
                <<"Action    --->: MDsteps[" << s_hmc_params_llr_in->MDsteps <<"] --->: "
                << C_MAGENTA << s_llrparams_in->S << C_RESET << std::endl;
        std::cout
                << Grid::GridLogLLR
                <<"Plaquette --->: MDsteps[" << s_hmc_params_llr_in->MDsteps <<"] --->: "
                << C_MAGENTA << s_llrparams_in->plaq << C_RESET << std::endl;

        std::cout << Grid::GridLogLLR << "--------------------------------------------------"<<std::endl;

        // End the if block
        std::cout<<C_CYAN<<"End of if block ...    with_llr ----->: "<< with_llr << C_RESET <<std::endl;
    } /* [end-if] with_llr */

    // Finalising Grid environment.
    Grid::Grid_finalize();
    // End statement
    std::cout<<C_RED<<"<---- End LLR_HMC_Sp2n_WilsonGauge.cc ---->"<<C_RESET<<std::endl;
    run_LLR_HMC_logfile <<"<---- End LLR_HMC_Sp2n_WilsonGauge.cc ---->"<<C_RESET<<std::endl;
    // Closing log file
    run_LLR_HMC_logfile.close();
    run_LLR_HMC_csvfile.close();

} /* end of main LLR_HMC_Sp2n_WilsonGauge.cc */
