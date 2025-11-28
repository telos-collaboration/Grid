/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./HMC/HMC_WilsonGauge.cc

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
#include <Grid/Grid.h>

#include <Grid/qcd/llr_hmc/llr_hmc.h>
//////////////////////////////////////////////////////////////////////////////
// Declaration
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//// [HMC-Structure]
//////////////////////////////////////////////////////////////////////////////
////////////////////////
/// [ActionLoggerObsParameters]
////////////////////////
struct ActionLoggerObsParameters: Grid::Serializable {
    Grid::TopologyObsParameters* QObsParameters_ = nullptr;
    Grid::RealD beta_ = 0.0;
    int MDsteps_ = 1;
    std::ofstream* logFile_ = nullptr;
    std::ofstream* csvFile_ = nullptr;
    ActionLoggerObsParameters(Grid::RealD beta = 6.0,
                              int MDsteps = 1,
                              std::ofstream* logFile = nullptr,
                              std::ofstream* csvFile = nullptr,
                              Grid::TopologyObsParameters* QObsParameters = nullptr)
            : beta_(beta), MDsteps_(MDsteps),
              logFile_(logFile), csvFile_(csvFile),
              QObsParameters_(QObsParameters){}

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
class HMCActionLogger : public Grid::HmcObservable<typename Impl::Field> {
    ActionLoggerObsParameters Pars;
public:
    // necessary for HmcObservable compatibility
    typedef typename Impl::Field Field;
    Grid::RealD computed_action = 0.0;
    Grid::RealD computed_plaquette = 0.0;
    //////////////////////////////////////////////////////////////////////////
    /// [Constructor]
    //////////////////////////////////////////////////////////////////////////
    HMCActionLogger(Grid::RealD beta = 0.4): Pars(beta){}
    HMCActionLogger(ActionLoggerObsParameters Pars_): Pars(Pars_){}
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
        std::cout << Grid::GridLogMessage << "HMC Unsmeared plaquette"<<std::endl;
        // Un-smeared
        TrajectoryComplete(traj,SmartConfig.get_U(false),sRNG,pRNG);
        std::cout << Grid::GridLogMessage << "HMC Smeared plaquette"<<std::endl;
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
		//////////////////////////////////////////////////////////////////////////
		/// [Action]
		//////////////////////////////////////////////////////////////////////////
        Grid::RealD action = Pars.beta_ * (1.0 - plaq) * (Grid::Nd * (Grid::Nd - 1.0)) * vol * 0.5;

        // putting the results in the setter
        set_action(action);
        set_plaquette(plaq);
		//////////////////////////////////////////////////////////////////////////
		/// [TopoQ]
		//////////////////////////////////////////////////////////////////////////
		int def_prec = std::cout.precision();
        // Now compute the topological charge
		Grid::RealD T0 = 0.0;
		Grid::RealD topoQ = 0.0;
		if (traj%Pars.QObsParameters_->interval == 0) {

			Field Usmear = U;

			if (Pars.QObsParameters_->do_smearing) {
				Grid::WilsonFlowAdaptive<Grid::PeriodicGimplR>
				WF(Pars.QObsParameters_->Smearing.init_step_size,
				Pars.QObsParameters_->Smearing.maxTau,
				Pars.QObsParameters_->Smearing.tolerance,
				Pars.QObsParameters_->Smearing.meas_interval);
				WF.smear(Usmear, U);
				T0 = WF.energyDensityPlaquette(Pars.QObsParameters_->Smearing.maxTau, Usmear);
			}
			topoQ =  Grid::WilsonLoops<Impl>::TopologicalCharge(Usmear);
		}
		//////////////////////////////////////////////////////////////////////////
		/// [Printing]
		//////////////////////////////////////////////////////////////////////////
        // Printing to standard output
        std::cout
				<< Grid::GridLogMessage
                << std::setprecision(std::numeric_limits<Grid::Real>::digits10 + 1)
                << "HMC beta: [ " << Pars.beta_<< " ] < -- > "
                << "MDsteps: [ " << Pars.MDsteps_<< " ] < -- > "
                << "HMC Plaquette: [ " << traj << " ] ["<< plaq << " ] < -- > "
                << "HMC Action: [ " << action << " ] < -- > "
				<< "T0 : [ " << T0 << " ] < -- > "
				<< "TopoQ : [ " << topoQ << " ] < -- > "
                << "Volume: [ " << vol << " ] "
                << std::endl;

        // Writing to log file after checking the file being open.
        std::ostream &out_log = (Pars.logFile_ && Pars.logFile_->is_open()) ? *Pars.logFile_ : std::cout;
        out_log
				<< Grid::GridLogMessage
                << std::setprecision(std::numeric_limits<Grid::Real>::digits10 + 1)
                << "HMC beta: [ " << Pars.beta_<< " ] < -- > "
                << "MDsteps: [ " << Pars.MDsteps_<< " ] < -- > "
                << "HMC Plaquette: [ " << traj << " ] [" << plaq << " ] < -- > "
                << "HMC Action: [ " << action << " ] < -- > "
				<< "T0 : [ " << T0 << " ] < -- > "
				<< "TopoQ : [ " << topoQ << " ] < -- > "
                << "Volume: [ " << vol << " ] "
                << std::endl;

        // Writing to csv file.
        std::ostream &out_csv = (Pars.csvFile_ && Pars.csvFile_->is_open()) ? *Pars.csvFile_ : std::cout;
        out_csv
				<< Grid::GridLogMessage
                << std::setprecision(std::numeric_limits<Grid::Real>::digits10 + 1)
                << "," << traj
                << "," << Pars.beta_
                << "," << Pars.MDsteps_
                << "," << plaq
				<< "," << action
				<< "," << T0
				<< "," << topoQ
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
class HMCPlaquetteMod: public Grid::ObservableModule<Grid::PlaquetteLogger<Impl>, Grid::NoParameters>{
typedef Grid::ObservableModule<Grid::PlaquetteLogger<Impl>, Grid::NoParameters> ObsBase;
using ObsBase::ObsBase; // for constructors

// acquire resource
virtual void initialize(){
    this->ObservablePtr.reset(new Grid::PlaquetteLogger<Impl>());
}
public:
HMCPlaquetteMod(): ObsBase(Grid::NoParameters()){}
};
////////////////////////
/// [Action]
////////////////////////
template < class Impl >
class HMCActionMod: public Grid::ObservableModule<HMCActionLogger<Impl>, ActionLoggerObsParameters>{
typedef Grid::ObservableModule< HMCActionLogger<Impl>, ActionLoggerObsParameters > ObsBase;
using ObsBase::ObsBase; // for constructors
// acquire resource
virtual void initialize(){
    //this->ObservablePtr.reset(new LLRActionLogger<Impl>(this->Par_.beta_, this->Par_.a_)); //this->Par_
    this->ObservablePtr.reset(new HMCActionLogger<Impl>(this->Par_)); //this->Par_
}
public:
HMCActionMod(ActionLoggerObsParameters Par): ObsBase(Par){}
HMCActionMod(): ObsBase(){}
Grid::RealD log_action = HMCActionLogger<Impl>(this->Par_).get_action();
Grid::RealD log_plaquette = HMCActionLogger<Impl>(this->Par_).get_plaquette();
};
/////////////////////////////////////////////////////////////
/// Helpers
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/// Main
/////////////////////////////////////////////////////////////
int main(int argc, char **argv) {
    std::cout << Grid::GridLogMessage <<"Start HMC_SUn_WilsonGauge.cc" <<std::endl;
    // Initializing Grid library

    Grid::Grid_init(&argc, &argv);
    Grid::GridLogLayout();

    // Create and initialise to some arbitrary values the hmc_parameters structure
    namespace_LLR::hmc_params_llr* s_hmc_params_in =
            (struct namespace_LLR::hmc_params_llr*) malloc(sizeof (struct namespace_LLR::hmc_params_llr));
    // Initialising the hmc_params_llr structure
    s_hmc_params_in->saveInterval = 1;
    s_hmc_params_in->StartTrajectory = 1;
    s_hmc_params_in->beta = 2.4;
    s_hmc_params_in->trajL = 1.0;
    s_hmc_params_in->MDsteps = 44;
    s_hmc_params_in->Trajectories = 100;
    s_hmc_params_in->Thermalizations = 20;

    // bringing the llr_hmc object class
    namespace_LLR::llr_hmc* p_hmc_main_o = new namespace_LLR::llr_hmc(s_hmc_params_in);

    // Get the command line reader parameters and update data structure
    s_hmc_params_in = p_hmc_main_o->ReadCommandLine(argc, argv,s_hmc_params_in);

    // printing the structure to see how it is constructed.
    p_hmc_main_o->print_s_hmc_params_llr(s_hmc_params_in);

    // Constructing the output files using the command line input parameters
    // Creating the output file
    std::cout << Grid::GridLogMessage <<C_RED<<"<---- Creating/Opening output csv/log file ... --->"<<C_RESET<<std::endl;
    int md = s_hmc_params_in->MDsteps;
    std::ofstream run_HMC_logfile("HMC_SUn_WilsonGauge_MDsteps-"+std::to_string(md)+".log");
    std::ofstream run_HMC_csvfile("HMC_SUn_WilsonGauge_MDsteps-"+std::to_string(md)+".csv");
    // Writing to csv file the header.
    run_HMC_csvfile
            << "Grid::GridLogMessage"
            << "," << "traj"
            << "," << "beta"
            << "," << "MDsteps"
            << "," << "plaq"
            << "," << "action"
			<< "," << "T0"
			<< "," << "topoQ"
            << "," << "vol"
            << std::endl;
    run_HMC_logfile << "<---- Start HMC_SUn_WilsonGauge.cc ---->" << C_RESET << std::endl;

    // SU(N) representation
    typedef Grid::GenericHMCRunner <Grid::MinimumNorm2> HMCWrapper;
    HMCWrapper TheHMC;

    // Grid from the command line
    TheHMC.ReadCommandLine(argc, argv);
    TheHMC.Resources.AddFourDimGrid("gauge");
    // Possible to create the module by hand
    // hardcoding parameters or using a Reader

    // Checking the parameters in the used, we will use the same as the
    // Standard wilson ones.
    Grid::CheckpointerParameters CPparams;
    CPparams.config_prefix = "ckpoint_lat";
    CPparams.rng_prefix = "ckpoint_rng";
    CPparams.saveInterval = 100;
    CPparams.format = "IEEE64BIG";

    TheHMC.Resources.LoadNerscCheckpointer(CPparams);

    Grid::RNGModuleParameters RNGpar;
    RNGpar.serial_seeds = "1 2 3 4 5";
    RNGpar.parallel_seeds = "6 7 8 9 10";
    TheHMC.Resources.SetRNGSeeds(RNGpar);

    // Construct observables
    // here there is too much indirection
    typedef Grid::PlaquetteMod<HMCWrapper::ImplPolicy> PlaqObs;
    typedef Grid::TopologicalChargeMod<HMCWrapper::ImplPolicy> QObs;

    TheHMC.Resources.AddObservable<PlaqObs>();
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
    typedef HMCPlaquetteMod<HMCWrapper::ImplPolicy> HMCPlaqObs;
    TheHMC.Resources.AddObservable<HMCPlaqObs>();

    ActionLoggerObsParameters HMCActParams;
    HMCActParams.beta_ = s_hmc_params_in->beta;
    HMCActParams.MDsteps_ = s_hmc_params_in->MDsteps;
    HMCActParams.logFile_ = &run_HMC_logfile;
    HMCActParams.csvFile_ = &run_HMC_csvfile;
    HMCActParams.QObsParameters_ = &TopParams;

    typedef HMCActionMod<HMCWrapper::ImplPolicy> HMCSUnActObs;
    TheHMC.Resources.AddObservable<HMCSUnActObs>(HMCActParams);

    /////////////////////////////////////////////////////////////
    // Collect actions, here use more encapsulation
    // standard Wilson Gauge action

    typedef Grid::WilsonGaugeAction<Grid::PeriodicGimplR> GaugeActionR;

    GaugeActionR action(s_hmc_params_in->beta);

    Grid::ActionLevel<HMCWrapper::Field> Level1(1);
    Level1.push_back(&action);
    //Level1.push_back(WGMod.getPtr());
    TheHMC.TheAction.push_back(Level1);
    /////////////////////////////////////////////////////////////

    // HMC parameters are serialisable
    TheHMC.Parameters.MD.MDsteps = s_hmc_params_in->MDsteps;
    TheHMC.Parameters.MD.trajL   = s_hmc_params_in->trajL;

    TheHMC.Run();  // no smearing

    // Finalizing the Grid library
    Grid::Grid_finalize();
    // End statement
    std::cout << Grid::GridLogMessage << "<---- End HMC_SUn_WilsonGauge.cc ---->" <<std::endl;
    run_HMC_logfile <<"<---- End HMC_SUn_WilsonGauge.cc ---->"<<C_RESET<<std::endl;
    // Closing log file
    run_HMC_logfile.close();
    run_HMC_csvfile.close();

}