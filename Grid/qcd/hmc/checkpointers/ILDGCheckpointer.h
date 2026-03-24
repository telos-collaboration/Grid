/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/hmc/ILDGCheckpointer.h

Copyright (C) 2016

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
#ifndef ILDG_CHECKPOINTER
#define ILDG_CHECKPOINTER

#ifdef HAVE_LIME

#include <iostream>
#include <sstream>
#include <string>

NAMESPACE_BEGIN(Grid);

// Only for Gauge fields
template <class Implementation>
class ILDGHmcCheckpointer : public BaseHmcCheckpointer<Implementation> {
private:
  CheckpointerParameters Params;

public:
  INHERIT_GIMPL_TYPES(Implementation);
  typedef GaugeStatistics<Implementation> GaugeStats;

  ILDGHmcCheckpointer(const CheckpointerParameters &Params_) { initialize(Params_); }

  void initialize(const CheckpointerParameters &Params_) {
    Params = Params_;

    // check here that the format is valid
    int ieee32big = (Params.format == std::string("IEEE32BIG"));
    int ieee64big = (Params.format == std::string("IEEE64BIG"));

    if (!(ieee64big || ieee32big)) {
      std::cout << GridLogError << "Unrecognized file format " << Params.format
                << std::endl;
      std::cout << GridLogError
                << "Allowed: IEEE32BIG | IEEE64BIG"
                << std::endl;
      exit(1);
    }

    if ( !((Params.group == std::string("su")) || (Params.group == std::string("sp"))) ) {
      std::cout << GridLogError << "Unrecognized gauge group " << Params.group << std::endl;
      std::cout << GridLogError << "Allowed: su | sp" << std::endl;
      exit(1);
    }
  }

  // choose appropriate template instantiation here since the 
  // non-const checkpointer parameters cannot be used directly as
  // template arguments to IldgWriter
  void chooseIldgWriter( std::string format, std::string group, bool reduced_matrix, 
                         std::string lat_obj, int traj, 
                         ConfigurationBase<GaugeField> &SmartConfig) {
 
      GridBase *grid = SmartConfig.get_U(false).Grid();

      IldgWriter _IldgWriter(grid->IsBoss());
      _IldgWriter.open(lat_obj);

      if(format=="IEEE64BIG") {
        if(group=="su" && reduced_matrix) {
          _IldgWriter.writeConfiguration<GaugeStats, GroupName::SU, MatrixFormat::REDUCED, FloatingPointFormat::IEEE64BIG>(SmartConfig.get_U(false), traj, lat_obj, lat_obj);
        }
        else if (group=="su" && !reduced_matrix) {
          _IldgWriter.writeConfiguration<GaugeStats, GroupName::SU, MatrixFormat::FULL, FloatingPointFormat::IEEE64BIG>(SmartConfig.get_U(false), traj, lat_obj, lat_obj);
        }
        else if (group=="sp" && reduced_matrix) {
          _IldgWriter.writeConfiguration<GaugeStats, GroupName::Sp, MatrixFormat::REDUCED, FloatingPointFormat::IEEE64BIG>(SmartConfig.get_U(false), traj, lat_obj, lat_obj);
        }
        else if (group=="sp" && !reduced_matrix) {
          _IldgWriter.writeConfiguration<GaugeStats, GroupName::Sp, MatrixFormat::FULL, FloatingPointFormat::IEEE64BIG>(SmartConfig.get_U(false), traj, lat_obj, lat_obj);
        }
      }
      else if (format=="IEEE32BIG") {
         if(group=="su" && reduced_matrix) {
          _IldgWriter.writeConfiguration<GaugeStats, GroupName::SU, MatrixFormat::REDUCED, FloatingPointFormat::IEEE32BIG>(SmartConfig.get_U(false), traj, lat_obj, lat_obj);
        }
        else if (group=="su" && !reduced_matrix) {
          _IldgWriter.writeConfiguration<GaugeStats, GroupName::SU, MatrixFormat::FULL, FloatingPointFormat::IEEE32BIG>(SmartConfig.get_U(false), traj, lat_obj, lat_obj);
        }
        else if (group=="sp" && reduced_matrix) {
          _IldgWriter.writeConfiguration<GaugeStats, GroupName::Sp, MatrixFormat::REDUCED, FloatingPointFormat::IEEE32BIG>(SmartConfig.get_U(false), traj, lat_obj, lat_obj);
        }
        else if (group=="sp" && !reduced_matrix) {
          _IldgWriter.writeConfiguration<GaugeStats, GroupName::Sp, MatrixFormat::FULL, FloatingPointFormat::IEEE32BIG>(SmartConfig.get_U(false), traj, lat_obj, lat_obj);
        }
      }

      _IldgWriter.close();
  }
   

  void TrajectoryComplete(int traj,
			  ConfigurationBase<GaugeField> &SmartConfig,
			  GridSerialRNG &sRNG,
                          GridParallelRNG &pRNG) {
    if ((traj % Params.saveInterval) == 0) {
      std::string config, rng, smr;
      this->build_filenames(traj, Params, config, smr, rng);
      uint32_t nersc_csum,scidac_csuma,scidac_csumb;
      BinaryIO::writeRNG(sRNG, pRNG, rng, 0,nersc_csum,scidac_csuma,scidac_csumb);
      std::cout << GridLogMessage << "Written BINARY RNG " << rng
                << " checksum " << std::hex 
		<< nersc_csum<<"/"
		<< scidac_csuma<<"/"
		<< scidac_csumb
		<< std::dec << std::endl;

    chooseIldgWriter( Params.format, Params.group, Params.reduced_matrix, config, traj,
                      SmartConfig );

      std::cout << GridLogMessage << "Written ILDG Configuration on " << config
                << " checksum " << std::hex 
		<< nersc_csum<<"/"
		<< scidac_csuma<<"/"
		<< scidac_csumb
		<< std::dec << std::endl;

      if ( Params.saveSmeared ) { 
        chooseIldgWriter( Params.format, Params.group, Params.reduced_matrix, smr, traj, 
                          SmartConfig );
          
        std::cout << GridLogMessage << "Written ILDG Configuration on " << smr
                    << " checksum " << std::hex 
        << nersc_csum<<"/"
        << scidac_csuma<<"/"
        << scidac_csumb
        << std::dec << std::endl;
      }

    }
  };

  void CheckpointRestore(int traj, GaugeField &U, GridSerialRNG &sRNG,
                         GridParallelRNG &pRNG) {
    std::string config, rng, smr;
    this->build_filenames(traj, Params, config, smr, rng);
    this->check_filename(rng);
    this->check_filename(config);

    uint32_t nersc_csum,scidac_csuma,scidac_csumb;
    BinaryIO::readRNG(sRNG, pRNG, rng, 0,nersc_csum,scidac_csuma,scidac_csumb);

    FieldMetaData header;
    IldgReader _IldgReader;
    _IldgReader.open(config);

    if(Params.unique_su) {
    _IldgReader.readConfiguration<GaugeStats, true>(U,header);// format from the header
    }
    else { _IldgReader.readConfiguration<GaugeStats, false>(U,header); }

    _IldgReader.close();

    std::cout << GridLogMessage << "Read ILDG Configuration from " << config
              << " checksum " << std::hex 
	      << nersc_csum<<"/"
	      << scidac_csuma<<"/"
	      << scidac_csumb
	      << std::dec << std::endl;
  };
};

NAMESPACE_END(Grid);

#endif  // HAVE_LIME
#endif  // ILDG_CHECKPOINTER
