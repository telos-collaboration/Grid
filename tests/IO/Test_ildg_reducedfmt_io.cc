    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_nersc_io.cc

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>
Author: Gaurav Ray <gaurav.sinharay@swansea.ac.uk>

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

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
#include <Grid/Grid.h>

using namespace std;
using namespace Grid;

//////////////////////////////////////////////
//
// this template function generates a
// LatticeGaugeField of the chosen gauge
// group passed as a template argument.
//
//////////////////////////////////////////////
template<class gaugeGroup>
LatticeGaugeField generateHotFieldConfiguration( GridCartesian &Fine, std::vector<int> seed ) {

  GridParallelRNG   pRNGa(&Fine);
  std::cout <<GridLogMessage<< " seeding... "<<std::endl;
  pRNGa.SeedFixedIntegers(seed);
  std::cout <<GridLogMessage<< " ...done "<<std::endl;

  LatticeGaugeField Umu(&Fine);

  if constexpr( std::is_same_v<gaugeGroup, GroupName::SU> == true ) {
    SU<Nc>::HotConfiguration(pRNGa,Umu);
  } else if constexpr ( std::is_same_v<gaugeGroup, GroupName::Sp> == true ) {
    Sp<Nc>::HotConfiguration(pRNGa,Umu);
  } else { static_assert(true, "Grid does not recognise the gauge group"); }

  return Umu; 
}

int main (int argc, char ** argv)
{
#ifdef HAVE_LIME
  Grid_init(&argc,&argv);

  std::cout <<GridLogMessage<< " main "<<std::endl;

  auto simd_layout = GridDefaultSimd(4,vComplex::Nsimd());
  auto mpi_layout  = GridDefaultMpi();
  Coordinate latt_size  ({8,8,8,16});

  std::vector<int> seed0 = {15,16,23,42};
  std::vector<int> seed1 = {7,10,19,95};
  std::vector<int> seed2 = {55,34,23,13};

  GridCartesian  Fine(latt_size,simd_layout,mpi_layout);

  // set the gauge group
  using grpName = GroupName::Sp; 

  LatticeGaugeField Umu       = generateHotFieldConfiguration<grpName>(Fine, seed0);
  LatticeGaugeField Umu_diff  = generateHotFieldConfiguration<grpName>(Fine, seed1);
  LatticeGaugeField Umu_saved = generateHotFieldConfiguration<grpName>(Fine, seed2);

  // define enums for different field formats
  FloatingPointFormat const fmt64 = FloatingPointFormat::IEEE64BIG; 
  FloatingPointFormat const fmt32 = FloatingPointFormat::IEEE32BIG; 
  MatrixFormat const noGrpRdc     = MatrixFormat::FULL;
  MatrixFormat const GrpRdc       = MatrixFormat::REDUCED;

  FieldMetaData header;

  std::cout <<GridLogMessage<<"**************************************"<<std::endl;
  std::cout <<GridLogMessage<<"** Writing out FULL FAT ILDG cfg  ****"<<std::endl;
  std::cout <<GridLogMessage<<"**************************************"<<std::endl;
  std::string file("./ckpoint_lat_test_FULLFAT");
  IldgWriter _IldgWriter(Fine.IsBoss());
  _IldgWriter.open(file);
  _IldgWriter.writeConfiguration<grpName,noGrpRdc,fmt64>(Umu,4000,std::string("dummy_ildg_LFN"),std::string("dummy_config"));
  _IldgWriter.close();
  Umu_saved = Umu;
  std::cout <<GridLogMessage<<"**************************************"<<std::endl;
  std::cout <<GridLogMessage<<"** Reading back ILDG conf    *********"<<std::endl;
  std::cout <<GridLogMessage<<"**************************************"<<std::endl;
  IldgReader _IldgReader;
  _IldgReader.open(file);
  _IldgReader.readConfiguration(Umu,header);
  _IldgReader.close();
  Umu_diff = Umu - Umu_saved;
  std::cout <<GridLogMessage<< "norm2 Gauge Diff = "<<norm2(Umu_diff)<<std::endl;

  std::cout <<GridLogMessage<<"**************************************"<<std::endl;
  std::cout <<GridLogMessage<<"** Writing out reduced matrix ILDG cfg **"<<std::endl;
  std::cout <<GridLogMessage<<"**************************************"<<std::endl;
  std::string file1("./ckpoint_lat_test_SEMISKIMMED");
  IldgWriter _IldgWriter1(Fine.IsBoss());
  _IldgWriter1.open(file1);
  _IldgWriter1.writeConfiguration<grpName,GrpRdc,fmt64>(Umu,4000,std::string("dummy_ildg_LFN"),std::string("dummy_config"));
  _IldgWriter1.close();
  Umu_saved = Umu;
  std::cout <<GridLogMessage<<"**************************************"<<std::endl;
  std::cout <<GridLogMessage<<"** Reading back ILDG conf    *********"<<std::endl;
  std::cout <<GridLogMessage<<"**************************************"<<std::endl;
  IldgReader _IldgReader1;
  _IldgReader1.open(file1);
  _IldgReader1.readConfiguration(Umu,header);
  _IldgReader1.close();
  Umu_diff = Umu - Umu_saved;
  std::cout <<GridLogMessage<< "norm2 Gauge Diff = "<<norm2(Umu_diff)<<std::endl;

  std::cout <<GridLogMessage<<"**************************************"<<std::endl;
  std::cout <<GridLogMessage<<"*Writing out reduced single precision ILDG cfg*"<<std::endl;
  std::cout <<GridLogMessage<<"**************************************"<<std::endl;
  std::string file2("./ckpoint_lat_test_SKIMMED");
  IldgWriter _IldgWriter2(Fine.IsBoss());
  _IldgWriter2.open(file2);
  _IldgWriter2.writeConfiguration<grpName,GrpRdc,fmt32>(Umu,4000,std::string("dummy_ildg_LFN"),std::string("dummy_config"));
  _IldgWriter2.close();
  Umu_saved = Umu;
  std::cout <<GridLogMessage<<"**************************************"<<std::endl;
  std::cout <<GridLogMessage<<"** Reading back ILDG conf    *********"<<std::endl;
  std::cout <<GridLogMessage<<"**************************************"<<std::endl;
  IldgReader _IldgReader2;
  _IldgReader2.open(file2);
  _IldgReader2.readConfiguration(Umu,header);
  _IldgReader2.close();
  Umu_diff = Umu - Umu_saved;
  std::cout <<GridLogMessage<< "norm2 Gauge Diff = "<<norm2(Umu_diff)<<std::endl;

  Grid_finalize();
#endif
}
