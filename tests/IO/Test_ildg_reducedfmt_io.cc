    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_nersc_io.cc

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>

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
 ;


int main (int argc, char ** argv)
{
#ifdef HAVE_LIME
  Grid_init(&argc,&argv);

  std::cout <<GridLogMessage<< " main "<<std::endl;

  auto simd_layout = GridDefaultSimd(4,vComplex::Nsimd());
  auto mpi_layout  = GridDefaultMpi();
  //std::vector<int> latt_size  ({48,48,48,96});
  //std::vector<int> latt_size  ({32,32,32,32});
  Coordinate latt_size  ({8,8,8,16});
  Coordinate clatt_size  ({2,2,2,4});
  int orthodir=3;
  int orthosz =latt_size[orthodir];
    
  GridCartesian     Fine(latt_size,simd_layout,mpi_layout);
  GridCartesian     Coarse(clatt_size,simd_layout,mpi_layout);

  GridParallelRNG   pRNGa(&Fine);
  GridParallelRNG   pRNGb(&Fine);
  GridSerialRNG     sRNGa;
  GridSerialRNG     sRNGb;

  std::cout <<GridLogMessage<< " seeding... "<<std::endl;
  pRNGa.SeedFixedIntegers(std::vector<int>({45,12,81,9}));
  sRNGa.SeedFixedIntegers(std::vector<int>({45,12,81,9}));
  std::cout <<GridLogMessage<< " ...done "<<std::endl;

  LatticeGaugeField Umu(&Fine);
  LatticeGaugeField Umu_diff(&Fine);
  LatticeGaugeField Umu_saved(&Fine);

  //std::vector<LatticeColourMatrix> U(4,&Fine);
  
  SU<Nc>::HotConfiguration(pRNGa,Umu);
  using grpName = GroupName::SU; 

  //Sp<Nc>::HotConfiguration(pRNGa,Umu);
  //using grpName = GroupName::Sp; // ideally want to automate this

  FP_FMT const fmt64 = FP_FMT::IEEE64BIG; 
  FP_FMT const fmt32 = FP_FMT::IEEE32BIG; 
  MATRIX_FMT const noGrpRdc = MATRIX_FMT::FULL;
  MATRIX_FMT const GrpRdc =   MATRIX_FMT::REDUCED;

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
