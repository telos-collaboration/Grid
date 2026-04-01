    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./Test_ildg_reducedfmt_io.cc

    Copyright (C) 2026

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
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

using namespace Grid;

// this test demonstrates and checks IldgWriter/Readers' ability to 
// read/write ildg 1.2 compliant SU and Sp lattices, including in reduced format.

//////////////////////////////////////////////
// this template function returns a
// LatticeGaugeField of the chosen gauge
// group (passed as a template argument)
//////////////////////////////////////////////
template<class gaugeGroup>
LatticeGaugeField generateHotFieldConfiguration( GridCartesian &Grid, std::vector<int> seed ) {

  GridParallelRNG   pRNGa(&Grid);
  std::cout <<GridLogMessage<< " seeding... "<<std::endl;
  pRNGa.SeedFixedIntegers(seed);
  std::cout <<GridLogMessage<< " ...done "<<std::endl;

  LatticeGaugeField Umu(&Grid);

  if constexpr( std::is_same_v<gaugeGroup, GroupName::SU> ) {
    SU<Nc>::HotConfiguration(pRNGa,Umu);
  } else if constexpr ( std::is_same_v<gaugeGroup, GroupName::Sp> ) {
    Sp<Nc>::HotConfiguration(pRNGa,Umu);
  } else { static_assert(true, "Grid does not recognise the gauge group"); }

  return Umu; 
}

/////////////////////////////////////////////////////////////
// this template function writes a lattice
// field of a given gaugeGroup to disk. It can write in
// reduced format and single/double precision depending on 
// the values of matrix_fmt and fp_fmt. unique_su toggles
// between different functions when reading reduced fmt lattices. 
/////////////////////////////////////////////////////////////
template<class gaugeGroup, int N, MatrixFormat matrix_fmt, FloatingPointFormat fp_fmt, bool unique_su = false>
void writeReadIldgConfiguration( LatticeGaugeField &Umu, GridCartesian &Grid, std::string file)  {

  if constexpr( std::is_same_v<gaugeGroup,GroupName::Sp> && N%2==1) {
    std::cout <<GridLogMessage<<"**************************************"<<std::endl;
    std::cout <<GridLogMessage<< "CAN NOT WRITE LATTICE" << std::endl;
    std::cout <<GridLogMessage<< "For Sp fields Nc must be even and >= 4" << std::endl;
    std::cout <<GridLogMessage<<"**************************************"<<std::endl;
    return;
  } 

  using stats = PeriodicGaugeStatistics;

  std::cout <<GridLogMessage<<"**************************************"<<std::endl;
  std::cout <<GridLogMessage<<"*** Writing out ILDG cfg ***"<<std::endl;
  std::cout <<GridLogMessage<<"**************************************"<<std::endl;
  IldgWriter _IldgWriter(Grid.IsBoss());
  _IldgWriter.open(file);
  _IldgWriter.writeConfiguration<stats, gaugeGroup, matrix_fmt, fp_fmt>(Umu,4000,std::string("dummy_ildg_LFN"),std::string("dummy_config"));
  _IldgWriter.close();

  LatticeGaugeField Umu_saved(&Grid);

  FieldMetaData header;

  std::cout <<GridLogMessage<<"**************************************"<<std::endl;
  std::cout <<GridLogMessage<<"** Reading back ILDG conf    *********"<<std::endl;
  std::cout <<GridLogMessage<<"**************************************"<<std::endl;
  IldgReader _IldgReader;
  _IldgReader.open(file);
  _IldgReader.readConfiguration<stats, unique_su>(Umu_saved, header);
  _IldgReader.close();
  std::cout <<GridLogMessage<< "norm2 Gauge Diff = "<<norm2((Umu_saved-Umu))<<std::endl;
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
  std::vector<int> seed1 = {95,63,71,66};

  GridCartesian  Grid(latt_size,simd_layout,mpi_layout);

  const bool not_unique_su = false; // use reconstructSU
  const bool unique_su = true;      // use unique_reconstructSU

  LatticeGaugeField Umu   = generateHotFieldConfiguration<GroupName::SU>(Grid, seed0);
  LatticeGaugeField UmuSp = generateHotFieldConfiguration<GroupName::Sp>(Grid, seed1);

  // write and read back SU lattices
  writeReadIldgConfiguration<GroupName::SU, Nc, MatrixFormat::FULL, FloatingPointFormat::IEEE64BIG>(Umu, Grid, "./ckpoint_su64_"+std::to_string(Nc)+"x"+std::to_string(Nc)+".4000");
  writeReadIldgConfiguration<GroupName::SU, Nc, MatrixFormat::REDUCED, FloatingPointFormat::IEEE64BIG, not_unique_su>(Umu, Grid, "./ckpoint_su64_"+std::to_string(Nc-1)+"x"+std::to_string(Nc)+".4000");
  writeReadIldgConfiguration<GroupName::SU, Nc, MatrixFormat::REDUCED, FloatingPointFormat::IEEE64BIG, unique_su>(Umu, Grid, "./ckpoint_su64_unique_"+std::to_string(Nc-1)+"x"+std::to_string(Nc)+".4000");

  writeReadIldgConfiguration<GroupName::SU, Nc, MatrixFormat::FULL, FloatingPointFormat::IEEE32BIG>(Umu, Grid, "./ckpoint_su32_"+std::to_string(Nc)+"x"+std::to_string(Nc)+".4000");
  writeReadIldgConfiguration<GroupName::SU, Nc, MatrixFormat::REDUCED, FloatingPointFormat::IEEE32BIG, not_unique_su>(Umu, Grid, "./ckpoint_su32_"+std::to_string(Nc-1)+"x"+std::to_string(Nc)+".4000");
  writeReadIldgConfiguration<GroupName::SU, Nc, MatrixFormat::REDUCED, FloatingPointFormat::IEEE32BIG, unique_su>(Umu, Grid, "./ckpoint_su32_unique_"+std::to_string(Nc-1)+"x"+std::to_string(Nc)+".4000");
  
  // write and read Sp lattices
  writeReadIldgConfiguration<GroupName::Sp, Nc, MatrixFormat::FULL, FloatingPointFormat::IEEE64BIG, not_unique_su>(UmuSp, Grid, "./ckpoint_sp64_"+std::to_string(Nc)+"x"+std::to_string(Nc)+".4000");
  writeReadIldgConfiguration<GroupName::Sp, Nc, MatrixFormat::REDUCED, FloatingPointFormat::IEEE64BIG>(UmuSp, Grid, "./ckpoint_sp64_"+std::to_string(Nc/2)+"x"+std::to_string(Nc)+".4000");

  writeReadIldgConfiguration<GroupName::Sp, Nc, MatrixFormat::FULL, FloatingPointFormat::IEEE32BIG>(UmuSp, Grid, "./ckpoint_sp32_"+std::to_string(Nc)+"x"+std::to_string(Nc)+".4000");
  writeReadIldgConfiguration<GroupName::Sp, Nc, MatrixFormat::REDUCED, FloatingPointFormat::IEEE32BIG>(UmuSp, Grid, "./ckpoint_sp32_"+std::to_string(Nc/2)+"x"+std::to_string(Nc)+".4000");

  Grid_finalize();
#endif
}
