    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_ildg_io.cc

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
/*! \file Test_ildg_io.cc
 *  \brief tests whether IldgWriter can write an ildg-format lattice to disk, checking by reading it back into Grid with IldgReader. 
 */

#include <Grid/Grid.h>

using namespace std;
using namespace Grid;
 ;

/*! \page test_ildg_io Testing the ILDG I/O classes
 *  This test demonstrates how to use the IldgReader and IldgWriter classes provided by Grid and defined in parallelIO/IldgIO.h
 *  \subsection subsec1 Setting up Grid
 *  Note that the LIME library is required when using the Ildg classes in %Grid. If you don't have it <a href="https://usqcd-software.github.io/c-lime/">install it</a> and let %Grid know its location with at configure time with the \code{.sh}../configure --with-lime=<path-to-lime>\endcode flag. 
 *  \snippet this setup
 *  \subsection sec2 Grid layout
 *  Some boilerplate code to setup a basic Grid using the defaults and a smaller lattice.
 *  \snippet this grid
 *  this is some interleaved text
 *  \snippet this gen_field 
 *  setup random number generators for generating SU fields 
 *  \snippet this umu
 *  generate SU field using pRNGa and write into Umu
 *  \subsection sec3 ILDG I/O
 *  create an instance of the IldgWriter object to write lattice into a file %ckpoint_ildg.4000 
 *  \snippet this ildg_write
 *  we store the lattice for later comparison. 
 *  \snippet this umu_save
 *  \snippet this ildg_read
 *  finally we need to compare the two lattices by computing the norm of the difference.
 *  \snippet this umu_diff
 *  \snippet this check_norm
 *  \subsection sec4 Scidac I/O
 *  \snippet this scidac_write
 *  \snippet this scidac_read
 */
int main (int argc, char ** argv)
{//! [setup]
#ifdef HAVE_LIME
  Grid_init(&argc,&argv);

  std::cout <<GridLogMessage<< " main "<<std::endl;
//! [setup]
//! [grid]
  auto simd_layout = GridDefaultSimd(4,vComplex::Nsimd());
  auto mpi_layout  = GridDefaultMpi();
  Coordinate latt_size  ({16,16,16,32});
  Coordinate clatt_size  ({4,4,4,8});
  /*! \todo delete these definitions. copied over from nersc_io test - not used here*/
  int orthodir=3;
  int orthosz =latt_size[orthodir];
    
  GridCartesian     Fine(latt_size,simd_layout,mpi_layout);
  GridCartesian     Coarse(clatt_size,simd_layout,mpi_layout);
//! [grid]
//! [gen_field]
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

  std::vector<LatticeColourMatrix> U(4,&Fine);
//! [gen_field]
//! [umu]
  SU<Nc>::HotConfiguration(pRNGa,Umu);
//! [umu]


  FieldMetaData header;

  std::cout <<GridLogMessage<<"**************************************"<<std::endl;
  std::cout <<GridLogMessage<<"** Writing out  ILDG conf    *********"<<std::endl;
  std::cout <<GridLogMessage<<"**************************************"<<std::endl;
//! [ildg_write]
  std::string file("./ckpoint_ildg.4000");
  IldgWriter _IldgWriter(Fine.IsBoss());
  _IldgWriter.open(file);
  _IldgWriter.writeConfiguration(Umu,4000,std::string("dummy_ildg_LFN"),std::string("dummy_config"));
  _IldgWriter.close();
//! [ildg_write]

//! [umu save]  
  Umu_saved = Umu;
//! [umu save]  
  std::cout <<GridLogMessage<<"**************************************"<<std::endl;
  std::cout <<GridLogMessage<<"** Reading back ILDG conf    *********"<<std::endl;
  std::cout <<GridLogMessage<<"**************************************"<<std::endl;
//! [ildg_read]
  IldgReader _IldgReader;
  _IldgReader.open(file);
  _IldgReader.readConfiguration(Umu,header);
  _IldgReader.close();
//! [ildg_read]
//! [umu diff]
  /*! compute difference between initially generated lattice and lattice read from disk */ 
  Umu_diff = Umu - Umu_saved;
//! [umu diff]

  /*! \todo convert this to an assert condition. */
//! [check_norm]
  std::cout <<GridLogMessage<< "norm2 Gauge Diff = "<<norm2(Umu_diff)<<std::endl;
//! [check_norm]

  std::cout <<GridLogMessage<<"**************************************"<<std::endl;
  std::cout <<GridLogMessage<<"** Writing out Scidac conf    *********"<<std::endl;
  std::cout <<GridLogMessage<<"**************************************"<<std::endl;
  file = std::string("./ckpoint_scidac.4000");
  /*! see the comment for emptyUserRecord */
  emptyUserRecord record;
//! [scidac_write]
  ScidacWriter _ScidacWriter(Fine.IsBoss());
  _ScidacWriter.open(file);
  _ScidacWriter.writeScidacFieldRecord(Umu,record);
  _ScidacWriter.close();
//! [scidac_write]

  Umu_saved = Umu;
  std::cout <<GridLogMessage<<"**************************************"<<std::endl;
  std::cout <<GridLogMessage<<"** Reading back Scidac  conf    *********"<<std::endl;
  std::cout <<GridLogMessage<<"**************************************"<<std::endl;
//! [scidac_read]
  ScidacReader _ScidacReader;
  _ScidacReader.open(file);
  _ScidacReader.readScidacFieldRecord(Umu,record);
  _ScidacReader.close();
//! [scidac_read]
  Umu_diff = Umu - Umu_saved;

  /*! \todo convert this to an assert condition. */
  std::cout <<GridLogMessage<< "norm2 Gauge Diff = "<<norm2(Umu_diff)<<std::endl;

  Grid_finalize();
#endif
}
