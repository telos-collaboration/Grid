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

/*! \page test_ildg_io Testing Grid's ILDG I/O functionality
 *
 *  This test demonstrates how to use the IldgReader and IldgWriter classes provided by Grid and
    defined in parallelIO/IldgIO.h. It also tests if they work correctly by generating and writing a lattice field to disk, before reading that field back and checking it matches, within tolerance, to the initially generated field.
 *  \section sec1 Setting up Grid
 *  The `C-LIME` library is required when using the Ildg classes in %Grid. If you don't have it [__install it__](https://usqcd-software.github.io/c-lime) and configure Grid with \code{.sh}../configure --with-lime=<path-to-lime>\endcode Alternatively one can append `-I<path-to-lime>` to `CXXFLAGS`. If this is done successfully then the `HAVE_LIME` symbol should be defined.
 *  \snippet{lineno} this setup
 *  \section sec2 Grid layout
 *  Start with some boilerplate code to setup a basic Grid using the default simd and mpi layouts and a relatively small `16^3` by `32` lattice.
 *  \snippet{lineno} this grid
 *  Now we can initialise random number generators and lattice fields on our Grid object. We need three field objects to make the comparison between our initial and final field objects.
 *  \snippet{lineno} this gen_field 
 *  Then we can populate one of the field objects with random `SU(3)` matrices.
 *  \snippet{lineno} this umu
 *  Now we have generated a random lattice that we can use to test the I/O routines. 
 *  \section sec3 ILDG I/O
 *  We create an instance of the IldgWriter object to write our lattice into a file `ckpoint_ildg.4000`. In this case this file is only about 73MB. The arguments `4000`, "dummy_ildg_LFN", and "dummy_config" are metadata that will be written into `ckpoint_ildg.4000`'s header. Note the use of `Fine.IsBoss()` in the constructor of the writer class. Some I/O methods, such as the methods for writing metadata into headers, cannot run in parallel meaning Grid needs a way of only letting one mpi rank execute these methods.
 *  \snippet{lineno} this ildg_write
 *  We store the lattice for later comparison.  
 *  \snippet{lineno} this umu_save
 *  Let's immediately read `ckpoint_ildg.4000` back into memory. We need three things; a filename, an object to hold the lattice field, and a FieldMetaData object (`header`) to store the parsed metadata of the field. 
 *  \snippet{lineno} this ildg_read
 *  \section sec4 Validating the result
 *  Finally we compare the two field objects by computing the norm squared of their difference. While this difference won't be exactly zero it will be numerically zero when computed on lattices whose underlying type is double. With single precision fields it should not be much larger than 10^-10.
 *  \snippet{lineno} this check_norm
 */
int main (int argc, char ** argv)
{//! [setup]
#ifdef HAVE_LIME
  Grid_init(&argc,&argv);
//! [setup]

  std::cout <<GridLogMessage<< " main "<<std::endl;
//! [grid]
  auto simd_layout = GridDefaultSimd(4,vComplex::Nsimd());
  auto mpi_layout  = GridDefaultMpi();
  Coordinate latt_size  ({16,16,16,32});
    
  GridCartesian     Fine(latt_size,simd_layout,mpi_layout);
//! [grid]
//! [gen_field]
  GridParallelRNG   pRNGa(&Fine);

  std::cout <<GridLogMessage<< " seeding... "<<std::endl;
  pRNGa.SeedFixedIntegers(std::vector<int>({45,12,81,9}));
  std::cout <<GridLogMessage<< " ...done "<<std::endl;

  LatticeGaugeField Umu(&Fine);
  LatticeGaugeField Umu_diff(&Fine);
  LatticeGaugeField Umu_saved(&Fine);
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

//! [umu_save]
  Umu_saved = Umu;
//! [umu_save]
  std::cout <<GridLogMessage<<"**************************************"<<std::endl;
  std::cout <<GridLogMessage<<"** Reading back ILDG conf    *********"<<std::endl;
  std::cout <<GridLogMessage<<"**************************************"<<std::endl;
//! [ildg_read]
  IldgReader _IldgReader;
  _IldgReader.open(file);
  _IldgReader.readConfiguration(Umu,header);
  _IldgReader.close();
//! [ildg_read]
  /*! compute difference between initially generated lattice and lattice read from disk */ 
//! [check_norm]
  Umu_diff = Umu - Umu_saved;

  /*! \todo convert this to an assert condition. */
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
