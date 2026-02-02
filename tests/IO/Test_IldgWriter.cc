    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./Test_IldgWriter.cc

    Copyright (C) 2015

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

using namespace std;
using namespace Grid;

// Unit tests to check the IldgWriter class.

void checkWriteLimeIldgLFN(std::string &test_string);
void checkWriteConfiguration(LatticeGaugeField &Umu);

// write a lime record and then read it back and check it matches.
void checkWriteLimeIldgLFN(std::string &test_string) {
  auto simd_layout = GridDefaultSimd(4,vComplex::Nsimd());
  auto mpi_layout  = GridDefaultMpi();
  Coordinate latt_size  ({8,8,8,16});

  GridCartesian  Grid(latt_size,simd_layout,mpi_layout);

  std::string test_file   = "test_file.lime";
  std::string check_string;
  // open lime file
  IldgWriter _ildgWriter(Grid.IsBoss());
  _ildgWriter.open(test_file);
  _ildgWriter.writeLimeIldgLFN( test_string );
  _ildgWriter.close();
  // read back record 
  IldgReader _ildgReader;
  _ildgReader.open(test_file);
  _ildgReader.readLimeObject(check_string, ILDG_DATA_LFN);
  std::cout <<GridLogMessage<< "test string was " << test_string << std::endl; 
  std::cout <<GridLogMessage<< "check string is " << check_string << std::endl; 
  assert( test_string==check_string );
}

int main (int argc, char ** argv)
{
#ifdef HAVE_LIME
  Grid_init(&argc,&argv);

  std::cout <<GridLogMessage<< " main "<<std::endl;

  std::string test_string = "TestIldgLFN";
  checkWriteLimeIldgLFN(test_string);

  Grid_finalize();
#endif
}
