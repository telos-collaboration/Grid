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

// Unit tests to check the various (un)munger classes in parallelIO/Metadata.h
void checkBinarySimpleMungers();
void checkGaugeSimpleMungers();
void checkGaugeDoubleStoredMungers();
void checkGauge3x2mungers();
void checkGaugeSpmungers();

// mungeing between the same type should yield the same result.
void checkBinarySimpleMungers() {
  auto simd_layout = GridDefaultSimd(4,vComplex::Nsimd());
  auto mpi_layout  = GridDefaultMpi();
  Coordinate latt_size = GridDefaultLatt();

  GridCartesian  Grid(latt_size,simd_layout,mpi_layout);

  LorentzColourMatrixF in_scalar_objectF;  // single precision
  LorentzColourMatrixF out_scalar_objectF;  // single precision
  LorentzColourMatrixD in_scalar_objectD;  // double precision
  LorentzColourMatrixD out_scalar_objectD;  // double precision
  
  out_scalar_objectF = Zero();
  out_scalar_objectD = Zero();

  in_scalar_objectF = 1.0;
  in_scalar_objectD = 1.0;
 
  // BinarySimpleUnmunger<out_type,in_type> 
  BinarySimpleUnmunger<LorentzColourMatrixF,LorentzColourMatrixF> unmungerFF;
  BinarySimpleUnmunger<LorentzColourMatrixD,LorentzColourMatrixF> unmungerDF;
  BinarySimpleUnmunger<LorentzColourMatrixD,LorentzColourMatrixD> unmungerDD;
  BinarySimpleUnmunger<LorentzColourMatrixF,LorentzColourMatrixD> unmungerFD;

  unmungerFF(in_scalar_objectF, out_scalar_objectF);
  assert(in_scalar_objectF==out_scalar_objectF);
  //assert(in_scalar_objectD==out_scalar_objectD);

  // BinarySimpleMunger<in_type,out_type> 
  BinarySimpleMunger<LorentzColourMatrixF,LorentzColourMatrixF> mungerFF;
  BinarySimpleMunger<LorentzColourMatrixD,LorentzColourMatrixF> mungerDF;
  BinarySimpleMunger<LorentzColourMatrixD,LorentzColourMatrixD> mungerDD;
  BinarySimpleMunger<LorentzColourMatrixF,LorentzColourMatrixD> mungerFD;
  
}

void checkGaugeSpmungers() {
  ColourMatrix Sp_field;

  const Complex  a(0.5, 0.5), abar(0.5, -0.5);
  const Complex  b(0.3, 0.9), bbar(0.3, -0.9);

  // fill top left
  for(int i=0;i<Nc/2;i++) {
    for(int j=0;j<Nc/2;j++) {
      Sp_field()()(i,j) = a;
    }
  }
  // fill top right
  for(int i=0;i<Nc/2;i++) {
    for(int j=2;j<Nc;j++) {
      Sp_field()()(i,j) = b;
    }
  }
  // fill bottom left
  for(int i=2;i<Nc;i++) {
    for(int j=0;j<Nc/2;j++) {
      Sp_field()()(i,j) = -bbar;
    }
  }
  // fill bottom right
  for(int i=2;i<Nc;i++) {
    for(int j=2;j<Nc;j++) {
      Sp_field()()(i,j) = abar;
    }
  }

  LorentzColourMatrixD test = Zero();

  for(int mu=0; mu<Nd; mu++){
    pokeLorentz(test,Sp_field,mu);
  }

  std::cout << test << std::endl;
  
  // reduce field. GaugeSpunmunger<out_type,in_type>
  GaugeSpunmunger<LorentzColourNx2ND,LorentzColourMatrixD> unmunger;
  LorentzColourNx2ND test_rdc = Zero();
  unmunger(test,test_rdc);

  //std::cout << test_rdc << std::endl;
 
  // reconstruct full field. GaugeSpmunger<in_type,out_type>
  GaugeSpmunger<LorentzColourNx2ND,LorentzColourMatrixD> munger;
  LorentzColourMatrixD test_recon;
  munger(test_rdc,test_recon);

  std::cout << test_recon << std::endl;

  // round-trip test
  assert(test==test_recon);


}


int main (int argc, char ** argv)
{
#ifdef HAVE_LIME
  Grid_init(&argc,&argv);

  std::cout <<GridLogMessage<< " main "<<std::endl;

  checkBinarySimpleMungers();
  //checkGaugeSimpleMungers();
  //checkGaugeDoubleStoredMungers();
  //checkGauge3x2mungers();
  checkGaugeSpmungers();

  Grid_finalize();
#endif
}
