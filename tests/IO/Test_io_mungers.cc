    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./Test_io_mungers.cc

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

// tests to check the various (un)munger classes 
// and functions in parallelIO/Metadata.h
bool check_reconstruct3();
bool check_reconstructSU();
bool checkGauge3x2mungers();
bool checkGaugeSUmungers();
bool check_reconstructSp();
bool checkGaugeSpmungers();
bool checkBinarySimpleMungers();
bool checkGaugeSimpleMungers();
bool checkGaugeDoubleStoredMungers();

void mock_SU_field();

//////////////////////////////////////////////
// helper function for generating SU matrices
//////////////////////////////////////////////
void mock_SU_field( ColourMatrixD &cm, std::vector<int> seed )  {
  GridSerialRNG sRNG;   sRNG.SeedFixedIntegers(seed);

  ComplexD ca;
  ColourMatrixD lie, la, ta;
  ComplexD ci(0.0, 1.0);

  lie = Zero(); 
    
  for (int a = 0; a < Nc*Nc-1; a++) {
    random(sRNG, ca);

    ca = (ca + conjugate(ca)) * 0.5;
    ca = ca - 0.5;

    SU<Nc>::generator(a, ta);

    la = timesI(ca * ta);

    lie = lie + la;  // e^{i la ta}
  }

  cm = Exponentiate(lie, 2.0);

}
////////////////////////////////////////////

bool check_reconstruct3() {

  ColourMatrixD SU3_field;
  std::vector<int> rng_seed = {1,2,3,4};
  mock_SU_field(SU3_field, rng_seed);

  //set last row equal to zero
  for(int j=0;j<Nc;j++) {
    SU3_field()()(Nc-1,j) = 0.0;
  }

  LorentzColourMatrixD test = Zero();

  for(int mu=0; mu<Nd; mu++){
    pokeLorentz(test,SU3_field,mu);
  }

  reconstruct3(test);

  // check result is unitary and det==1
  for(int mu=0;mu<Nd;mu++) {
    auto new_cm = peekIndex<LorentzIndex>(test, mu);
    auto det = Determinant(new_cm);
    assert( abs(norm2(det)-1.0) < 1e-15 );
    assert( norm2( new_cm*adj(new_cm)-1.0 ) < 1e-15 );
  }

  return true;

}

bool check_reconstructSU() {
  
  ColourMatrixD SUN_field;
  std::vector<int> rnd_seed = {4,8,15,16};
  mock_SU_field(SUN_field, rnd_seed);

  //set last row equal to zero
  for(int j=0;j<Nc;j++) {
    SUN_field()()(Nc-1,j) = 0.0;
  }

  LorentzColourMatrixD test = Zero();

  for(int mu=0; mu<Nd; mu++){
    pokeLorentz(test,SUN_field,mu);
  }

  reconstructSU(test);

  // check result is unitary and det==1
  for(int mu=0;mu<Nd;mu++) {
    auto new_cm = peekIndex<LorentzIndex>(test, mu);
    auto det = Determinant(new_cm);
    assert( abs( norm2(det)-1.0 ) < 1e-15 );
    assert( norm2( new_cm*adj(new_cm)-1.0 ) < 1e-15 );
  }

  return true;

}

bool checkGauge3x2mungers() {

  ColourMatrixD SU3_field;
  std::vector<int> rnd_seed{20,6,15,8};
  mock_SU_field(SU3_field, rnd_seed);

  LorentzColourMatrixD test=Zero(); 
  for(int mu=0; mu<Nd; mu++){
    pokeLorentz(test,SU3_field,mu);
  }

  // reduce field. Gauge3x2unmunger<out_type,in_type>
  Gauge3x2unmunger<LorentzColour2x3D,LorentzColourMatrixD> unmunger;
  LorentzColour2x3D test_rdc;
  unmunger(test,test_rdc);

  // reconstruct full field. Gauge3x2munger<in_type,out_type>
  Gauge3x2munger<LorentzColour2x3D,LorentzColourMatrixD> munger;
  LorentzColourMatrixD test_recon;
  munger(test_rdc,test_recon);

  // round-trip test
  //std::cout << "checkGauge3x2mungers: " << norm2(test_recon-test) << std::endl;
  assert(norm2(test_recon-test)<1e-10);

  return true;

}

bool checkGaugeSUmungers() {
  ColourMatrixD SUN_field;
  std::vector<int> rnd_seed{23,42,66,98};
  mock_SU_field(SUN_field, rnd_seed);

  LorentzColourMatrixD test=Zero(); 

  for(int mu=0; mu<Nd; mu++){
    pokeLorentz(test,SUN_field,mu);
  }

  // reduce field. GaugeSUunmunger<out_type,in_type>
  GaugeSUunmunger<LorentzColour2x3D,LorentzColourMatrixD> unmunger;
  LorentzColour2x3D test_rdc = Zero();
  unmunger(test,test_rdc);

  // reconstruct full field. GaugeSUmunger<in_type,out_type>
  GaugeSUmunger<LorentzColour2x3D,LorentzColourMatrixD> munger;
  LorentzColourMatrixD test_recon;
  munger(test_rdc,test_recon);
  
  // round-trip test
  //std::cout << "checkGaugeSUmungers: " << norm2(test_recon-test) << std::endl;
  assert(norm2(test_recon-test)<1e-15);

  // check result is unitary and det==1
  for(int mu=0;mu<Nd;mu++) {
    auto new_cm = peekIndex<LorentzIndex>(test_recon, mu);
    auto det = Determinant(new_cm);
    assert( abs(norm2(det)-1.0) < 1e-15 );
    assert( norm2( new_cm*adj(new_cm)-1.0 ) < 1e-15 );
  }

  return true;

}


bool check_reconstructSp() {

  LorentzColourMatrix cm = Zero();
  LorentzColourMatrix cm_zero = Zero();

  // should be zero matrix after reconstruction.
  reconstructSp(cm);
  assert(cm==cm_zero);
  
  ColourMatrix Sp_field = Zero();

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

  LorentzColourMatrixD test = Zero();

  for(int mu=0; mu<Nd; mu++){
    pokeLorentz(test,Sp_field,mu);
  }

  reconstructSp(test);

  // test Sp block structure 
  //  A  B
  // -B* A*
  for(int mu=0; mu<Nd; mu++) {
    auto Sp_cm = peekIndex<LorentzIndex>(test,mu);
    for(int i=0;i<Nc/2;i++) {
      for(int j=0;j<Nc/2;j++) {
        assert( Sp_cm()()(i,j) == conjugate(Sp_cm()()(i+Nc/2,j+Nc/2)) );
        assert( Sp_cm()()(i,j+Nc/2) == -conjugate(Sp_cm()()(i+Nc/2,j)) );
      }
    }
  }

  return true;

}

bool checkGaugeSpmungers() {
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

  // reduce field. GaugeSpunmunger<out_type,in_type>
  GaugeSpunmunger<LorentzColourNx2ND,LorentzColourMatrixD> unmunger;
  LorentzColourNx2ND test_rdc = Zero();
  unmunger(test,test_rdc);

  // reconstruct full field. GaugeSpmunger<in_type,out_type>
  GaugeSpmunger<LorentzColourNx2ND,LorentzColourMatrixD> munger;
  LorentzColourMatrixD test_recon;
  munger(test_rdc,test_recon);

  // round-trip test
  assert(test==test_recon);

  return true;

}


// mungeing between the same type should yield the same result.
bool checkBinarySimpleMungers() {

  LorentzColourMatrixF in_scalar_objectF;  // single precision
  LorentzColourMatrixF out_scalar_objectF = Zero();  
  LorentzColourMatrixD in_scalar_objectD;  // double precision
  LorentzColourMatrixD out_scalar_objectD = Zero();  
  
  in_scalar_objectF = 1.0;
  in_scalar_objectD = 1.0;
 
  // BinarySimpleUnmunger<out_type,in_type> 
  BinarySimpleUnmunger<LorentzColourMatrixF,LorentzColourMatrixF> unmungerFF;
  BinarySimpleUnmunger<LorentzColourMatrixD,LorentzColourMatrixF> unmungerDF;
  BinarySimpleUnmunger<LorentzColourMatrixD,LorentzColourMatrixD> unmungerDD;
  BinarySimpleUnmunger<LorentzColourMatrixF,LorentzColourMatrixD> unmungerFD;

  unmungerFF(in_scalar_objectF, out_scalar_objectF);
  assert(in_scalar_objectF==out_scalar_objectF);

  // BinarySimpleMunger<in_type,out_type> 
  BinarySimpleMunger<LorentzColourMatrixF,LorentzColourMatrixF> mungerFF;
  BinarySimpleMunger<LorentzColourMatrixD,LorentzColourMatrixF> mungerDF;
  BinarySimpleMunger<LorentzColourMatrixD,LorentzColourMatrixD> mungerDD;
  BinarySimpleMunger<LorentzColourMatrixF,LorentzColourMatrixD> mungerFD;

  mungerDD(in_scalar_objectD, out_scalar_objectD);
  assert(in_scalar_objectD==out_scalar_objectD);
 
  return true;

}
// these are used in NerscIO.h
bool checkGaugeSimpleMungers() {

  LorentzColourMatrixF in_scalar_objectF;  // single precision
  LorentzColourMatrixF out_scalar_objectF = Zero(); 
  LorentzColourMatrixD in_scalar_objectD;  // double precision
  LorentzColourMatrixD out_scalar_objectD = Zero(); 
  
  in_scalar_objectF = 1.0;
  in_scalar_objectD = 1.0;

  // GaugeSimpleUnmunger<out_type,in_type> 
  GaugeSimpleUnmunger<LorentzColourMatrixF,LorentzColourMatrixF> unmungerFF;
  GaugeSimpleUnmunger<LorentzColourMatrixD,LorentzColourMatrixF> unmungerDF;
  GaugeSimpleUnmunger<LorentzColourMatrixD,LorentzColourMatrixD> unmungerDD;
  GaugeSimpleUnmunger<LorentzColourMatrixF,LorentzColourMatrixD> unmungerFD;

  unmungerFF(in_scalar_objectF, out_scalar_objectF);
  assert(in_scalar_objectF==out_scalar_objectF);

  // GaugeSimpleMunger<in_type,out_type> 
  GaugeSimpleMunger<LorentzColourMatrixF,LorentzColourMatrixF> mungerFF;
  GaugeSimpleMunger<LorentzColourMatrixD,LorentzColourMatrixF> mungerDF;
  GaugeSimpleMunger<LorentzColourMatrixD,LorentzColourMatrixD> mungerDD;
  GaugeSimpleMunger<LorentzColourMatrixF,LorentzColourMatrixD> mungerFD;
 
  mungerDD(in_scalar_objectD, out_scalar_objectD);
  assert(in_scalar_objectD==out_scalar_objectD);

  return true;

}

bool checkGaugeDoubleStoredMungers() {

  DoubleStoredColourMatrixF in_scalar_objectF;  // single precision
  DoubleStoredColourMatrixF out_scalar_objectF = Zero(); 
  DoubleStoredColourMatrixD in_scalar_objectD;  // double precision
  DoubleStoredColourMatrixD out_scalar_objectD = Zero(); 
  
  in_scalar_objectF = 1.0;
  in_scalar_objectD = 1.0;

  // GaugeDoubleStoredUnmunger<out_type,in_type> 
  GaugeDoubleStoredUnmunger<DoubleStoredColourMatrixF,DoubleStoredColourMatrixF> unmungerFF;
  GaugeDoubleStoredUnmunger<DoubleStoredColourMatrixD,DoubleStoredColourMatrixF> unmungerDF;
  GaugeDoubleStoredUnmunger<DoubleStoredColourMatrixD,DoubleStoredColourMatrixD> unmungerDD;
  GaugeDoubleStoredUnmunger<DoubleStoredColourMatrixF,DoubleStoredColourMatrixD> unmungerFD;

  unmungerFF(in_scalar_objectF, out_scalar_objectF);
  assert(in_scalar_objectF==out_scalar_objectF);

  // GaugeDoubleStoredMunger<in_type,out_type> 
  GaugeDoubleStoredMunger<DoubleStoredColourMatrixF,DoubleStoredColourMatrixF> mungerFF;
  GaugeDoubleStoredMunger<DoubleStoredColourMatrixD,DoubleStoredColourMatrixF> mungerDF;
  GaugeDoubleStoredMunger<DoubleStoredColourMatrixD,DoubleStoredColourMatrixD> mungerDD;
  GaugeDoubleStoredMunger<DoubleStoredColourMatrixF,DoubleStoredColourMatrixD> mungerFD;
 
  mungerDD(in_scalar_objectD, out_scalar_objectD);
  assert(in_scalar_objectD==out_scalar_objectD);

  return true;

}


int main (int argc, char ** argv)
{
#ifdef HAVE_LIME
  Grid_init(&argc,&argv);
  
  std::cout <<GridLogMessage<< " main "<<std::endl;

  std::cout <<GridLogMessage<< "Testing SU(" << Nc << ") mungers..." << std::endl;
  if constexpr(Nc>1 && Nc<4) {  
    if( check_reconstruct3() ) {
      std::cout << GridLogMessage << "reconstruct3: PASS" << std::endl;
    } 
    
    if( checkGauge3x2mungers() ) {
      std::cout << GridLogMessage << "Gauge3x2mungers: PASS" << std::endl;
    }
  }

  if( check_reconstructSU() ) {
    std::cout << GridLogMessage << "reconstructSU: PASS" << std::endl;
  }
  if( checkGaugeSUmungers() ) {
    std::cout << GridLogMessage << "GaugeSUmungers: PASS" << std::endl;
  }

  if constexpr(Nc>2 && Nc%2==0) {
    std::cout <<GridLogMessage<< "Testing Sp(" << Nc << ") mungers..." << std::endl;

    if( check_reconstructSp() ) {
      std::cout << GridLogMessage << "reconstructSp: PASS" << std::endl;
    }

    if ( checkGaugeSpmungers() ) {
      std::cout << GridLogMessage << "GaugeSpmungers: PASS" << std::endl;
    }

  }

  std::cout << GridLogMessage << "Testing BinarySimple mungers..." << std::endl;
  if( checkBinarySimpleMungers() ) {
    std::cout << GridLogMessage << "BinarySimpleMungers: PASS" << std::endl;
  }

  std::cout << GridLogMessage << "Testing GaugeSimple mungers..." << std::endl;
  if( checkGaugeSimpleMungers() ) {
    std::cout << GridLogMessage << "GaugeSimpleMungers: PASS" << std::endl;
  }

  std::cout << GridLogMessage << "Testing GaugeDoubleStored mungers..." << std::endl;
  if( checkGaugeDoubleStoredMungers() ) {
    std::cout << GridLogMessage << "GaugeDoubleStoredMungers: PASS" << std::endl;
  }

  Grid_finalize();
#endif
}
