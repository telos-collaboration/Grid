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

using namespace Grid;

// tests to check the various (un)munger classes 
// and functions in parallelIO/Metadata.h
bool check_reconstruct3();
bool check_reconstructSU();
bool check_is_perm_even();
bool check_unique_reconstructSU();
bool checkGauge3x2mungers();
bool checkGaugeSUmungers();
bool check_reconstructSp();
bool checkGaugeSpmungers();
bool checkBinarySimpleMungers();
bool checkGaugeSimpleMungers();
bool checkGaugeDoubleStoredMungers();

void mock_SU_field();
void mock_Sp_field();

//////////////////////////////////////////////
// helper functions for generating (mock) SU/Sp matrices
//////////////////////////////////////////////
void mock_SU_field( ColourMatrixD &cm, std::vector<int> seed )  {

  GridSerialRNG sRNG;   sRNG.SeedFixedIntegers(seed);

  ColourMatrixD lie, la, ta;
  ComplexD ca;

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

void mock_Sp_field( ColourMatrixD &cm, std::vector<int> seed ) {

  GridSerialRNG sRNG;   sRNG.SeedFixedIntegers(seed);

  std::vector<ComplexD> elem(Nc*Nc/2);
  for(auto&& e: elem) {
    random(sRNG,e);
  }

  // fill upper left and upper right blocks only
  for(int i=0;i<Nc/2;i++) {
    for(int j=0;j<Nc/2;j++) {
      cm()()(i,j) = elem[i*Nc+j];
      cm()()(i,j+Nc/2) = elem[i*Nc+j+Nc/2];
    }
  }

}
////////////////////////////////////////////

bool check_reconstruct3() {

  ColourMatrixD SU3_xfield, SU3_yfield, SU3_zfield, SU3_tfield;

  mock_SU_field(SU3_xfield, std::vector<int>({30}));
  mock_SU_field(SU3_yfield, std::vector<int>({95}));
  mock_SU_field(SU3_zfield, std::vector<int>({7}));
  mock_SU_field(SU3_tfield, std::vector<int>({10}));

  //set last row equal to zero
  for(int j=0;j<Nc;j++) {
    SU3_xfield()()(Nc-1,j) = 0.0;
    SU3_yfield()()(Nc-1,j) = 0.0;
    SU3_zfield()()(Nc-1,j) = 0.0;
    SU3_tfield()()(Nc-1,j) = 0.0;
  }

  LorentzColourMatrixD scalar = Zero();

  pokeLorentz(scalar, SU3_xfield, 0);
  pokeLorentz(scalar, SU3_yfield, 1);
  pokeLorentz(scalar, SU3_zfield, 2);
  pokeLorentz(scalar, SU3_tfield, 3);

  reconstruct3(scalar);

  // check result is unitary and det==1
  for(int mu=0;mu<Nd;mu++) {
    auto new_cm = peekIndex<LorentzIndex>(scalar, mu);
    auto det = Determinant(new_cm);
    assert( abs(norm2(det)-1.0) < 1e-12 );
    assert( norm2( new_cm*adj(new_cm)-1.0 ) < 1e-12 );
  }

  return true;

}

bool check_reconstructSU() {

  ColourMatrixD SUN_xfield, SUN_yfield, SUN_zfield, SUN_tfield;

  mock_SU_field(SUN_xfield, std::vector<int>({4}));
  mock_SU_field(SUN_yfield, std::vector<int>({8}));
  mock_SU_field(SUN_zfield, std::vector<int>({15}));
  mock_SU_field(SUN_tfield, std::vector<int>({16}));

  //set last row equal to zero
  for(int j=0;j<Nc;j++) {
    SUN_xfield()()(Nc-1,j) = 0.0;
    SUN_yfield()()(Nc-1,j) = 0.0;
    SUN_zfield()()(Nc-1,j) = 0.0;
    SUN_tfield()()(Nc-1,j) = 0.0;
  }

  LorentzColourMatrixD scalar = Zero();

  pokeLorentz(scalar, SUN_xfield, 0);
  pokeLorentz(scalar, SUN_yfield, 1);
  pokeLorentz(scalar, SUN_zfield, 2);
  pokeLorentz(scalar, SUN_tfield, 3);

  reconstructSU(scalar);

  // check result is unitary and det==1
  for(int mu=0;mu<Nd;mu++) {
    auto new_cm = peekIndex<LorentzIndex>(scalar, mu);
    auto det = Determinant(new_cm);
    assert( abs( norm2(det)-1.0 ) < 1e-12 );
    assert( norm2( new_cm*adj(new_cm)-1.0 ) < 1e-12 );
  }

  return true;

}

bool check_is_perm_even() {

  std::vector<int> v = {0,1,2,3};
  int count = 0;

  // swap two elements and check if parity flipped.
  for(auto& e: v) {
    for(auto& d: v) {
      if(e!=d) {
        if(count%2==0) { assert( is_perm_even(v)==true ); }
        else { assert( is_perm_even(v)==false ); }
        std::swap(e, d);
        count += 1;
      }
    }
  }      
 
  return true;

}

bool check_unique_reconstructSU() {

  ColourMatrixD SUN_xfield, SUN_yfield, SUN_zfield, SUN_tfield;

  mock_SU_field(SUN_xfield, std::vector<int>({2}));
  mock_SU_field(SUN_yfield, std::vector<int>({40}));
  mock_SU_field(SUN_zfield, std::vector<int>({28}));
  mock_SU_field(SUN_tfield, std::vector<int>({25}));

  //set last row equal to zero
  for(int j=0;j<Nc;j++) {
    SUN_xfield()()(Nc-1,j) = 0.0;
    SUN_yfield()()(Nc-1,j) = 0.0;
    SUN_zfield()()(Nc-1,j) = 0.0;
    SUN_tfield()()(Nc-1,j) = 0.0;
  }

  LorentzColourMatrixD scalar = Zero();

  pokeLorentz(scalar, SUN_xfield, 0);
  pokeLorentz(scalar, SUN_yfield, 1);
  pokeLorentz(scalar, SUN_zfield, 2);
  pokeLorentz(scalar, SUN_tfield, 3);

  unique_reconstructSU(scalar);

  // check result is unitary and det==1
  for(int mu=0;mu<Nd;mu++) {
    auto new_cm = peekIndex<LorentzIndex>(scalar, mu);
    auto det = Determinant(new_cm);
    assert( abs( norm2(det)-1.0 ) < 1e-12 );
    assert( norm2( new_cm*adj(new_cm)-1.0 ) < 1e-12 );
  }

  return true;

}

bool checkGauge3x2mungers() {

  ColourMatrixD SU3_xfield, SU3_yfield, SU3_zfield, SU3_tfield;

  mock_SU_field(SU3_xfield, std::vector<int>({66}));
  mock_SU_field(SU3_yfield, std::vector<int>({63}));
  mock_SU_field(SU3_zfield, std::vector<int>({71}));
  mock_SU_field(SU3_tfield, std::vector<int>({98}));

  LorentzColourMatrixD scalar = Zero();

  pokeLorentz(scalar, SU3_xfield, 0);
  pokeLorentz(scalar, SU3_yfield, 1);
  pokeLorentz(scalar, SU3_zfield, 2);
  pokeLorentz(scalar, SU3_tfield, 3);

  // reduce field. Gauge3x2unmunger<out_type,in_type>
  Gauge3x2unmunger<LorentzColour2x3D,LorentzColourMatrixD> unmunger;
  LorentzColour2x3D scalar_rdc;
  unmunger(scalar,scalar_rdc);

  // reconstruct full field. Gauge3x2munger<in_type,out_type>
  Gauge3x2munger<LorentzColour2x3D,LorentzColourMatrixD> munger;
  LorentzColourMatrixD scalar_recon;
  munger(scalar_rdc,scalar_recon);

  // round-trip test
  assert(norm2(scalar_recon-scalar)<1e-10);

  return true;

}

template<bool unique_su>
bool checkGaugeSUmungers() {

  ColourMatrixD SUN_xfield, SUN_yfield, SUN_zfield, SUN_tfield;

  mock_SU_field(SUN_xfield, std::vector<int>({360}));
  mock_SU_field(SUN_yfield, std::vector<int>({804}));
  mock_SU_field(SUN_zfield, std::vector<int>({77}));
  mock_SU_field(SUN_tfield, std::vector<int>({65}));

  LorentzColourMatrixD scalar = Zero();

  pokeLorentz(scalar, SUN_xfield, 0);
  pokeLorentz(scalar, SUN_yfield, 1);
  pokeLorentz(scalar, SUN_zfield, 2);
  pokeLorentz(scalar, SUN_tfield, 3);

  // reduce field. GaugeSUunmunger<out_type,in_type>
  GaugeSUunmunger<LorentzColour2x3D,LorentzColourMatrixD> unmunger;
  LorentzColour2x3D scalar_rdc = Zero();
  unmunger(scalar,scalar_rdc);

  // reconstruct full field. GaugeSUmunger<in_type,out_type>
  GaugeSUmunger<LorentzColour2x3D,LorentzColourMatrixD,unique_su> munger;
  LorentzColourMatrixD scalar_recon;
  munger(scalar_rdc, scalar_recon);
  
  // round-trip test
  assert(norm2(scalar_recon-scalar)<1e-12);

  return true;

}


bool check_reconstructSp() {

  ColourMatrixD Sp_xfield, Sp_yfield, Sp_zfield, Sp_tfield;
  Sp_xfield = Sp_yfield = Sp_zfield = Sp_tfield = Zero();

  mock_Sp_field(Sp_xfield, std::vector<int>({1066}));
  mock_Sp_field(Sp_yfield, std::vector<int>({1995}));
  mock_Sp_field(Sp_zfield, std::vector<int>({997}));
  mock_Sp_field(Sp_tfield, std::vector<int>({58}));

  LorentzColourMatrixD scalar = Zero();

  pokeLorentz(scalar, Sp_xfield, 0);
  pokeLorentz(scalar, Sp_yfield, 1);
  pokeLorentz(scalar, Sp_zfield, 2);
  pokeLorentz(scalar, Sp_tfield, 3);

  reconstructSp(scalar);

  // test Sp block structure 
  //  A  B
  // -B* A*
  for(int mu=0; mu<Nd; mu++) {
    auto Sp_cm = peekIndex<LorentzIndex>(scalar,mu);
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

  ColourMatrixD Sp_xfield, Sp_yfield, Sp_zfield, Sp_tfield;
  Sp_xfield = Sp_yfield = Sp_zfield = Sp_tfield = Zero();

  mock_Sp_field(Sp_xfield, std::vector<int>({23}));
  mock_Sp_field(Sp_yfield, std::vector<int>({42}));
  mock_Sp_field(Sp_zfield, std::vector<int>({93}));
  mock_Sp_field(Sp_tfield, std::vector<int>({49}));

  LorentzColourMatrixD scalar = Zero();

  pokeLorentz(scalar, Sp_xfield, 0);
  pokeLorentz(scalar, Sp_yfield, 1);
  pokeLorentz(scalar, Sp_zfield, 2);
  pokeLorentz(scalar, Sp_tfield, 3);

  reconstructSp(scalar);

  // reduce field. GaugeSpunmunger<out_type,in_type>
  GaugeSpunmunger<LorentzColourNx2ND,LorentzColourMatrixD> unmunger;
  LorentzColourNx2ND scalar_rdc = Zero();
  unmunger(scalar,scalar_rdc);

  // reconstruct full field. GaugeSpmunger<in_type,out_type>
  GaugeSpmunger<LorentzColourNx2ND,LorentzColourMatrixD> munger;
  LorentzColourMatrixD scalar_recon;
  munger(scalar_rdc,scalar_recon);

  // round-trip test
  assert(scalar==scalar_recon);

  return true;

}


bool checkBinarySimpleMungers() {

  LorentzColourMatrixF in_scalar_objectF;  // single precision
  LorentzColourMatrixF out_scalar_objectF = Zero();  
  LorentzColourMatrixD in_scalar_objectD;  // double precision
  LorentzColourMatrixD out_scalar_objectD = Zero();  

  GridSerialRNG sRNG;   sRNG.SeedFixedIntegers(std::vector<int>({117}));

  random(sRNG,in_scalar_objectF);
  random(sRNG,in_scalar_objectD);

  // BinarySimpleUnmunger<out_type,in_type> 
  BinarySimpleUnmunger<LorentzColourMatrixF,LorentzColourMatrixF> unmungerFF;
  BinarySimpleUnmunger<LorentzColourMatrixD,LorentzColourMatrixF> unmungerFD;
  BinarySimpleUnmunger<LorentzColourMatrixD,LorentzColourMatrixD> unmungerDD;
  BinarySimpleUnmunger<LorentzColourMatrixF,LorentzColourMatrixD> unmungerDF;

  unmungerFF(in_scalar_objectF, out_scalar_objectF);
  unmungerDD(in_scalar_objectD, out_scalar_objectD);
  assert(in_scalar_objectF==out_scalar_objectF);
  assert(in_scalar_objectD==out_scalar_objectD);

  random(sRNG,in_scalar_objectF);
  unmungerFD(in_scalar_objectF, out_scalar_objectD);
  unmungerDF(out_scalar_objectD, out_scalar_objectF);
  assert(in_scalar_objectF==out_scalar_objectF);

  random(sRNG,in_scalar_objectD);
  unmungerDF(in_scalar_objectD, out_scalar_objectF);
  unmungerFD(out_scalar_objectF, out_scalar_objectD);
  // lose exactness when going from double-->single precision
  assert(norm2( in_scalar_objectD - out_scalar_objectD ) < 1e-10);

  // BinarySimpleMunger<in_type,out_type> 
  BinarySimpleMunger<LorentzColourMatrixF,LorentzColourMatrixF> mungerFF;
  BinarySimpleMunger<LorentzColourMatrixD,LorentzColourMatrixF> mungerDF;
  BinarySimpleMunger<LorentzColourMatrixD,LorentzColourMatrixD> mungerDD;
  BinarySimpleMunger<LorentzColourMatrixF,LorentzColourMatrixD> mungerFD;

  random(sRNG,in_scalar_objectF);
  random(sRNG,in_scalar_objectD);

  mungerFF(in_scalar_objectF, out_scalar_objectF);
  mungerDD(in_scalar_objectD, out_scalar_objectD);
  assert(in_scalar_objectD==out_scalar_objectD);
  assert(in_scalar_objectF==out_scalar_objectF);
 
  random(sRNG,in_scalar_objectF);
  mungerFD(in_scalar_objectF, out_scalar_objectD);
  mungerDF(out_scalar_objectD, out_scalar_objectF);
  assert(in_scalar_objectF==out_scalar_objectF);

  random(sRNG,in_scalar_objectD);
  mungerDF(in_scalar_objectD, out_scalar_objectF);
  mungerFD(out_scalar_objectF, out_scalar_objectD);
  // lose exactness when going from double-->single precision
  assert(norm2( in_scalar_objectD - out_scalar_objectD ) < 1e-10);

  return true;

}

// these are used in NerscIO.h
bool checkGaugeSimpleMungers() {

  LorentzColourMatrixF in_scalar_objectF;  // single precision
  LorentzColourMatrixF out_scalar_objectF = Zero(); 
  LorentzColourMatrixD in_scalar_objectD;  // double precision
  LorentzColourMatrixD out_scalar_objectD = Zero(); 
  
  GridSerialRNG sRNG;   sRNG.SeedFixedIntegers(std::vector<int>({57}));

  random(sRNG,in_scalar_objectF);
  random(sRNG,in_scalar_objectD);

  // GaugeSimpleUnmunger<out_type,in_type> 
  GaugeSimpleUnmunger<LorentzColourMatrixF,LorentzColourMatrixF> unmungerFF;
  GaugeSimpleUnmunger<LorentzColourMatrixD,LorentzColourMatrixF> unmungerFD;
  GaugeSimpleUnmunger<LorentzColourMatrixD,LorentzColourMatrixD> unmungerDD;
  GaugeSimpleUnmunger<LorentzColourMatrixF,LorentzColourMatrixD> unmungerDF;

  unmungerFF(in_scalar_objectF, out_scalar_objectF);
  unmungerDD(in_scalar_objectD, out_scalar_objectD);
  assert(in_scalar_objectF==out_scalar_objectF);
  assert(in_scalar_objectD==out_scalar_objectD);

  random(sRNG,in_scalar_objectF);
  unmungerFD(in_scalar_objectF, out_scalar_objectD);
  unmungerDF(out_scalar_objectD, out_scalar_objectF);
  assert(in_scalar_objectF==out_scalar_objectF);

  random(sRNG,in_scalar_objectD);
  unmungerDF(in_scalar_objectD, out_scalar_objectF);
  unmungerFD(out_scalar_objectF, out_scalar_objectD);
  // lose exactness when going from double-->single precision
  assert(norm2( in_scalar_objectD - out_scalar_objectD ) < 1e-10);

  // GaugeSimpleMunger<in_type,out_type> 
  GaugeSimpleMunger<LorentzColourMatrixF,LorentzColourMatrixF> mungerFF;
  GaugeSimpleMunger<LorentzColourMatrixD,LorentzColourMatrixF> mungerDF;
  GaugeSimpleMunger<LorentzColourMatrixD,LorentzColourMatrixD> mungerDD;
  GaugeSimpleMunger<LorentzColourMatrixF,LorentzColourMatrixD> mungerFD;

  random(sRNG,in_scalar_objectF);
  random(sRNG,in_scalar_objectD);

  mungerFF(in_scalar_objectF, out_scalar_objectF);
  mungerDD(in_scalar_objectD, out_scalar_objectD);
  assert(in_scalar_objectD==out_scalar_objectD);
  assert(in_scalar_objectF==out_scalar_objectF);
 
  random(sRNG,in_scalar_objectF);
  mungerFD(in_scalar_objectF, out_scalar_objectD);
  mungerDF(out_scalar_objectD, out_scalar_objectF);
  assert(in_scalar_objectF==out_scalar_objectF);

  random(sRNG,in_scalar_objectD);
  mungerDF(in_scalar_objectD, out_scalar_objectF);
  mungerFD(out_scalar_objectF, out_scalar_objectD);
  // lose exactness when going from double-->single precision
  assert(norm2( in_scalar_objectD - out_scalar_objectD ) < 1e-10);

  return true;

}

bool checkGaugeDoubleStoredMungers() {

  DoubleStoredColourMatrixF in_scalar_objectF;  // single precision
  DoubleStoredColourMatrixF out_scalar_objectF = Zero(); 
  DoubleStoredColourMatrixD in_scalar_objectD;  // double precision
  DoubleStoredColourMatrixD out_scalar_objectD = Zero(); 
 
  GridSerialRNG sRNG;   sRNG.SeedFixedIntegers(std::vector<int>({169}));

  random(sRNG,in_scalar_objectF);
  random(sRNG,in_scalar_objectD);

  // GaugeDoubleStoredUnmunger<out_type,in_type> 
  GaugeDoubleStoredUnmunger<DoubleStoredColourMatrixF,DoubleStoredColourMatrixF> unmungerFF;
  GaugeDoubleStoredUnmunger<DoubleStoredColourMatrixD,DoubleStoredColourMatrixF> unmungerFD;
  GaugeDoubleStoredUnmunger<DoubleStoredColourMatrixD,DoubleStoredColourMatrixD> unmungerDD;
  GaugeDoubleStoredUnmunger<DoubleStoredColourMatrixF,DoubleStoredColourMatrixD> unmungerDF;

  unmungerFF(in_scalar_objectF, out_scalar_objectF);
  unmungerDD(in_scalar_objectD, out_scalar_objectD);
  assert(in_scalar_objectF==out_scalar_objectF);
  assert(in_scalar_objectD==out_scalar_objectD);

  random(sRNG,in_scalar_objectF);
  unmungerFD(in_scalar_objectF, out_scalar_objectD);
  unmungerDF(out_scalar_objectD, out_scalar_objectF);
  assert(in_scalar_objectF==out_scalar_objectF);

  random(sRNG,in_scalar_objectD);
  unmungerDF(in_scalar_objectD, out_scalar_objectF);
  unmungerFD(out_scalar_objectF, out_scalar_objectD);
  // lose exactness when going from double-->single precision
  assert(norm2( in_scalar_objectD - out_scalar_objectD ) < 1e-10);

  // GaugeDoubleStoredMunger<in_type,out_type> 
  GaugeDoubleStoredMunger<DoubleStoredColourMatrixF,DoubleStoredColourMatrixF> mungerFF;
  GaugeDoubleStoredMunger<DoubleStoredColourMatrixD,DoubleStoredColourMatrixF> mungerDF;
  GaugeDoubleStoredMunger<DoubleStoredColourMatrixD,DoubleStoredColourMatrixD> mungerDD;
  GaugeDoubleStoredMunger<DoubleStoredColourMatrixF,DoubleStoredColourMatrixD> mungerFD;
 
  random(sRNG,in_scalar_objectF);
  random(sRNG,in_scalar_objectD);

  mungerFF(in_scalar_objectF, out_scalar_objectF);
  mungerDD(in_scalar_objectD, out_scalar_objectD);
  assert(in_scalar_objectD==out_scalar_objectD);
  assert(in_scalar_objectF==out_scalar_objectF);
 
  random(sRNG,in_scalar_objectF);
  mungerFD(in_scalar_objectF, out_scalar_objectD);
  mungerDF(out_scalar_objectD, out_scalar_objectF);
  assert(in_scalar_objectF==out_scalar_objectF);

  random(sRNG,in_scalar_objectD);
  mungerDF(in_scalar_objectD, out_scalar_objectF);
  mungerFD(out_scalar_objectF, out_scalar_objectD);
  // lose exactness when going from double-->single precision
  assert(norm2( in_scalar_objectD - out_scalar_objectD ) < 1e-10);

  return true;

}


int main (int argc, char ** argv)
{
#ifdef HAVE_LIME
  Grid_init(&argc,&argv);
  
  std::cout <<GridLogMessage<< " main "<<std::endl;

  if( check_is_perm_even() ) {
    std::cout << GridLogMessage << "is_perm_even: PASS" << std::endl;
  }

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
  if( checkGaugeSUmungers<false>() ) {
    std::cout << GridLogMessage << "(unique_su=false) GaugeSUmungers: PASS" << std::endl;
  }

  if( check_unique_reconstructSU() ) {
    std::cout << GridLogMessage << "unique_reconstructSU: PASS" << std::endl;
  }
  if( checkGaugeSUmungers<true>() ) {
    std::cout << GridLogMessage << "(unique_su=true) GaugeSUmungers: PASS" << std::endl;
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
