/*
MRGComponent.cc for ISO C++
version 1.00

authors: Hicham Wahbi
         Frédérik Rozon
         Richard Simard
*/

#include "latmrg/MRGComponent.h"
#include "latcommon/Util.h"
#include "latcommon/Const.h"
#include "latmrg/IntPrimitivity.h"
#include "latmrg/PolyPE.h"
#include "latcommon/IntFactor.h"

#include <fstream>

using namespace std;
using namespace LatCommon;


namespace LatMRG
{


MRGComponent::MRGComponent (const MScal & m, const MVect & a0, int k0)
{
   PolyPE::setM(m);
   module.init(m);
   k = k0;
   a.SetLength(k + 1);
   CopyVect(a0, a, k);
   orbitSeed.SetLength(k + 1);
}


//===========================================================================


MRGComponent::MRGComponent (long p, long e, long c, const MVect & a0, int k0)
{
   module.init(p, e, c);
   PolyPE::setM(getM());
   k = k0;
   a.SetLength(k + 1);
   CopyVect(a0, a, k);
   orbitSeed.SetLength(k + 1);
}


//===========================================================================

MRGComponent::MRGComponent (const MRGComponent & lat) : k(lat.k)
//   , ifm1(lat.ifm1), ifr(lat.ifr)
{
   module = lat.module;
   nj = lat.nj;
   rho = lat.rho;
//   a.kill();
   a.SetLength(k + 1);
   CopyVect(lat.a, a, k);
//   orbitSeed.kill();
   orbitSeed.SetLength(k + 1);
   CopyVect(lat.orbitSeed, orbitSeed, k);
}


//===========================================================================

MRGComponent & MRGComponent::operator= (const MRGComponent & lat)
{
   if (this != &lat) {
      k = lat.k;
      module = lat.module;
      nj = lat.nj;
      rho = lat.rho;
//    a.kill();
      a.SetLength(k + 1);
      CopyVect(lat.a, a, k);
 //     orbitSeed.kill();
      orbitSeed.SetLength(k + 1);
      CopyVect(lat.orbitSeed, orbitSeed, k);
//      ifm1 = lat.ifm1;
//      ifr = lat.ifr;
   }
   return *this;
}


//===========================================================================

void MRGComponent::init (const MScal & m0, int k0, DecompType decom1,
          const char *filem1, DecompType decor, const char *filer)
{
   PolyPE::setM(m0);
   module.init(m0);
   k = k0;
   a.SetLength(k + 1);
   orbitSeed.SetLength(k + 1);

   MScal m1;
   m1 = getM() - 1;
   ifm1.setNumber (m1);

   if (decom1 == DECOMP_READ)
      ifm1.read (filem1);
   else if (decom1 == DECOMP)
      ifm1.factorize();
   else if (decom1 == DECOMP_WRITE) {
      ifm1.factorize();
      ofstream fout(filem1);
      fout << ifm1.toString();
   }
   ifm1.calcInvFactors();

   MScal r;
   r = (power(m0, k) - 1) / (m0 - 1);
   ifr.setNumber(r);

   if (decor == DECOMP_READ)
      ifr.read (filer);
   else if (decor == DECOMP)
      ifr.factorize();
   else if (decor == DECOMP_WRITE) {
      ifr.factorize();
      ofstream fout(filer);
      fout << ifr.toString();
   } else if (decor == DECOMP_PRIME)
      ifr.setStatus (PRIME);

   ifr.calcInvFactors();
}


//===========================================================================

MRGComponent::MRGComponent (const MScal & m, int k, DecompType decom1,
      const char *filem1, DecompType decor, const char *filer)
{
   init (m, k, decom1, filem1, decor, filer);
}


//===========================================================================

MRGComponent::MRGComponent (Modulus & modu, int k, DecompType decom1,
      const char *filem1, DecompType decor, const char *filer)
{
   init (modu.m, k, decom1, filem1, decor, filer);

   PrimeType status = IntFactor::isPrime (modu.m, 100);
   if (status == PRIME || PROB_PRIME == status) {
      modu.primeF = true;
   } else {
      cout << " WARNING:  m is NOT prime" << endl;
      modu.primeF = false;
   }
   module = modu;
}


//===========================================================================

MRGComponent::~MRGComponent()
{
 //  a.kill();
 //  orbitSeed.kill();
}


//===========================================================================

void MRGComponent::setA (const MVect & b)
{
   a = b;
   //  CopyVect(b, a, k);
}


//===========================================================================

bool MRGComponent::maxPeriod (const MVect & a0)
{
   PolyPE::setM(getM());
   a = a0;
   PolyPE::reverse (a, k, 2);
   PolyPE::setF(a);
   PolyPE pol;
   IntPrimitivity privfm(ifm1, getM(), 1);
   return pol.isPrimitive(privfm, ifr);
}


//===========================================================================

bool MRGComponent::maxPeriod23 (const MVect & a0)
{
   PolyPE::setM(getM());
   a = a0;
   PolyPE::reverse (a, k, 2);
   PolyPE::setF(a);
   PolyPE pol;
   // La condition 1 a déjà été vérifiée dans SeekMain
   return pol.isPrimitive(ifr);
}


//===========================================================================

string MRGComponent::toString () 
{
   ostringstream os;
   os << "MRGComponent:";
   MScal mm = getM();
   os << "\n   m = " << mm;
   os << "\n   k = " << k;
   os << "\n   a = ";
   string str (os.str ());
   string s2 = LatCommon::toString(a, k);
   str += s2;
   str += "\n";
   return str;
}


//===========================================================================

}
