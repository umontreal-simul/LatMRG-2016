#include <sstream>
#include <iostream>
#include <string>
#include <stdexcept>

#include "NTL/ZZ.h"
#include "NTL/ZZ_p.h"
#include "NTL/ZZ_pX.h"
#include "NTL/ZZ_pXFactoring.h"

#include "NTL/lzz_p.h"
#include "NTL/lzz_pE.h"
#include "NTL/lzz_pX.h"
#include "NTL/lzz_pXFactoring.h"

#include "latcommon/Util.h"
#include "latcommon/IntFactor.h"
#include "latmrg/PolyPE.h"

using namespace NTL;
using namespace std;


namespace LatMRG
{

/*=========================================================================*/

MScal PolyPE::m_m;
long PolyPE::m_k;
PolX PolyPE::m_x;



/*=========================================================================*/
//PolyPE::PolyPE (const MVect & C, long n)
//{
//   init (C, n);
//}


/*=========================================================================*/

PolyPE::PolyPE ()
{}


/*=========================================================================*/

void PolyPE::setM (const MScal & m)
{
   MScalP::init (m);
   m_m = m;
}


/*=========================================================================*/

void PolyPE::setF (const MVectP & c)
{
   m_k = c.length () - 1;
   PolX f;
   for (long i = 0; i <= m_k; i++)
      SetCoeff (f, i, c[i]);
   f.normalize ();
   PolE::init (f);
   string str = "[0 1]";
   istringstream in (str);
   m_x.kill();
   in >> m_x;
}

/*=========================================================================*/

void PolyPE::setF (const MVect & V)
{
   MVectP vP;
   conv (vP, V);
   setF (vP);
}


/*=========================================================================*/

void PolyPE::setVal (long j)
{
   MScalP coeff;
   conv (coeff, 1);
   PolX val;
   SetCoeff (val, j, coeff);
   val.normalize ();
   init (val);
}


/*=========================================================================*/

void PolyPE::setVal (const MVect & C)
{
   ostringstream out;
   out << "[";
   for (int i = 0; i < C.length (); i++) {
      out << C[i] << " ";
   }
   out << "]";
   std::string str = out.str ();
   setVal (str);
}


/*=========================================================================*/

void PolyPE::powerMod (const MScal & j)
{
   string str = "[0 1]";
   istringstream in (str);
   PolyPE A;
   in >> A;
   power (*this, A, j);
}


/*=========================================================================*/

void PolyPE::setVal (std::string & str)
{
   istringstream in (str);
   in >> *this;
}


//===========================================================================

std::string PolyPE::toString ()const
{
   std::ostringstream sortie;

   PolX pX = PolE::modulus ().val ();
   sortie << "m = " << PolyPE::m_m << endl;
   sortie << "k = " << getK () << endl;
   sortie << "f = " << getF () << endl;
   sortie << "v = " << *this << endl;
   sortie << "    " << endl;
   return sortie.str ();
}


/*=========================================================================*/

void PolyPE::toVector (MVect & C)
{
   // C = rep(*this);
   for (int i = 0; i <= getK (); ++i) {
      C[i] = rep (coeff (rep (*this), i));
   }
}


/*=========================================================================*/
#if 0
bool PolyPE::isIrreducible ()
{
   // Méthode de crandall-Pomerance. La vitesse de cette méthode est
   // pratiquement identique à celle de NTL::DetIrredTest, mais plus
   // lente que NTL::IterIrredTest, spécialement pour grand k.
   PolE g;
   PolX d;
   string str = "[0 1]";
   istringstream in (str);
   in >> g;
   for (int i = 1; i <= m_k/2; ++i) {
      power(g, g, m_m);
      GCD (d, PolE::modulus(), rep(g) - m_x);
      if (1 != IsOne(d))
         return false;
   }
   return true;
}
#endif

/*=========================================================================*/

bool PolyPE::isPrimitive (const IntFactorization & r)
{
   if (1 == getK())
      return true;

   // Is f irreducible ?
   PolX Q;
   Q = getF ();
   //  if (!isIrreducible())      // slow
   //  if (0 == DetIrredTest(Q))   // medium slow
   if (0 == IterIrredTest(Q))   // fastest
      return false;

   // ---- Condition 2
   MScal r0;
   r0 = r.getNumber ();
   Q = PowerXMod (r0, getF ());
   if (0 != deg (Q))
      return false;

   MScal T1;
   T1 = rep (ConstTerm (getF ()));
   if ((getK () & 1) == 1)
      T1 = -T1;
   if (T1 < 0)
      T1 += getM ();
   if (rep (ConstTerm (Q)) != T1)
      return false;

   // ---- Condition 3
   if (r.getStatus () == LatCommon::PRIME)
      return true;

   std::vector <MScal> invFactorList = r.getInvFactorList ();
   //   assert (!invFactorList.empty ());
   std::vector <MScal>::const_iterator it = invFactorList.begin ();

   while (it != invFactorList.end ()) {
      Q = PowerXMod (*it, getF());
      if (0 == deg (Q))
         return false;
      ++it;
   }
   return true;
}


/*=========================================================================*/

bool PolyPE::isPrimitive (const IntPrimitivity & fm,
                          const IntFactorization & fr)
{
   MScal a0;
   a0 = -rep (ConstTerm (getF ()));
   if ((getK () & 1) == 0)
      a0 = -a0;
   if (!fm.isPrimitiveElement (a0))
      return false;
   return isPrimitive (fr);
}


/*=========================================================================*/

void PolyPE::reverse (MVect & C, long n, int kind)
{
   long i;
   MScal temp;
   if ((kind == 2) || (kind == 1)) {
      for (i = 0; i < (n + 1) / 2; i++) {
         temp = C[n - i];
         C[n - i] = C[i];
         C[i] = temp;
      }
   }
   if (kind == 2) {
      for (i = 0; i < n; i++)
         C[i] = -C[i];
      C[n] = 1;
   }
}

}
