#include "latmrg/Modulus.h"
#include "latcommon/Util.h"

#include <cassert>


using namespace std;
using namespace LatCommon;


namespace LatMRG
{

Modulus::Modulus()
{}


//===========================================================================

Modulus::Modulus(const MScal & j)
{
   init(j);
}


//===========================================================================

Modulus::Modulus (long m1, long m2, long m3)
{
   init(m1, m2, m3);
}


//===========================================================================

Modulus::~Modulus()
{}


//===========================================================================

void Modulus::init (const MScal & j)
{
   mRed = m = j;
   threeF = false;
   b = 0;
   e = 0;
   c = 0;
   mRac = (MScal) SqrRoot (m);
   mRacNeg = -mRac;
   Eight = 8;
   Four = 4;
}


//===========================================================================

void Modulus::init (long m1, long m2, long m3)
{
   assert (m1 > 1);
   m = m1;
   if (m2 <= 0) {
      init (m);
      return;
   }

   m = power (m, m2) + m3;
   init (m);
   b = m1;
   e = m2;
   c = m3;
   threeF = true;
   bm1 = b - 1;
   conv (Y, b);
   b2 = Y*Y;
}


//===========================================================================

bool Modulus::perMaxPowPrime (const MScal & A)
{
   if (!threeF || c != 0) {
      MyExit(1, "perMaxPowPrime:   m must be a power of a prime");
      return false;
   }

   if (b == 2) {
      Modulo (A, Eight, Y);
      return (Y == 5) || (Y == 3);
   } else {
      Y = A % b2;
      Y = PowerMod (Y, bm1, b2);
      return (Y != 1);
   }
}


//===========================================================================

void Modulus::reduceM (const MScal & a)
{
   if (!threeF || c != 0)
      return;

   MScal Y2, Y3, Y4;

   // We now assume that m is a power of a prime b
   if (b == 2) {
      // m is a power of 2
      Modulo (a, Eight, Y);
      if (Y == 5)
         Quotient (m, Four, mRed);
      else if (Y == 3)
         Quotient (m, Eight, mRed);
      else {
         Y2 = a;
         Y3 = Eight;
         if (Y == 1)
            --Y2;
         else if (Y == 7)
            ++Y2;
         else
            MyExit (1, "Power of 2 modulus with even multiplier");
         do {
            Y3 = Y3 + Y3;
            Modulo (Y2, Y3, Y4);
         } while (Y4 == 0);
         Quotient (m, Y3, mRed);
      }

   } else {
      // mj is a power of a prime b > 2
      conv (Y3, b);
      if (perMaxPowPrime (a))
         Divide (mRed, Y4, m, Y3);
      else {
         Y2 = b2;
         do {
            Y2 *= b;        // Y = b^e
            Y = a % Y2;
            Y = PowerMod (Y, bm1, Y2);
         } while (Y == 1);
         Divide (Y2, Y4, Y2, Y3);
         Quotient (m, Y2, mRed);
      }
   }
}


//===========================================================================

}
