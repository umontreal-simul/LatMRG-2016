#include "latmrg/Primes.h"
#include "latcommon/Types.h"
#include "latcommon/Util.h"
#include "latcommon/IntFactor.h"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <climits>


using namespace std;
using namespace NTL;
using namespace LatCommon;

namespace
{
MScal One, Two;
MScal m, m1, m1s2, r, Sm0, Sm1, Sm2, Sdiff;
}


namespace LatMRG
{

//===========================================================================

Primes::Primes ()
{
   set9(One);
   conv (Two, 2);
   timer.init();
}


//===========================================================================

Primes::~Primes ()
{
}


//===========================================================================
/*
inline void Primes::nextM (MScal & m)
{
   m -= 2;
   if ((0 == m % 5) && (m > 5))
      m -= 2;
}
*/

//===========================================================================

void Primes::find (int k, int e, int s, const MScal & S1, const MScal & S2,
     bool safe, bool facto, ofstream & fout)
{
   set9 (r);
   if (IsOdd (S2))
      m = S2;
   else
      m = S2 - One;
   int i = 0;
   const long KTRIALS = 200;

   while (i < s && m >= S1) {
      // m1 = m-1
      // m1s2 = (m-1)/2
      // r = (m^k-1)/(m-1)
      PrimeType status = IntFactor::isPrime (m, KTRIALS);
      if (status == PRIME || status == PROB_PRIME) {
         m1 = m - One;
         if (safe) {
            if (1 == m % 4) {
               nextM(m);
               continue;
            }
            Quotient (m1, Two, m1s2);
            status = IntFactor::isPrime (m1s2, KTRIALS);
            if (status != PRIME && status != PROB_PRIME) {
               nextM(m);
               continue;
            }
         }
         if (k > 1) {
            r = power (m, k);
            --r;
            Quotient (r, m1, r);
            status = IntFactor::isPrime (r, KTRIALS);
         }
         if (k == 1 || status == PRIME || status == PROB_PRIME) {
            i++;
            fout << "   m = " << m << endl;
            Sdiff = m - Sm0;
            fout << "     = 2^" << e;
            if (Sdiff >= 0) {
               fout << " + ";
               fout << Sdiff << endl;
            } else {
               fout << " - ";
               fout << (-Sdiff) << endl;
            }
            if (facto) {
               ifac.clear();
               ifac.setNumber (m1);
               ifac.factorize ();
               fout << "   Factors of m - 1:\n";
               fout << ifac.toString ();
            }
         }
      }
      nextM(m);
   }
}


//===========================================================================

void Primes::find (int k, int e, int s, bool safe, bool facto,
                   std::ofstream & fout)
{
   writeHeader (k, e, INT_MAX, INT_MAX, safe, facto, fout);
   timer.init();
   Sm0 = One << e;
   Sm2 = Sm0 - 1;
   Sm1 = Sm0 >> 1;
   find (k, e, s, Sm1, Sm2, safe, facto, fout);
   writeFooter (fout);
}


//===========================================================================

void Primes::find (int e, int s, bool facto, std::ofstream & fout)
{
   writeHeader (1, e, INT_MAX, INT_MAX, false, facto, fout);
   timer.init();
   Sm0 = One << e;
   Sm2 = Sm0 - 1;
   Sm1 = Sm0 >> 1;
   find (1, e, s, Sm1, Sm2, false, facto, fout);
   writeFooter (fout);
}


//===========================================================================

void Primes::find (int k, int e, long c1, long c2, bool safe,
                   bool facto, ofstream & fout)
{
   writeHeader (k, e, c1, c2, safe, facto, fout);
   timer.init();
   Sm0 = One << e;
   Sm1 = Sm0 + c1;
   Sm2 = Sm0 + c2;
   find (k, e, INT_MAX, Sm1, Sm2, safe, facto, fout);
   writeFooter (fout);
}


//===========================================================================

void Primes::writeHeader (int k, int e, long c1, long c2, bool safe,
                   bool facto, ofstream & fout)
{
   fout << "-----------------------------------------------------" << endl;
   fout << "Values such that m";
   if (safe)
      fout << ", (m-1)/2,";
   if (k > 1)
      fout << " and (m^k-1)/(m-1) are prime.\n";
   else
      fout << " is prime.\n";
   fout << "k  = " << k << endl;
   fout << "e  = " << e << endl;
   if (c1 < INT_MAX)
      fout << "c1 = " << c1 << endl;
   if (c2 < INT_MAX)
      fout << "c2 = " << c2 << endl;
   fout << "safe = " << boolalpha << safe << endl;
   fout << "facto = " << facto << endl;

   if (c1 < INT_MAX) {
      fout << "\nm from 2^" << e;
      if (c1 >= 0)
         fout << " + ";
      else
         fout << " ";
      fout << c1;
      fout << " to 2^" << e;
      if (c2 >= 0)
         fout << " + ";
      else
         fout << " ";
      fout << c2 << endl << endl;

   } else {
      fout << "\nm < 2^" << e << endl << endl;
   }
}


//===========================================================================

void Primes::writeFooter (ofstream & fout)
{
   fout << "\nCPU time: ";
   fout << timer.toString () << endl;
   char *machine;
   machine = getenv ("HOST");
   fout << "Tests made on " << machine << endl << endl;
}

//===========================================================================

}
