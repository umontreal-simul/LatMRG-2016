#include "latcommon/Util.h"
#include "latmrg/SearcherKorobov.h"

#include <cfloat>

using namespace std;
using namespace LatCommon;


namespace LatMRG
{

//===========================================================================

void SearcherKorobov::print (long y)
{
   // For debugging
   cout << "  a = " << y << endl;
}


//===========================================================================

inline void SearcherKorobov::calcAs (long n, int s, long a, std::vector<long>& z)
{
   // assert (n <= 2147483647);
   z.assign(s + 1, 1L);
   z[1] = 1;
   for (int j = 2; j <= s; ++j)
      z[j] = (a * z[j - 1]) % n;
}


//===========================================================================

SearcherKorobov::SearcherKorobov (const FigureOfMerit* merit) : Searcher (merit)
{
}


//===========================================================================

SearcherKorobov::~SearcherKorobov ()
{
}


//===========================================================================

double SearcherKorobov::exhaustPrimePower2 (int s, long n)
{
   std::vector<long> z(s + 1, 1L);
   double err;
   double bestVal = DBL_MAX;
   long bestA = -1;

   for (long i = 1; i < n; i += 2) {
      calcAs (n, s, i, z);
      err = m_merit->compute (z, n, bestVal);
      if (err < bestVal) {
         bestVal = err;
         bestA = i;
      }
   }

   calcAs (n, s, bestA, m_bestAs);
   m_bestVal = bestVal;
   m_bestA = bestA;

   return bestVal;
}


//===========================================================================

double SearcherKorobov::exhaust (int s, long n, bool relPrime)
{
   if (relPrime && isPower2(n))
      return exhaustPrimePower2 (s, n);

   std::vector<long> z(s + 1, 1L);
   double err;
   double bestVal = DBL_MAX;
   long bestA = -1;

   for (int i = 1; i < n; i++) {
      // Consider only values i relatively prime to n
      if (relPrime && (gcd (n, i) != 1))
         continue;
      calcAs (n, s, i, z);
      err = m_merit->compute (z, n, bestVal);
      if (err < bestVal) {
         bestVal = err;
         bestA = i;
      }
   }

   calcAs (n, s, bestA, m_bestAs);
   m_bestVal = bestVal;
   m_bestA = bestA;
   return bestVal;
}


//===========================================================================

double SearcherKorobov::random (int s, long n, int k, bool relPrime)
{
   if (k >= n)
      return exhaust (s, n, relPrime);

   bool power2F = isPower2(n);
   const int t = (int) (Lg (n) + 0.5); // n = 2^t
   const int nm1 = n - 1;
   double err;
   double bestVal = DBL_MAX;
   long bestA = -1;
   std::vector<long> z(s + 1, 1L);
   long a;
   int i = 0;

   while (i < k) {
      if (power2F) {
         a = RandBits (t);
         //  a = unif01_StripB (m_gen, 0, t);
         if (relPrime || (a == 0))
            a |= 1;
      } else {
         do {
            a = RandInt (1, nm1);
            // a = 1 + unif01_StripL (m_gen, 0, nm1);
         } while (relPrime && (gcd (n, a) != 1));
      }

      i++;
      calcAs (n, s, a, z);
      err = m_merit->compute (z, n, bestVal);
      if (err < bestVal) {
         bestVal = err;
         bestA = a;
      }
   }

   m_bestVal = bestVal;
   m_bestA = bestA;
   calcAs (n, s, bestA, m_bestAs);
   return bestVal;
}


//===========================================================================

long SearcherKorobov::getBestA ()
{
   return m_bestA;
}

}
