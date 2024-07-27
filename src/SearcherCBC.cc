#include "latcommon/Util.h"
#include "latmrg/SearcherCBC.h"

#include <cfloat>


using namespace LatCommon;


//===========================================================================

namespace LatMRG
{

//===========================================================================

double SearcherCBC::exhaust (int s, long n, bool relPrime)
{
   /* Exhaustive search for each dimension j <= s. If relPrime is true,
      only values i relatively prime to n are examined */
   bool power2F = isPower2(n);
   int pos = -1;
   double err;
   double best = -1;
   int i;
   m_bestVals.assign(s + 1, -1);
   std::vector<long> y(1, 0L);  // first element is unused
   m_bestAs.assign(1, 0L);
   // increase the size of y by one element
   y.push_back(1L);
   m_bestAs.push_back(1L);

   for (int j = 2; j <= s; j++) {
      // increase the size of y by one element
      y.push_back(1L);
      m_bestAs.push_back(1L);

      best = DBL_MAX;
      pos = -1;
      for (i = 1; i < n; i++) {
         if (relPrime) {
            // Test only values i relatively prime to n
            if (power2F) {
               if ((i & 1) == 0)
                  continue;
            } else {
               if (gcd (n, i) != 1)
                  continue;
            }
         }
         y[j] = i;
         err = m_merit->compute (y, n, best);
         if (err < best) {
            best = err;
            pos = i;
         }
      }
      y[j] = m_bestAs[j] = pos;
      m_bestVals[j] = best;
   }

   m_bestVal = m_bestVals[s];
   return m_bestVal;
}


//===========================================================================

double SearcherCBC::random (int s, long n, int k, bool relPrime)
{
   /* Random search for each dimension j <= s. k values are examined for each
      dimension j. If relPrime is true, only values i relatively prime to n
      are examined */
   if (k >= n)
      return exhaust (s, n, relPrime);

   bool power2F = isPower2(n);
   const int t = (int) (Lg(n) + 0.5);    // n = 2^t
   const int nm1 = n - 1;
   int pos = -1;
   double err;
   double best = 0;
   int i;
   long x;
   m_bestVals.assign(s + 1, -1);
   std::vector<long> y(1, 0L);  // first element is unused
   m_bestAs.assign(1, 0L);
   // increase the size of y by one element
   y.push_back(1L);
   m_bestAs.push_back(1L);

   for (int j = 2; j <= s; j++) {
      // increase the size of y by one element
      y.push_back(1L);
      m_bestAs.push_back(1L);

      best = DBL_MAX;
      i = 0;
      while (i < k) {
         if (power2F) {
            x = RandBits (t);
            //   x = unif01_StripB (m_gen, 0, t);
            if (relPrime || (x == 0))
               x |= 1;
         } else {
            do {
               x = RandInt (1, nm1);
               // x = 1 + unif01_StripL (m_gen, 0, nm1);
            } while (relPrime && (gcd (n, x) != 1));
         }
         y[j] = x;
         err = m_merit->compute (y, n, best);
         if (err < best) {
            best = err;
            pos = y[j];
         }
         i++;
      }
      y[j] = m_bestAs[j] = pos;
      m_bestVals[j] = best;
   }

   m_bestVal = m_bestVals[s];
   return m_bestVal;
}


//===========================================================================

SearcherCBC::SearcherCBC (const FigureOfMerit* merit) : Searcher (merit)
{
}


//===========================================================================

SearcherCBC::~SearcherCBC ()
{
}


//===========================================================================

const std::vector<double>& SearcherCBC::getBestVals ()
{
   return m_bestVals;
}


//===========================================================================

}
