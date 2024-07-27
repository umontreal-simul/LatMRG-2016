#include "latcommon/Util.h"
#include "latmrg/Searcher.h"
#include "latcommon/IntFactor.h"

#include <cmath>
#include <cstdlib>
#include <cfloat>
#include <iostream>
#include <iomanip>


using namespace std;
using namespace LatCommon;


//===========================================================================

namespace LatMRG
{

//===========================================================================

void Searcher::print (std::vector<long>& y, int s, double val)
{
   // For debugging: print vector y and associated merit value
   cout << "  a = [ ";
   for (int j = 1; j <= s; ++j) {
      cout << y[j] << " ";
   }
   cout << "],    val = " << val << endl;
}


//===========================================================================

void Searcher::incr (std::vector<long>& y, long n, int s)
{
   // increment y[i]; if y[i] == n, set y[i] = 1 and increment y[i-1];
   // recursively. Thus we consider all y[i], each component taking all
   // values in [1, 2, ... , n-1]
   // If y[2] >= n, stop. The exhaustive search of all y[j] is done.

   for (int i = s; i > 1; --i) {
      ++y[i];
      if (y[i] < n)
         return ;
      if (i > 2)
         y[i] = 1;
   }
}


//===========================================================================

void Searcher::incrPrime (std::vector<long>& y, long n, int s, bool power2F)
{
   // increment y[i] as in incr; but keeps only y[i] relatively prime to n.

   for (int i = s; i > 1; --i) {
      ++y[i];
      if (power2F) {
         ++y[i];
      } else {
         while (gcd (n, y[i]) != 1)
            ++y[i];
      }
      if (y[i] < n)
         return ;
      if (i > 2)
         y[i] = 1;
   }
}


//===========================================================================

double Searcher::exhaust (int s, long n, bool relPrime)
{
   // If relPrime is true, only the multipliers y_j relatively prime to n
   // are considered; otherwise, all multipliers are.

   bool power2F = isPower2(n);
   double err;
   std::vector<long> y(s + 1, 1L);
   m_bestVal = DBL_MAX;
   m_bestAs.assign(s + 1, 1L);

   while (y[2] < n) {
      err = m_merit->compute (y, n, m_bestVal);
      if (err < m_bestVal) {
         m_bestVal = err;
         // keep the best y in m_bestAs
         for (int j = 2; j <= s; j++)
            m_bestAs[j] = y[j];
      }
      if (relPrime)
         incrPrime (y, n, s, power2F);
      else
         incr (y, n, s);
   }

   return m_bestVal;
}


//===========================================================================

double Searcher::randomPower2 (int s, long n, int k, bool relPrime)
{
   // Random search when n is a power of 2.
   // If relPrime is true, only multipliers y_j relatively prime to n
   // are considered; otherwise, all multipliers are.

   const int t = (int) (Lg(n) + 0.5);    // n = 2^t
   double err;
   int j;
   m_bestVal = DBL_MAX;
   std::vector<long> y(s + 1, 1L);
   m_bestAs.assign(s + 1, 1L);
   y[1] = m_bestAs[1] = 1;

   for (int i = 0; i < k; i++) {
      for (j = 2; j <= s; j++) {
         //         y[j] = unif01_StripB (m_gen, 0, t);
         y[j] = RandBits (t);
         if (relPrime || (y[j] == 0))
            y[j] |= 1;
      }
      err = m_merit->compute (y, n, m_bestVal);
      if (err < m_bestVal) {
         m_bestVal = err;
         for (j = 2; j <= s; j++)
            m_bestAs[j] = y[j];
      }
   }

   return m_bestVal;
}


//===========================================================================

double Searcher::random (int s, long n, int k, bool relPrime)
{
   // If relPrime is true, only multipliers y_j relatively prime to n
   // are considered; otherwise, all multipliers are.

   if (isPower2(n))
      return randomPower2 (s, n, k, relPrime);

   const long nm1 = n - 1;
   double err;
   int j;
   m_bestVal = DBL_MAX;
   std::vector<long> y(s + 1, 1L);
   m_bestAs.assign(s + 1, 1L);
   y[1] = m_bestAs[1] = 1;

   for (int i = 0; i < k; i++) {
      for (j = 2; j <= s; j++) {
         do {
            // y[j] = 1 + unif01_StripL (m_gen, 0, nm1);
            y[j] = RandInt (1, nm1);
         } while (relPrime && (gcd (n, y[j]) != 1L));
      }
      err = m_merit->compute (y, n, m_bestVal);
      if (err < m_bestVal) {
         m_bestVal = err;
         for (j = 2; j <= s; j++)
            m_bestAs[j] = y[j];
      }
   }

   return m_bestVal;
}


//===========================================================================

Searcher::Searcher (const FigureOfMerit* merit)
{
   m_merit = merit;
}


//===========================================================================

Searcher::~Searcher ()
{
}

//===========================================================================

double Searcher::exhaust (int s, long n)
{
   return exhaust (s, n, false);
}


//===========================================================================

double Searcher::exhaustPrime (int s, long n)
{
   // When n is prime, all a_j are implicitly relatively prime with n,
   // so there is no need to constrain them.
   return exhaust (s, n, !isPrime(n));
}


//===========================================================================

double Searcher::random (int s, long n, int k)
{
   return random (s, n, k, false);
}


//===========================================================================

double Searcher::randomPrime (int s, long n, int k)
{
   // When n is prime, all a_j are implicitly relatively prime with n,
   // so there is no need to constrain them.
   return random (s, n, k, !isPrime(n));
}


//===========================================================================

double Searcher::getBestVal ()
{
   return m_bestVal;
}


//===========================================================================

const std::vector<long>& Searcher::getBestAs ()
{
   return m_bestAs;
}

//===========================================================================

void Searcher::initRand (unsigned long seed)
{
   SetSeed (seed);
}

//===========================================================================

bool Searcher::isPrime (long n)
{
   MScal big_n;
   conv (big_n, n);
   PrimeType p = IntFactor::isPrime (big_n, 50L);
   return p == PRIME || p == PROB_PRIME;
}

//===========================================================================

bool Searcher::isPower2 (long n)
{
   return (n & (n - 1)) == 0;
}

//===========================================================================
}
