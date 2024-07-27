#include "latmrg/BoundJSstar.h"
#include "latcommon/Num.h"

using namespace std;
using namespace LatCommon;


namespace LatMRG
{

//===========================================================================

void BoundJSstar::init (double gamma[], int d)
{
   m_N = m_n * m_r;
   m_Beta = new double[m_dim + 1];
   setBeta (gamma, d);
   setFourier1 (m_C, m_n, m_N);
}


//===========================================================================

BoundJSstar::BoundJSstar (long n, bool prime, double gamma[], int d):
      Discrepancy::Discrepancy(n, prime, gamma, d)
{
   init (gamma, d);
}


//===========================================================================

BoundJSstar::BoundJSstar (long n, long r, bool prime, double gamma[], int d):
      Discrepancy::Discrepancy(n, r, prime, gamma, d)
{
   init (gamma, d);
}


//===========================================================================

BoundJSstar::~BoundJSstar()
{
   delete[] m_Beta;
}


//===========================================================================

const string BoundJSstar::toString () const
{
   return "BoundJSstar:\n" + Discrepancy::toString ();
}


//===========================================================================

void BoundJSstar::setBeta (double gam[], int d)
{
   for (int i = 1; i <= d; i++) {
      m_Beta[i] = 1.0 + gam[i];
   }
   m_Beta[0] = 0.0;
}


//===========================================================================
#if 1

double BoundJSstar::compute (long z[], int d)
{
   double bound = 0;
   double prod;
   long s;
   int j;

   for (long k = 0; k < m_n; k++) {
      prod = 1.0;
      for (j = 1; j <= d; j++) {
         s = (k * z[j]) % m_n;
         prod *= (m_Beta[j] + m_Gamma[j] * m_C[s]);
      }
      bound += prod;
   }

   prod = 1.0;
   for (j = 1; j <= d; j++)
      prod *= m_Beta[j];

   return bound / m_n - prod;
}

#else

double BoundJSstar::compute (long z[], int d)
{
   // This method is 3 times faster than the method compute above.
   // In principle, this method should give the same result as above, but
   // there are small differences sometimes: the reason
   // is that small numerical inaccuracies cause a different z[j] to be the
   // best sometimes, even though they gives the same result mathematically.
   // In dimension j, for example, these two give the same bound (within
   // epsilon), but for dimensions greater than j, the best CBC component
   // may not be the same in the two cases, and will give different bounds.
   double bound = 0;
   double prod;
   long s;
   int j;

   for (long k = 1; k <= (m_n - 1) / 2; k++) {
      prod = 1.0;
      for (j = 1; j <= d; j++) {
         s = (k * z[j]) % m_n;
         prod *= (m_Beta[j] + m_Gamma[j] * m_C[s]);
      }
      bound += prod;
   }
   bound *= 2.0;

   prod = 1.0;
   for (j = 1; j <= d; j++)
      prod *= (m_Beta[j] + m_Gamma[j] * m_C[0]);
   bound += prod;

   if ((m_n & 1) == 0) {
      prod = 1.0;
      for (j = 1; j <= d; j++) {
         s = ((m_n/2) * z[j]) % m_n;
         prod *= (m_Beta[j] + m_Gamma[j] * m_C[s]);
      }
      bound += prod;
   }

   prod = 1.0;
   for (j = 1; j <= d; j++)
      prod *= m_Beta[j];

   return bound / m_n - prod;
}
#endif
//===========================================================================

}
