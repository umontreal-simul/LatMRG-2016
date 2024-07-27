#include "latmrg/OrdepBound.h"
#include "latcommon/Util.h"
#include "latcommon/Num.h"

#include <stdexcept>

using namespace std;
using namespace LatCommon;


namespace LatMRG
{

//===========================================================================

void OrdepBound::init()
{
   CreateMatr (m_Sigma, m_dim);
}


//===========================================================================

OrdepBound::OrdepBound (long n, bool prime, double gamma[], int d):
      BoundJSstar::BoundJSstar(n, prime, gamma, d)
{
   init();
}


//===========================================================================

OrdepBound::OrdepBound (long n, bool prime, int d, double gamma[], int dimg):
      BoundJSstar::BoundJSstar(n, prime, gamma, d)
{
   init();
   m_dimGam = dimg;
}


//===========================================================================

OrdepBound::OrdepBound (long n, long r, bool prime, double gamma[], int d):
      BoundJSstar::BoundJSstar(n, r, prime, gamma, d)
{
   init();
}


//===========================================================================

OrdepBound::~OrdepBound()
{
   DeleteMatr(m_Sigma, m_dim);
}


//===========================================================================

const string OrdepBound::toString () const
{
   return "OrdepBound:\n" + Discrepancy::toString ();
}


//===========================================================================

double OrdepBound::compute (const std::vector<long> & z, int n, double abortThreshold) const
{
   // FIXME
   if (n != m_n)
      throw std::invalid_argument("OrdepBound:  compute called with value of n different "
            "from that used during object instantiation");

   double sum = 0.0;
   int d = z.size() - 1;

   for (long k = 0; k < m_n; k++) {
     computeSigma(z, m_C, k);
//    for (int j = 1; j <= d; j++) {
      // m_Gamma[j] = 0 for j > m_dimGam
      for (int j = 1; j <= m_dimGam; j++)
         sum += m_Gamma[j] * m_Sigma[d][j];
   }

   return sum / m_n;
}


//===========================================================================

}
