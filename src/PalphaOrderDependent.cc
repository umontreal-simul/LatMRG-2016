#include "latmrg/PalphaOrderDependent.h"
#include "latcommon/Util.h"

#include <sstream>
#include <stdexcept>

using namespace std;
using namespace LatCommon;


namespace LatMRG
{

//===========================================================================

void PalphaOrderDependent::init(int alpha)
{
   m_alpha = alpha;
   setFourier(m_C, m_n, m_alpha);
   CreateMatr (m_Sigma, m_dim);
}


//===========================================================================

PalphaOrderDependent::PalphaOrderDependent (int alpha, long n, bool prime, const OrderDependentWeights& weights, int d):
      Discrepancy::Discrepancy(n, prime, weights, d)
{
   init(alpha);
}


//===========================================================================

PalphaOrderDependent::PalphaOrderDependent (int alpha, long n, bool prime, double gamma[], int d):
      Discrepancy::Discrepancy(n, prime, gamma, d)
{
   init(alpha);
}


//===========================================================================

PalphaOrderDependent::PalphaOrderDependent (int alpha, long n, bool prime, int d,
                          double gamma[], int dimg):
   Discrepancy::Discrepancy(n, prime, d, gamma, dimg)
{
   init(alpha);
}


//===========================================================================

PalphaOrderDependent::~PalphaOrderDependent()
{
   DeleteMatr (m_Sigma, m_dim);
}


//===========================================================================

const string PalphaOrderDependent::toString () const
{
   ostringstream os;
   os << " alpha = " << m_alpha << endl;
   return "PalphaOrderDependent:\n" + os.str () + Discrepancy::toString ();
}


//===========================================================================

double PalphaOrderDependent::compute (const std::vector< long > & z, int n, double abortThreshold) const
{
   // FIXME
   if (n != m_n)
      throw std::invalid_argument("PalphaProduct: compute called with value of n different "
            "from that used during object instantiation");

   double sum = 0.0;
   int d = z.size() - 1;

   for (long k = 0; k < m_n; k++) {
      computeSigma(z, m_C, k);
//    for (j = 1; j <= d; j++) {
      // m_Gamma[j] = 0 for j > m_dimGam
      for (int j = 1; j <= m_dimGam; j++)
         sum += m_Gamma[j] * m_Sigma[d][j];
   }

   return sum / m_n;
}


//===========================================================================

}
