#include "latmrg/PalphaProduct.h"
#include "latcommon/Num.h"

extern "C" {
#include "num.h"
}

#include <cmath>
#include <sstream>
#include <stdexcept>

#define JOE


using namespace std;
using namespace LatCommon;

namespace LatMRG
{

//===========================================================================

void PalphaProduct::init (int alpha)
{
   m_alpha = alpha;
   m_G = new double[m_dim + 1];
   setG (alpha, m_Gamma, m_dim);
   setBernoulli (alpha, m_n);
}


//===========================================================================

PalphaProduct::PalphaProduct (int alpha, long n, bool prime, const ProductWeights& weights, int d):
      Discrepancy(n, prime, weights, d)
{
   init (alpha);
}


//===========================================================================

PalphaProduct::PalphaProduct (int alpha, long n, bool prime, double beta[], int d):
      Discrepancy(n, prime, beta, d)
{
   init (alpha);
}


//===========================================================================

void PalphaProduct::setBernoulli (int alpha, long n)
{
   const double nd = n;
   double u;

   for (long j = 0; j < n; j++) {
      u = j / nd;
      m_C[j] = BernoulliPoly (alpha, u);
   }
}


//===========================================================================

PalphaProduct::~PalphaProduct()
{
   delete[] m_G;
}


//===========================================================================

const string PalphaProduct::toString () const
{
   string S(Discrepancy::toString ());
   string::size_type i = S.find("gamma");
   S.replace(i, 5, "beta ");
   ostringstream os;
   os << " alpha = " << m_alpha << endl;
   return "PalphaProduct:\n" + os.str () + S;
}


//===========================================================================

void PalphaProduct::setG (int alpha, double beta[], int d)
{
   for (int i = 1; i <= d; i++) {
#ifdef JOE
      m_G[i] = beta[i] * pow(2.0 * num_Pi, (double)alpha) / Factorial(alpha);
#else
      m_G[i] = pow(2.0 * num_Pi * beta[i], (double)alpha) / Factorial(alpha);
#endif
      if ((alpha / 2) & 1)
         m_G[i] = -m_G[i];
   }
   m_G[0] = 0.0;
}


//===========================================================================

double PalphaProduct::compute (const std::vector< long > & z, int n, double abortThreshold) const
{
   // FIXME
   if (n != m_n)
      throw std::invalid_argument("PalphaProduct: compute called with value of n different "
            "from that used during object instantiation");

   // z est le vecteur générateur du réseau de rang 1
   int d = z.size() - 1;
   double sum = 0;
   double prod;

   for (long k = 0; k < m_n; k++) {
      prod = 1.0;
      for (int j = 1; j <= d; j++) {
         long s = (k * z[j]) % m_n;
         prod *= 1.0 - m_G[j] * m_C[s];
      }
      sum += prod;
   }

   return (sum / m_n - 1.0);
}


//===========================================================================

}
