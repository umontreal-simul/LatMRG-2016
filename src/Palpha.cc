#include "latcommon/Num.h"
#include "latmrg/Palpha.h"

extern "C" {
#include "num.h"
}

//!#define DEBUG

#include <cmath>
#include <sstream>
#include <stdexcept>

#ifdef DEBUG
#include <iostream>
#endif


using namespace std;
using namespace LatCommon;

namespace LatMRG
{

//===========================================================================

Palpha::Palpha (int alpha)
{
   if (alpha % 2 != 0)
      throw std::invalid_argument("Palpha: alpha must be an even integer");

   m_alpha = alpha;
}


//===========================================================================

const string Palpha::toString () const
{
   ostringstream os;
   os << "P" << m_alpha << " (general implementation)" << endl;
   os << "  (P_alpha with alpha=" << m_alpha << ")";
   return os.str ();
}


//===========================================================================

double Palpha::compute (const std::vector<long>& a, int n, const std::vector<int>& projection) const
{
   if (projection.size() == 0)
      return 1.0;

#ifdef DEBUG
   std::cout << "computing general Palpha in dimension " << projection.size();
   std::cout << " with generating vector [";
   for (std::vector<long>::const_iterator it = a.begin() + 1; it != a.end(); ++it)
      std::cout << " " << *it;
   std::cout << " ]" << std::endl;
#endif

   m_bern.init(m_alpha, n);

   double sum = 0.0;
   for (int i = 0; i < n; i++) {
      double prod = 1.0;
      for (std::vector<int>::const_iterator it = projection.begin(); it != projection.end(); ++it)
         prod *= m_bern[i * a[*it]];
      sum += prod;
   }
   sum /= n;

#ifdef DEBUG
   std::cout << "  Palpha value: " << sum << std::endl;
#endif

   return sum;
}

//===========================================================================

void Palpha::BernCache::init (int alpha, int n) {
   // if nothing has changed, do nothing
   if (alpha == m_alpha && n == m_n)
      return;
   // update the scaled values of the Bernouilli polynomial
   m_alpha = alpha;
   m_n = n;
   m_b.clear();
   m_b.reserve(n);
   double scal = - (alpha % 4 ? -1 : 1) * pow(2.0 * num_Pi, (double)alpha) / Factorial(alpha);
   for (int i = 0; i < n; i++)
      m_b[i] = scal * BernoulliPoly (alpha, i / (double) n);
}

//===========================================================================

}
