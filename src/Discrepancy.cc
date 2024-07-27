#include "latmrg/Discrepancy.h"
#include "latcommon/Util.h"
#include "latcommon/Num.h"
#include "NTL/ZZ.h"

extern "C" {
#include "num.h"
}

#include <cfloat>
#include <cmath>
#include <cassert>

using namespace std;
using namespace LatCommon;


namespace LatMRG
{

//===========================================================================

void Discrepancy::init (long n, bool prime, int dim, long r, double gamma[], int dimg)
{
   m_n = n;
   m_prime = prime;
   if ((n & (n - 1)) == 0) {   // true if n = 2^s
      m_power2 = true;
      if (n > 2)
         m_prime = false;
   } else
      m_power2 = false;
   m_dim = dim;

   m_N = -1;
   m_r = r;
   m_dimGam = dimg;
   m_Gamma = new double[m_dim + 1];
   setWeights (gamma, m_dim);
   m_C = new double[m_n];
   ytemp = new long[m_dim + 1];
}


//===========================================================================

Discrepancy::Discrepancy (long n, bool prime, double gamma[], int dim)
{
   init (n, prime, dim, 1, gamma, dim);
}

//===========================================================================

Discrepancy::Discrepancy (long n, bool prime, const ProductWeights & gamma, int dim)
{
   init (n, prime, dim, 1, 0, dim);
   setWeights (gamma);
}

//===========================================================================

Discrepancy::Discrepancy (long n, bool prime, const OrderDependentWeights & gamma, int dim)
{
   init (n, prime, dim, 1, 0, dim);
   setWeights (gamma);
}

//===========================================================================

Discrepancy::Discrepancy (long n, bool prime, int dim, double gamma[],
                          int dimg)
{
   init (n, prime, dim, 1, gamma, dimg);
}


//===========================================================================

Discrepancy::Discrepancy (long n, long r, bool prime, double gamma[], int dim)
{
   init (n, prime, dim, r, gamma, dim);
}

//===========================================================================

Discrepancy::Discrepancy (long n, long r, bool prime, const ProductWeights & gamma, int dim)
{
   init (n, prime, dim, r, 0, dim);
   setWeights (gamma);
}

//===========================================================================

Discrepancy::Discrepancy (long n, long r, bool prime, const OrderDependentWeights & gamma, int dim)
{
   init (n, prime, dim, r, 0, dim);
   setWeights (gamma);
}

//===========================================================================

Discrepancy::~Discrepancy ()
{
   delete[] ytemp;
   delete[] m_C;
   delete[] m_Gamma;
}


//===========================================================================

long Discrepancy::getNumPoints ()
{
   return m_n;
}


//===========================================================================

int Discrepancy::getDimension ()
{
   return m_dim;
}


//===========================================================================

bool Discrepancy::isPower2N ()
{
   return m_power2;
}

//===========================================================================

bool Discrepancy::isPrimeN ()
{
   return m_prime;
}


//===========================================================================

void Discrepancy::setWeights (double gam[], int d)
{
   if (gam != 0)
      for (int i = 0; i <= d; i++)
         m_Gamma[i] = gam[i];
   else
      for (int i = 0; i <= d; i++)
         m_Gamma[i] = 1;
}


//===========================================================================

void Discrepancy::setWeights (const ProductWeights & gamma)
{
   m_Gamma[0] = gamma.getWeight(LatCommon::Coordinates());
   for (int i = 1; i <= m_dim; i++) {
      LatCommon::Coordinates proj;
      proj.insert(i);
      m_Gamma[i] = gamma.getWeight(proj);
   }
}


//===========================================================================

void Discrepancy::setWeights (const OrderDependentWeights & gamma)
{
   LatCommon::Coordinates proj;
   m_Gamma[0] = gamma.getWeight(proj);
   for (int i = 1; i <= m_dim; i++) {
      proj.insert(i);
      m_Gamma[i] = gamma.getWeight(proj);
   }
}


//===========================================================================

const string Discrepancy::toString () const
{
   //string S(FigureOfMerit::toString ());
   std::ostringstream os;
   os << " n       = " << m_n
      << "\n prime   = " << boolalpha << m_prime;
   if (m_power2)
      os << "\n power of 2 = " << boolalpha << m_power2;
   os << "\n dim     = " << m_dim;
   os << "\n r       = " << m_r
      << "\n gamma   = [ ";
   for (int i = 1; i < m_dim; i++)
      os << m_Gamma[i] << "  ";
   os << m_Gamma[m_dim] << " ]\n";
   return os.str ();
}


//===========================================================================

void Discrepancy::setFourier1 (double *C, long n, long N)
{
   const double nd = n;
   double x;

   for (long k = 0; k < n; k++) {
      x = k / nd;
      C[k] = FourierE1 (x, N);
   }
}


//===========================================================================

void Discrepancy::setFourier (double *C, long n, int alpha)
{
   assert ((alpha & 1) == 0);
   const double nd = n;
   int sign = 1;
   if ((alpha/2 & 1) == 0)
      sign = -1;
   double D = pow(2.0 * num_Pi, alpha);
   double fact = sign * D / Factorial(alpha);
   double x;

   for (long k = 0; k < n; k++) {
      x = k / nd;
      C[k] = fact * BernoulliPoly (alpha, x);
   }
}


//===========================================================================

void Discrepancy::computeSigma(const std::vector<long> & z, double* C, long k) const
{
   int d = z.size() - 1;
   long s = 0;
   m_Sigma[1][1] = C[k];

   for (int m = 2; m <= d; m++) {
      double sum = 0.0;
      double prod = 1.0;
      for (int j = 1; j <= m; j++) {
         s = (k * z[j]) % m_n;
         sum += C[s];
         prod *= C[s];
      }
      m_Sigma[m][1] = sum;

//    for (j = 2; j < m; j++)
      for (int j = 2; j <= m_dimGam; j++)  // m_Gamma[j] = 0 for j > m_dimGam
         m_Sigma[m][j] = m_Sigma[m - 1][j] + C[s] * m_Sigma[m - 1][j - 1];
      m_Sigma[m][m] = prod;
   }
}

//===========================================================================

}
