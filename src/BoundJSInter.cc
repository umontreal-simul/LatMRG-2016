#include "latmrg/BoundJSInter.h"
#include "latcommon/Num.h"

#include <cassert>
#include <sstream>

using namespace std;
using namespace LatCommon;


namespace LatMRG
{

//===========================================================================

void BoundJSInter::init (long el, double gamma[], int d)
{
   m_el = el;
   setN();

   m_zTil = new long[d + 1];
   m_GammaTil = new double[d + 1];
   m_CTil = new double[m_n];
   m_Beta = new double[m_dim + 1];
   setBeta (gamma, d);

   int j;
   for (j = 1; j <= m_r; j++)
      m_GammaTil[j] = gamma[j]/el;

   for (j = m_r + 1; j <= d; j++)
      m_GammaTil[j] = gamma[j];

   setFourier1 (m_C, m_n, m_N);
   setFourier1 (m_CTil, m_n, m_NTil);
/*
   m_zTil[0] = 0;
   m_Beta[0] = 0;
   m_Gamma[0] = 0;
   m_GammaTil[0] = 0;
*/
}


//===========================================================================

BoundJSInter::BoundJSInter (long n, long r, long el, double gamma[], int d):
   Discrepancy::Discrepancy(n, r, true, gamma, d)
{
   assert (r <= d);
   init(el, gamma, d);
}


//===========================================================================

BoundJSInter::~BoundJSInter()
{
   delete[] m_Beta;
   delete[] m_CTil;
   delete[] m_GammaTil;
   delete[] m_zTil;
}


//===========================================================================

const string BoundJSInter::toString () const
{
   std::ostringstream os;
   os << " n     = " << m_n
      << "\n r     = " << m_r
      << "\n el    = " << m_el
      << "\n dim   = " << m_dim
      << "\n gamma = [ ";
   for (int i = 1; i < m_dim; i++)
      os << m_Gamma[i] << "  ";
   os << m_Gamma[m_dim] << " ]\n";
   return "BoundJSInter:\n" + os.str ();
}


//===========================================================================

void BoundJSInter::setBeta (double gam[], int d)
{
   for (int i = 1; i <= d; i++)
      m_Beta[i] = 1.0 + gam[i];
//   m_Beta[0] = 0.0;
}


//===========================================================================

void BoundJSInter::setN()
{
   m_N = m_n;
   for (int j = 0; j < m_r; j++)
      m_N *= m_el;
   m_NTil = m_N/m_el;
}


//===========================================================================

long BoundJSInter::getNumPoints()
{
   return m_N;
}


//===========================================================================

double BoundJSInter::compute (long z[], int d)
{
   double bound = 0;
   double prod;
   long s;
   int j;

   for (j = 1; j <= m_r; j++)
      m_zTil[j] = z[j] * m_el;

   for (long k = 0; k < m_n; k++) {
      prod = 1.0;
      for (j = 1; j <= m_r; j++) {
         s = (k*m_zTil[j]) % m_n;
         prod *= (m_Beta[j] + m_GammaTil[j] * m_CTil[s]);
      }
      for (j = 1 + m_r; j <= d; j++) {
         s = (k*z[j]) % m_n;
         prod *= (m_Beta[j] + m_GammaTil[j] * m_C[s]);
      }
      bound += prod;
   }

   prod = 1.0;
   for (j = 1; j <= d; j++)
      prod *= m_Beta[j];

   return bound / m_n - prod;
}


//===========================================================================

}
