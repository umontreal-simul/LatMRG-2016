#include "latcommon/Util.h"
#include "latmrg/KorobovLattice.h"
#include <cassert>

#ifdef WITH_NTL
using namespace NTL;
#else
using namespace boost::numeric::ublas;
#endif

using namespace std;
using namespace LatCommon;


namespace LatMRG
{

KorobovLattice::KorobovLattice (const MScal & n, const MScal & a, int maxDim,
                                NormType norm) :
      LatMRG::IntLattice::IntLattice(n, 1, maxDim, norm)
{
   m_a = a;
   m_shift = 0;
   init();
}


//===========================================================================

KorobovLattice::KorobovLattice (const MScal & n, const MScal & a, int maxDim,
                                int t, NormType norm) :
      IntLattice::IntLattice(n, 1, maxDim, norm)
{
   m_a = a;
   m_shift = t;
   init();
}


//===========================================================================

KorobovLattice::~KorobovLattice()
{
}  


//=========================================================================

KorobovLattice::KorobovLattice (const KorobovLattice & lat):
            IntLattice::IntLattice(lat.m_m, 1, lat.getMaxDim (), lat.getNorm ())
{
   m_a = lat.m_a;
   m_shift = lat.m_shift;
}


//=========================================================================

KorobovLattice & KorobovLattice::operator= (const KorobovLattice & lat)
{
   if (this == &lat)
      return *this;
   copy (lat);
   m_a = lat.m_a;
   m_shift = lat.m_shift;
   return *this;
}


//===========================================================================

void KorobovLattice::init()
{
   IntLattice::init();   
//   double temp;
//   conv (temp, m_m);
//   m_lgVolDual2[1] = 2.0 * Lg(temp);

   for (int r = m_order + 1; r <= getMaxDim(); r++)
      m_lgVolDual2[r] = m_lgVolDual2[r - 1];
}


//===========================================================================

std::string KorobovLattice::toStringCoef () const
{
   std::ostringstream out;
   out << m_a;
   return out.str ();
}


//===========================================================================

void KorobovLattice::buildBasis(int d)
{
   assert(d <= getMaxDim());
   setDim(d);
   MScal tmp;
   conv(tmp, 1);

   for (int i = 0; i < m_shift; i++) {
      tmp *= m_a;
      tmp %= m_m;
   }

   for (int j = 1; j <= d; j++) {
      // V[1][j] = tmp % m;
      m_v (1, j) = tmp;
      tmp *= m_a;
      tmp %= m_m;
   }

   for (int i = 2; i <= d; i++) {
      for (int j = 1; j <= d; j++) {
         if (i == j)
            m_v (i, j) = m_m;
         else
            m_v (i, j) = 0;
      }
   }

   // Build dual basis
   CalcDual<Base>(m_v, m_w, d, m_m);
}


//===========================================================================

void KorobovLattice::incDimSlow()
{
   // Temporaire: très lent. Reprogrammer.
   int d = getDim();
   buildBasis(d + 1);
   m_v.setNegativeNorm(true);
   m_w.setNegativeNorm(true);
}


//===========================================================================

void KorobovLattice::incDim()
{
   const int dim = 1 + getDim();
   m_v.setDim(dim);
   m_w.setDim(dim);

   for (int i = 1; i < dim; i++) {
      conv (m_t2, m_v (i, dim - 1));
      m_t1 = m_t2 * m_a;
      Modulo (m_t1, m_m, m_t1);
      m_v (i, dim) = m_vSI (0, i) = m_t1;
   }

   matrix_row<Base> row1(m_v, dim - 1);
   SetZero (row1, dim - 1);

   for (int i = 1; i < dim; i++)
      m_v (dim, i) = 0;
   m_v (dim, dim) = m_m;

   for (int i = 1; i < dim; i++)
      m_w (i, dim) = 0;
   m_w (dim, dim) = 1;

   for (int j = 1; j < dim; j++) {
      clear (m_t1);
      for (int i = 1; i < dim; i++) {
         m_t2 = m_w (i, j);
         m_t2 *= m_vSI (0, i);
         m_t1 -= m_t2;
      }
      Quotient (m_t1, m_m, m_t3);
      m_w (dim, j) = m_t3;
   }

   m_v.setNegativeNorm(true);
   m_w.setNegativeNorm(true);
}


//===========================================================================

}
