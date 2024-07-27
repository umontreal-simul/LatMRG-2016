#include "latmrg/Rank1Lattice.h"
#include "latcommon/Util.h"
#include <cassert> 

#ifdef WITH_NTL
#else
using namespace boost::numeric::ublas;
#endif

using namespace std;
using namespace LatCommon;


namespace LatMRG
{

Rank1Lattice::Rank1Lattice (const MScal & n, const MVect & a, int maxDim,
                           NormType norm):
                  LatMRG::IntLattice::IntLattice (n, maxDim, maxDim, norm)
{
   m_a.resize(1 + maxDim);
   m_a = a;
   init();
}


//=========================================================================

Rank1Lattice::~Rank1Lattice()
{
   m_a.clear ();
}


//=========================================================================

void Rank1Lattice::init()
{
   IntLattice::init();   
   for (int r = 2; r <= getMaxDim(); r++)
      m_lgVolDual2[r] = m_lgVolDual2[r - 1];
}


//=========================================================================

Rank1Lattice & Rank1Lattice::operator= (const Rank1Lattice & lat)
{
   if (this == &lat)
      return * this;
   copy (lat);
   init ();
   int maxDim = lat.getMaxDim ();
   m_a.resize (1 + maxDim);
   m_a = lat.m_a;
   return *this;
}


//=========================================================================

Rank1Lattice::Rank1Lattice (const Rank1Lattice & lat):
      IntLattice::IntLattice (lat.m_m, lat.getOrder (),
                              lat.getMaxDim (), lat.getNorm ())
{
   // MyExit (1, "Rank1Lattice:: constructeur n'est pas terminé " );
   init ();
   int maxDim = lat.getMaxDim ();
   m_a.resize (1 + maxDim);
   m_a = lat.m_a;
}


//===========================================================================

std::string Rank1Lattice::toStringCoef ()const
{
   return toString (m_a, 1, getMaxDim ());
}


//=========================================================================

void Rank1Lattice::incDim ()
{
   // kill();
   buildBasis (1 + getDim ());
   m_v.setNegativeNorm (true);
   m_w.setNegativeNorm (true);
}


//===========================================================================

void Rank1Lattice::buildBasis (int d)
{
   // assert(d <= getMaxDim());
   setDim (d);

   // conv(m_v[1][1], 1);

   for (int j = 1; j <= d; j++) {
      m_v (1, j) = m_a[j];
   }

   for (int i = 2; i <= d; i++) {
      for (int j = 1; j <= d; j++) {
         if (i == j) {
            m_v (i, j) = m_m;
         } else {
            m_v (i, j) = 0;
         }
      }
   }

   // if a[1] != 1, the basis must be triangularized
   if (m_v (1, 1) != 1) {
      Triangularization < Base > (m_v, m_w, d, d, m_m);
      dualize ();
   }
   CalcDual < Base > (m_v, m_w, d, m_m);
   m_v.setNegativeNorm (true);
   m_w.setNegativeNorm (true);
}

}
