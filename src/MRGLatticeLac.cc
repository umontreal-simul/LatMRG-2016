#include "latmrg/MRGLatticeLac.h"
#include "latmrg/PolyPE.h"
#include "latcommon/Util.h"

#include <sstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;
using namespace NTL;
using namespace LatCommon;

//============================================================================

namespace LatMRG
{

//===========================================================================

MRGLatticeLac::MRGLatticeLac (const MRGLatticeLac & lat): MRGLattice::
      MRGLattice (lat.m_m, lat.m_aCoef, lat.getMaxDim (), lat.getOrder (),
                  lat.m_latType, lat.getNorm ()), m_lac (lat.m_lac)
{
   m_lacunaryFlag = true;
   MyExit (1, " MRGLatticeLac::  copy constructeur n'est pas terminé   ");
}


//=========================================================================

MRGLatticeLac & MRGLatticeLac::operator= (const MRGLatticeLac & lat)
{
   if (this == &lat)
      return * this;
   MyExit (1, " MRGLatticeLac::operator= n'est pas terminé   ");
   copy (lat);
   return *this;
}


//===========================================================================

MRGLatticeLac::MRGLatticeLac (const MScal & m, const MVect & a, int maxDim,
                              int k, BVect & Lac, LatticeType lat,
                              NormType norm):
      MRGLattice::MRGLattice (m, a, maxDim, k, lat, norm),
      m_lac (Lac, maxDim)
{
   m_lacunaryFlag = true;
   m_sta.SetDims(1, 1);
   for (int i = 1; i <= m_order; i++)
      m_aCoef[i] = a[i];
}


//============================================================================

MRGLatticeLac::~MRGLatticeLac ()
{
}


//=========================================================================

BScal & MRGLatticeLac::getLac (int j)
{
   if (j <= m_lac.getSize () && j > 0)
      return m_lac.getLac (j);
   throw std::out_of_range ("MRGLatticeLac::getLac");
}


//===========================================================================

void MRGLatticeLac::setLac (const Lacunary & lac)
{
   m_lac = lac;
   m_lacunaryFlag = true;
}


//===========================================================================

void MRGLatticeLac::buildBasis (int d)
{
   int ord = m_order;

   initStates ();
   int IMax = m_lac.getSize ();

   MVect b;
   b.SetLength (m_order + 1);
   Invert (m_aCoef, b, m_order);

   // b is the characteristic polynomial
   PolyPE::setM (m_m);
   PolyPE::setF (b);
   PolyPE pol;

   // Construction d'un systeme generateur modulo m.
   for (int k = 1; k <= IMax; k++) {
      // pour chaque indice lacunaire
      conv (m_e, m_lac[k]);

      // x^m_e Mod f(x) Mod m
      pol.powerMod (m_e);
      pol.toVector (m_xi);

      for (int i = 1; i <= m_order; i++) {
           m_wSI[i][k] = m_xi[i - 1];
      }
   }

   /* On veut s'assurer que la base m_v soit triangulaire (pour satisfaire
      les conditions de l'article \cite{rLEC94e} [sec. 3, conditions sur
      V_i >= i]) et de plein rang (on remplace les lignes = 0 par lignes
      avec m sur la diagonale).
    */
   /* Il serait possible de réserver m_wSI, m_vSI avec seulement IMax
      lignes et colonnes quand order >> IMax. Mais il faudrait un nouveau
      Triangularization qui pourrait nécessiter un appel pour chaque ligne,
      mais qui sauverait beaucoup de mémoire. 
      Il n'est pas certain que cela en vaille la peine. */
   Triangularization <BMat> (m_wSI, m_vSI, ord, IMax, m_m);
   CalcDual <BMat> (m_vSI, m_wSI, IMax, m_m);

   // Construire la base de dimension 1
   m_v[1][1] = m_vSI[1][1];
   m_w[1][1] = m_wSI[1][1];
   setDim (1);

   m_v.setNegativeNorm (true);
   m_w.setNegativeNorm (true);
   for (int i = 2; i <= d; i++)
      incDimBasis (IMax);

   // for debugging
   // trace("ESPION_2", i);
}


//===========================================================================

void MRGLatticeLac::incDimBasis (int IMax)
{
   const int dim = getDim ();

   if (dim >= IMax) {
      MyExit (0,
    "Dimension of the basis is too big:\nDim > Number of lacunary indices.");
   }

   for (int i = 1; i <= dim; i++) {
      // v[i] -> VSI[0].
      for (int j = 1; j <= dim; j++)
         m_vSI[0][j] = m_v[i][j];
      clear (m_vSI[i][0]);

      for (int i1 = 1; i1 <= dim; i1++) {
         ProdScal (m_vSI[0], m_wSI[i1], dim, m_wSI[i1][0]);
         Quotient (m_wSI[i1][0], m_m, m_wSI[i1][0]);
         m_t1 = m_wSI[i1][0] * m_vSI[i1][dim + 1];
         m_vSI[i][0] += m_t1;
      }
      Modulo (m_vSI[i][0], m_m, m_vSI[i][0]);
      m_v[i][dim + 1] = m_vSI[i][0];
   }

   for (int j = 1; j <= dim; j++)
      m_v[dim + 1][j] = 0;
   m_v[dim + 1][dim + 1] = m_vSI[dim + 1][dim + 1];

   for (int i = 1; i <= dim; i++)
      m_w[i][dim + 1] = 0;

   for (int j = 1; j <= dim; j++) {
      clear (m_wSI[0][j]);
      for (int i = 1; i <= dim; i++) {
         m_t1 = m_w[i][j];
         m_t1 *= m_vSI[i][0];
         m_wSI[0][j] += m_t1;
      }
      if (m_wSI[0][j] != 0)
         m_wSI[0][j] = -m_wSI[0][j];
      Quotient (m_wSI[0][j], m_vSI[dim + 1][dim + 1], m_wSI[0][j]);
      m_w[dim + 1][j] = m_wSI[0][j];
   }

   Quotient (m_m, m_vSI[dim + 1][dim + 1], m_t1);
   m_w[dim + 1][dim + 1] = m_t1;

   setDim (dim + 1);
   m_v.setNegativeNorm (true);
   m_w.setNegativeNorm (true);
//    trace("ESPION_3", dim);
}


//===========================================================================

void MRGLatticeLac::initStates ()
/*
 * Initialise la matrice carrée Sta. La matrice Sta est d'ordre égal à l'ordre du
 * générateur. Elle contient un système de générateurs pour le groupe d'états
 * considérés.
 */
{
   clear (m_t2);

   if (m_latType == RECURRENT) {
      // check if a_k is relatively prime with m ==> m_t1 = 1
      m_t1 = GCD (m_aCoef[m_order], m_m);
      m_t1 = abs (m_t1);
      set9 (m_t2);
   }

   if (m_latType == FULL || m_latType == PRIMEPOWER || (m_t1 == m_t2)) {
      // m_sta is set to identity matrix
      calcLgVolDual2 (m_lgm2);

   } else {
      if (m_latType == ORBIT) {
         MyExit (1, "case ORBIT ne fonctionne pas");
      } else if (m_latType == RECURRENT) {
         MyExit (1, "case RECURRENT ne fonctionne pas");
      }
   }
}

//===========================================================================
}
