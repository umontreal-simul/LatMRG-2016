#include "latmrg/MRGLattice.h"
#include "latmrg/MRGComponent.h"
#include "latcommon/Util.h"
#include "latmrg/PolyPE.h"

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

/* Max order for lacunary case in this class; takes too much memory.
 For order > ORDERMAX, use subclass MRGLatticeLac instead */
#define ORDERMAX 100

//===========================================================================

void MRGLattice::trace (char *mess, int d)
{
   cout << "-------------------------------------------------------------------" << endl;
   cout << mess << endl;
   m_v.setNegativeNorm(true);
   m_w.setNegativeNorm(true);
   m_v.updateVecNorm ();
   m_w.updateVecNorm ();
   m_v.write();
   m_w.write();
   /*
   for (int i = 0; i <= d; i++)
      cout << " VSI " << i << "    " << m_vSI[i] << endl;
   cout << endl;
   for (int i = 0; i <= d; i++)
      cout << " WSI " << i << "    " << m_wSI[i] << endl;
      */
   checkDuality ();
   d = -1;  // compiler warning
}


//===========================================================================

MRGLattice::MRGLattice(const MRGLattice &lat):
  IntLattice::IntLattice (lat.m_m, lat.getOrder(),
                          lat.getMaxDim (), lat.getNorm ()), m_lac(lat.m_lac)
{
   m_lossRho = lat.m_lossRho;
   m_rho = lat.m_rho;
   m_latType = lat.m_latType;
   m_lacunaryFlag = lat.m_lacunaryFlag;

   m_ip = new bool[m_order + 1];
   m_xi.SetLength (m_order + 1);
   m_aCoef.SetLength (m_order + 1);
   m_sta.SetDims (m_order + 1, m_order + 1);

   int maxDim = getMaxDim();
   int rmax = max(m_order, maxDim);
   m_wSI.SetDims (rmax + 1, maxDim + 1);

   int i;
   for (i = 0; i <= m_order; i++)
      m_aCoef[i] = lat.m_aCoef[i];
/*
   for (i = 0; i <= m_order; i++)
      m_xi[i] = lat.m_xi[i];
   for (i = 0; i <= m_order; i++)
      m_ip[i] = lat.m_ip[i];

   int j;
   for (i = 0; i <= m_order; i++)
      for (j = 0; j <= m_order; j++)
         m_sta[i][j] = lat.m_sta[i][j];

   for (i = 0; i <= maxDim; i++)
      for (j = 0; j <= maxDim; j++)
         m_wSI[i][j] = lat.m_wSI[i][j];
*/
}


//=========================================================================

MRGLattice & MRGLattice::operator= (const MRGLattice & lat)
{
   if (this == &lat)
      return *this;
   MyExit (1, " MRGLattice::operator= n'est pas terminé   " );
   copy (lat);
   return *this;
}


//===========================================================================

MRGLattice::MRGLattice(const MScal & m, const MVect & a, int maxDim, int k,
                       LatticeType lat, NormType norm):
      IntLattice::IntLattice(m, k, maxDim, norm)
{
   m_latType = lat;
   m_lacunaryFlag = false;
   m_ip = new bool[1];
   init();

   for (int i = 1; i <= m_order; i++)
      m_aCoef[i] = a[i];
}


//===========================================================================

MRGLattice::MRGLattice(const MScal & m, const MVect & a, int maxDim, int k,
                       BVect & I, LatticeType lat, NormType norm):
      IntLattice::IntLattice (m, k, maxDim, norm), m_lac(I, maxDim), m_ip(0)
{
   m_latType = lat;
   m_lacunaryFlag = true;
   init();

   for (int i = 1; i <= m_order; i++)
      m_aCoef[i] = a[i];
}

//===========================================================================

void MRGLattice::init()
{
   kill();
   IntLattice::init();
   m_xi.SetLength(m_order + 1);
   m_aCoef.SetLength(m_order + 1);
   if (m_order > ORDERMAX) {
      m_ip = new bool[1];
      m_sta.SetDims(1, 1);
   } else {
      m_ip = new bool[m_order + 1];
      m_sta.SetDims(m_order + 1, m_order + 1);
   }
   int rmax = max(m_order, getMaxDim());
   m_wSI.SetDims(rmax + 1, getMaxDim() + 1);

   if (m_latType == ORBIT)
      initOrbit();
}

//============================================================================

MRGLattice::~MRGLattice ()
{
   kill();
}


//============================================================================

void MRGLattice::kill()
{
   IntLattice::kill();
   if (0 != m_ip)
      delete[] m_ip;
   m_ip = 0;
   m_xi.kill();
   m_aCoef.kill();
   m_sta.kill();
   m_wSI.kill();
}


//=========================================================================

BScal & MRGLattice::getLac (int j)
{
   if (isLacunary() && j <= m_lac.getSize() && j > 0)
      return m_lac.getLac(j);
   throw std::out_of_range("MRGLattice::getLac");
}


//===========================================================================

void MRGLattice::setLac(const Lacunary & lac)
{
   m_lac = lac;
   m_lacunaryFlag = true;
}


//===========================================================================

string MRGLattice::toStringCoef () const
{
   std::ostringstream out;
   out << "[ ";
   for (int i = 1; i <= m_order; i++)
      out << m_aCoef[i] << "  ";
   out << "]";
   return out.str ();
}


//===========================================================================

void MRGLattice::buildBasis (int d)
{
   if (m_lacunaryFlag) {
      buildLaBasis(d);
   } else {
      buildNaBasis(d);
   }

}


//===========================================================================

void MRGLattice::buildNaBasis (int d)
// La base est construite en dimension d.
{
 // trace( "=====================================AVANT buildNaBasis", -10);
   initStates();

   int dk = d;
   if (dk > m_order)
      dk = m_order;

   int i, j;
   for (i = 1; i <= dk; i++) {
      if (m_ip[i]) {
         for (j = 1; j <= dk; j++)
            m_v[i][j] = m_sta[i][j];

      } else {
         for (j = 1; j <= dk; j++) {
            if (i != j)
               m_v[i][j] = 0;
            else
               m_v[i][j] = m_m;
         }
      }
   }

   CalcDual<BMat>(m_v, m_w, dk, m_m);

   setDim(dk);
   if (d > m_order) {
      for (i = m_order + 1; i <= d; i++)
         incDimBasis ();
   }
 // trace( "=================================APRES buildNaBasis", -10);
}


//===========================================================================

void MRGLattice::incDim()
{
   if (m_lacunaryFlag) {
      incDimLaBasis (getMaxDim());
   } else {
      incDimBasis ();
   }
//   write (1);
}


//===========================================================================

void MRGLattice::incDimBasis()
// x_n = (a_1 x_{n-1} + a_2 x_{n-2} +...+ a_k x_{n-k}) mod m. On a Dim >= Order.
{
// trace( "=================================AVANT incDimBasis", -10);

   const int dim = getDim() + 1;
   m_v.setDim(dim);
   m_w.setDim(dim);

   for (int i = 1; i < dim; i++) {
      clear (m_vSI[0][i]);
      for (int j = 1; j <= m_order; j++) {
         conv (m_t1, m_v[i][dim - j]);
         m_t1 = m_t1 * m_aCoef[j];
         m_vSI[0][i] = m_vSI[0][i] + m_t1;
      }
      Modulo (m_vSI[0][i], m_m, m_vSI[0][i]);
      m_v[i][dim] = m_vSI[0][i];
   }

   for (int i = 1; i < dim; i++)
      m_v[dim][i] = 0;
   m_v[dim][dim] = m_m;

   for (int i = 1; i < dim; i++)
      m_w[i][dim] = 0;
   m_w[dim][dim] = 1;

   for (int j = 1; j < dim; j++) {
      clear (m_t1);
      for (int i = 1; i < dim; i++) {
         m_t2 = m_w[i][j];
         m_t2 *= m_vSI[0][i];
         m_t1 -= m_t2;
      }
      Quotient (m_t1, m_m, m_t1);
      m_w[dim][j] = m_t1;
   }

   m_v.setNegativeNorm(true);
   m_w.setNegativeNorm(true);
/*
 if (!checkDuality ())
      MyExit (1, "BUG");
      */
 // trace("=================================APRES incDimBasis", -10);
}


//===========================================================================

void MRGLattice::buildLaBasis (int d) {
   if (m_order > ORDERMAX)
      MyExit (1, "MRGLattice::buildLaBasis:   k > ORDERMAX");
   initStates();
   int IMax = m_lac.getSize();

   MVect b;
   b.SetLength(m_order + 1);
   Invert(m_aCoef, b, m_order);

   // b is the characteristic polynomial
   PolyPE::setM (m_m);
   PolyPE::setF(b);
   PolyPE pol;
   int ord = 0;

   // Construction d'un systeme generateur modulo m.
   for (int k = 1; k <= IMax; k++) {
      // pour chaque indice lacunaire
      conv (m_e, m_lac[k]);

      // x^m_e Mod f(x) Mod m
      pol.powerMod(m_e);
      pol.toVector (m_xi);

      ord = 0;
      for (int i = 1; i <= m_order; i++) {
         if (m_ip[i]) {
            ++ord;
            m_t5 = 0;
            for (int j = 1; j <= m_order; j++)
               m_t5 += m_sta[i][j] * m_xi[j - 1];
            m_wSI[ord][k] = m_t5;
         }
      }
   }

   /* On veut s'assurer que la base m_v soit triangulaire (pour satisfaire les
   conditions de l'article \cite{rLEC94e} [sec. 3, conditions sur V_i >= i])
   et de plein rang (on remplace les lignes = 0 par lignes avec m sur la
   diagonale). */
   Triangularization<BMat>(m_wSI, m_vSI, ord, IMax, m_m);
   CalcDual<BMat>(m_vSI, m_wSI, IMax, m_m);

   // Construire la base de dimension 1
   m_v[1][1] = m_vSI[1][1];
   m_w[1][1] = m_wSI[1][1];
   setDim(1);

   m_v.setNegativeNorm(true);
   m_w.setNegativeNorm(true);

   for (int i = 2; i <= d; i++)
      incDimLaBasis (IMax);

   // for debugging
   // trace("ESPION_1", 1);
}


//===========================================================================

void MRGLattice::incDimLaBasis(int IMax)
{
   const int dim = getDim();

   if (dim >= IMax) {
 /*     cout << " Dimension of the basis is too big:\n";
      cout << " Dim > Number of lacunary indices.";
      cout << endl; */
      MyExit(0, " Dimension of the basis is too big:\nDim > Number of lacunary indices.");
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

   setDim(dim + 1);
   m_v.setNegativeNorm(true);
   m_w.setNegativeNorm(true);
}


//===========================================================================

void MRGLattice::initStates ()
/*
 * Initialise la matrice carrée Sta. La matrice Sta est d'ordre égal à 
 * l'ordre du générateur. Elle contient un système de générateurs pour le 
 * groupe d'états considérés.
 */
{
   int maxDim = getMaxDim();
   clear (m_t2);

   if (m_latType == RECURRENT) {
      // check if a_k is relatively prime to m --> m_t1 = 1
      m_t1 = GCD (m_aCoef[m_order], m_m);
      m_t1 = abs(m_t1);
      set9 (m_t2);
   }

   if (m_latType == FULL || m_latType == PRIMEPOWER || (m_t1 == m_t2)) {
      // m_sta is set to identity matrix
      for (int i = 1; i <= m_order; i++) {
         for (int j = 1; j <= m_order; j++) {
            if (i != j)
               clear (m_sta[i][j]);
            else
               set9 (m_sta[i][j]);
         }
         m_ip[i] = true;
      }
      calcLgVolDual2 (m_lgm2);

   } else {
      if (m_latType == ORBIT) {
MyExit (1, "case ORBIT is not finished");

         MVect InSta;
         InSta.SetLength (1 + m_order);
         clear (m_sta[0][m_order]);
         for (int i = 1; i <= m_order; i++) {
            InSta[0] = m_aCoef[i] * InSta[m_order + 1 - i];
            m_sta[0][m_order] += InSta[0];
         }
         m_sta[0][m_order] -= InSta[m_order];
         for (int i = 1; i < m_order; i++)
            m_sta[0][i] = InSta[i + 1] - InSta[i];
         InSta.kill();

      } else if (m_latType == RECURRENT) {
         PolyPE::setM (m_m);
/* Je crois que la version sunos devait fonctionner correctement.
 * Je crois qu'Ajmal a créé des bugs dans la version mcs, qui se sont propagés
 * à mcs2, xds98, ..., c++. Je n'ai pas réussi à trouver l'erreur: les
 * résultats de plusieurs exemples dans sunos ne concordent pas avec les
 * résultats des versions subséquentes. Voir l'exemple 4 dans l'article
 *      AUTHOR="P. L'Ecuyer and R. Couture",
 *      TITLE="An Implementation of the Lattice and Spectral Tests for
 *             Multiple Recursive Linear Random Number Generators"
 *   A CORRIGER: comparer avec /u/lecuyer/stochas/latmrg/sunos/
 *
 * Voir la déf du m effectif comparé au vrai m dans LATIO. Je soupçonne
 * que cela pourrait être l'origine des erreurs.
 */

 MyExit (1, "case RECURRENT ne fonctionne pas");
         printf("ESPION_RECURRENT\n");
         MVect b;
         b.SetLength(m_order + 1);
         CopyVect (m_aCoef, b, m_order);
         PolyPE::reverse (b, m_order, 2);
         // b is the characteristic polynomial
         PolyPE::setF(b);
         PolyPE pol;

         // Must have 2^m_e > m^k to be sure to reach a recurrent state
         m_e = 3 + (int) (m_order * 0.5 * m_lgm2);
         pol.powerMod(m_e);
         pol.toVector (m_xi);

         m_sta[0][1] = m_xi[m_order - 1];
         for (int i = 2; i <= m_order; i++) {
            // Multiplier m_xi par X et reduire mod X^k - a1 X^{k-1} - ....

            m_xi[m_order] = m_xi[m_order - 1];
            for (int j = 1; j < m_order; j++) {
               // Coeff. de X^{m_order-j}.
               m_xi[m_order - j] = m_xi[m_order - j - 1];
               // ********* ATTENTION: MulMod (a, b, c, d) est très différent
               // pour les types long et ZZ (voir ZZ.txt). C'est pourquoi on
               // utilise MulMod (a, b, c) même s'il est plus lent. *********
               m_xi[m_order - j - 1] = MulMod (m_xi[m_order], m_aCoef[j], m_m);
               m_xi[m_order - j] += m_xi[m_order - j - 1];
            }
            // Coeff. constant.
            m_xi[0] = MulMod (m_xi[m_order], m_aCoef[m_order], m_m);
            m_sta[0][i] = m_xi[m_order - 1];
         }
      }

      for (int i = 1; i <= m_order; i++) {
         for (int j = 1; j <= m_order; j++)
            clear (m_sta[i][j]);
         m_ip[i] = false;
         m_sta[i][0] = m_sta[0][i];
      }
      insertion (m_sta);

      for (int k = 2; k <= m_order; k++) {
         // On passe a l'etat suivant.
         for (int j = 1; j < m_order; j++)
            m_sta[0][j] = m_sta[j + 1][0];
         clear (m_sta[0][m_order]);
         for (int i = 1; i <= m_order; i++) {
            m_t1 = m_aCoef[i] * m_sta[m_order - i + 1][0];
            m_sta[0][m_order] += m_t1;
         }
         Modulo(m_sta[0][m_order], m_m, m_sta[0][m_order]);
         // On memorise l'etat suivant.
         for (int i = 1; i < m_order; i++)
            swap (m_sta[i][0], m_sta[i + 1][0]);

         m_sta[m_order][0] = m_sta[0][m_order];
         insertion (m_sta);
      }

      lemme2 (m_sta);

      // Calcul de lgVolDual2
      double x;
      if (m_ip[1]) {
         conv(x, m_m / m_sta[1][1]);
         m_lgVolDual2[1] = 2.0 * Lg (x);
      } else
         m_lgVolDual2[1] = 0.0;

      int rmax = min(m_order, maxDim);
      for (int r = 2; r <= rmax; r++) {
         if (m_ip[r]) {
            conv(x, m_m / m_sta[r][r]);
            m_lgVolDual2[r] = m_lgVolDual2[r-1] + 2.0 * Lg (x);
         } else
            m_lgVolDual2[r] = m_lgVolDual2[r - 1];
      }

      for (int r = m_order + 1; r <= maxDim; r++)
         m_lgVolDual2[r] = m_lgVolDual2[r - 1];
   }
}


//===========================================================================

void MRGLattice::insertion (BMat & sta)
/*
 * Cette procedure insere le vecteur Sta[0] dans la matrice triangulaire
 * Sta. Si IP[i] = TRUE, l'entree diagonale sur la i-ieme ligne de Sta est
 * non-nulle (modulo m) et divise m. Sinon, la i-ieme ligne est
 * identiquement nulle. L'insertion doit preserver ces proprietes.
 * Le vecteur Sta[0] est altere au cours de l'operation.
 */
{
   for (int j = 1; j <= m_order; j++) {
      Modulo (sta[0][j], m_m, sta[0][j]);
      if (!IsZero (sta[0][j])) {
         if (!m_ip[j]) {
            Euclide (sta[0][j], m_m, m_t1, m_t2, m_t3, m_t4, sta[j][j]);
            for (int i = j + 1; i <= m_order; i++) {
               sta[j][i] = m_t1 * sta[0][i];
               Modulo (sta[j][i], m_m, sta[j][i]);
            }
            m_ip[j] = true;
            return;

         } else {
            Euclide (sta[j][j], sta[0][j], m_t1, m_t2, m_t3, m_t4, sta[j][j]);
            clear (sta[0][j]);
            for (int i = j + 1; i <= m_order; i++) {
               m_t5 = m_t1 * sta[j][i];
               m_t6 = m_t2 * sta[0][i];
               m_t7 = m_t3 * sta[j][i];
               m_t8 = m_t4 * sta[0][i];
               sta[j][i] = m_t5 + m_t6;
               Modulo (sta[j][i], m_m, sta[j][i]);
               sta[0][i] = m_t7 + m_t8;
            }
         }
      }
   }
}


//===========================================================================

void MRGLattice::lemme2 (BMat & sta)
/*
 * Cette procedure suppose que la matrice Sta est triangulaire. Si
 * IP[i] = TRUE, l'entree diagonale sur la i-ieme ligne de Sta est
 * non-nulle (modulo m) et divise m. Sinon, la i-ieme ligne est
 * identiquement nulle.
 */
{
   for (int i = 1; i <= m_order; i++) {
      if (m_ip[i]) {
         Quotient (m_m, sta[i][i], m_t1);
         m_t1 = abs (m_t1);
         if (m_t1 < m_m) {
            for (int j = 1; j < i; j++)
               sta[0][j] = sta[i][j];
            clear (sta[0][i]);
            for (int j = i + 1; j <= m_order; j++)
               sta[0][j] = m_t1 * sta[i][j];
            insertion (sta);
         }
      }
   }
}


//===========================================================================

void MRGLattice::initOrbit()
{
   MyExit (1, "MRGLattice::initOrbit n\'est pas terminée.");

/*
   for (int j = 0; j < J; j++) {
      for (int i = 1; i <= k; i++) {
         if (j == 0)
            clear(InSta[i]);
            for (int i1 = 1; i1 <= kj; i1++) {
               Multiply (aj [j,i1], VectSup [i-i1], SupT3);
               Add (SupT3, VectSup [i], VectSup [i])
               }
         Multiply (nj [j], VectSup [i], SupT3);
         Add (InSta [i], SupT3, InSta [i]);
         if (j == J - 1) {  Modulo (InSta [i], mm, InSta [i])  }
         }
      }
*/
}


//===========================================================================

}
