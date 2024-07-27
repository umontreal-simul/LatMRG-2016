/**
* MRGLatticeFactory.cc for ISO C++
* version 1.00
* authors: Hicham Wahbi
*          Frédérik Rozon
*          Richard Simard
*/

#include "latmrg/MRGLatticeFactory.h"
#include "latmrg/MRGLattice.h"

#include <cassert>

using namespace LatCommon;

namespace LatMRG
{

MRGLattice *MRGLatticeFactory::fromCombMRG (MRGComponent ** comp, int J,
      int maxDim, BVect * I, LatticeType type, NormType norm)
{
   MScal _m;
   int _k = 1;
   MScal _n[J];
   MScal d, e, f, g;

   conv (_m, 1);
   for (int j = 0; j < J; j++) {
      _k = std::max (_k, comp[j]->k);
      _m *= comp[j]->getM();
   }

   // Calcul de nj
   MScal tmp;
   for (int j = 0; j < J; j++) {
      Quotient (_m, comp[j]->getM(), tmp);
      Euclide (tmp, comp[j]->getM(), _n[j], d, e, f, g);
      if (g < 0)
         _n[j] = -_n[j];
      _n[j] *= tmp;
      comp[j]->nj = _n[j];
   }

   MVect _a;
   CreateVect (_a, _k);

   // Calcul de ai
   for (int i = 1; i <= _k; ++i) {
      clear (_a[i]);
      for (int j = 0; j < J; ++j) {
         if (comp[j]->k >= i) {
            tmp = comp[j]->a[i] * _n[j];
            _a[i] += tmp;
         }
      }
      Modulo (_a[i], _m, _a[i]);
   }
   for (int j = 0; j < J; ++j) {
      tmp = _a[1] - comp[j]->a[1];
      Modulo (tmp, comp[j]->getM(), tmp);
      assert (0 == tmp);
   }

   // Calcul de rho, lossrho et rhoj
   MScal rho;
   MScal lossRho;
   rho = 1;
   lossRho = 1;

   for (int j = 0; j < J; j++) {
      comp[j]->rho = power (comp[j]->getM(), comp[j]->k) - 1;
//      Euclide (g, rho, tmp, d, e, f, comp[j]->rho);
      Euclide (comp[j]->rho, rho, tmp, d, e, f, g);
      lossRho *= g;
      rho *= comp[j]->rho;
      Quotient (rho, g, rho);
   }

   MRGLattice *lat;

   if (I != 0)
      lat = new MRGLattice (_m, _a, maxDim, _k, *I, type, norm);
   else
      lat = new MRGLattice (_m, _a, maxDim, _k, type, norm);

   lat->setRho (rho);
   lat->setLossRho (lossRho);
   DeleteVect (_a);
   lat->comp.clear();
   lat->comp.reserve(J);
   MRGComponent *mycomp;
   for (int j = 0; j < J; j++) {
      mycomp = new MRGComponent (*comp[j]);
      lat->comp.push_back(mycomp);
   }
   return lat;
}


//=========================================================================
#if 0

MRGLattice *MRGLatticeFactory::fromCombMRG (const MScal * m,
      const MMat & coef, int MaxDim, int *k, int J, BVect * I, LatticeType lat,
      NormType norm)
{
   MScal _m;
   int _k = 0;
   MScal _n[J];
   MVect _a;

   MScal d, e, f, g;

   conv (_m, 1);

   // Calcul de m et k;
   for (int i = 0; i < J; ++i) {
      _k = std::max (_k, k[i]);
      _m *= m[i];
   }

   // Calcul de nj
   MScal t1;
   for (int i = 0; i < J; ++i) {
      t1 = _m / m[i];

      Euclide (t1, m[i], _n[i], d, e, f, g);

      if (g < 0) {
         _n[i] = -_n[i];
      }
      _n[i] *= t1;

   }

   CreateVect (_a, _k);
   MScal tmp;
   // Calcul de ai
   for (int i = 1; i <= _k; ++i) {
      clear (_a[i]);
      for (int j = 0; j < J; ++j) {
         if (k[j] >= i) {
            tmp = coef[j][i] * _n[j];
            _a[i] += tmp;
         }
      }
      Modulo (_a[i], _m, _a[i]);
   }

   if (I != 0) {
      return new MRGLattice (_m, _a, MaxDim, _k, *I, lat, norm);
   } else {
      return new MRGLattice (_m, _a, MaxDim, _k, lat, norm);
   }
}


//=========================================================================

MRGLattice *MRGLatticeFactory::fromCombMRG (const MScal m[],
      const MMat & coef, int MaxDim, int k[], int J, LatticeType lat,
      NormType norm)
{
   return fromCombMRG (m, coef, MaxDim, k, J, 0, lat, norm);
}
#endif


//=========================================================================

MRGLattice *MRGLatticeFactory::fromMWC (const MVect & a, const MScal & b,
       int MaxDim, int k, LatticeType lat, NormType norm)
{
   return fromMWC (a, b, MaxDim, k, 0, lat, norm);
}


//=========================================================================

MRGLattice *MRGLatticeFactory::fromMWC (const MVect & a, const MScal & b,
          int MaxDim, int k, BVect * I, LatticeType lat_t, NormType norm)
{
   MScal _m;
   MScal _b;
   MVect _a;

   MScal d, e, f, g;

   conv (_m, 0);
   conv (_b, 1);
   CreateVect (_a, 1);

   // Calcul de m
   for (int i = 0; i <= k; i++) {
      _m += a[i] * _b;
      _b *= b;
   }

   // Calcul de a
   Euclide (b, _m, _a[1], d, e, f, g);

   if (g < 0) {
      _a[1] = -_a[1];
   }

   MRGLattice *lat;
   if (I) {
      lat = new MRGLattice (_m, _a, MaxDim, 1, *I, lat_t, norm);
   } else {
      lat = new MRGLattice (_m, _a, MaxDim, 1, lat_t, norm);
   }
   DeleteVect (_a);
   return lat;
}


//=========================================================================

}
