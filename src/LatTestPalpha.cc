#include "latmrg/LatTestPalpha.h"
#include "latmrg/LatticeTest.h"
#include "latmrg/PalphaLCG.h"
#include "latmrg/MRGComponent.h"
#include "latmrg/IntLattice.h"
#include "latcommon/IntFactor.h"
#include "latcommon/Const.h"
#include "latmrg/Merit.h"
#include "latcommon/Util.h"
#include "latmrg/PolyPE.h"


using namespace std;
using namespace NTL;
using namespace LatCommon;

namespace LatMRG
{

//===========================================================================

void LatTestPalpha::init ()
{
   const int N = 3;
   string header[N];
   if (m_config->calcPalpha == NORMPAL) {
      header[0] = "P_a";
      header[1] = "P_a/B_a";
      header[2] = "Cumul CPU t(sec)";
      dispatchTestInit ("Palpha", header, N);
   } else {
      header[0] = "P_a";
      header[1] = "Cumul CPU t(sec)";
      dispatchTestInit ("Palpha", header, N - 1);
   }
   timer.init();
}


//===========================================================================

void LatTestPalpha::prepAndDisp (int dim)
{
   const int N = 3;
   double results[N];
   if (m_config->calcPalpha == NORMPAL) {
      results[0] = m_merit.getMerit (dim);
      results[1] = m_merit[dim];
      results[2] = timer.val(Chrono::SEC);
      dispatchResultUpdate (results, N);
   } else {
      results[0] = m_merit[dim];
      results[1] = timer.val(Chrono::SEC);
      dispatchResultUpdate (results, N - 1);
   }
}


//===========================================================================

LatTestPalpha::LatTestPalpha (Normalizer * normal, LatMRG::IntLattice * lat): LatticeTest (lat)
{
   m_criter = PALPHA;
   m_bound = normal;
//   m_config = config;
//   m_fromDim = m_config->td[0];
//   m_toDim = m_config->td[1];
   m_dualF = false;
   m_maxAllDimFlag = false;
}


//===========================================================================
/*
bool LatTestPalpha::test (double MinVal[])
{
   return test (m_config->fromDim, m_config->toDim, MinVal);
}
*/

//===========================================================================

bool LatTestPalpha::test (int fromDim, int toDim, double minVal[])
{
   m_fromDim = fromDim;
   m_toDim = toDim;
   init ();

   if (m_config->verifyM) {
      if (PRIME == IntFactor::isPrime (m_config->comp[0]->getM(), 50)) {
         cout << "Verify prime m:   true" << endl;
         m_config->primeM = true;
      } else {
         cout << "Verify prime m:   false" << endl;
         m_config->primeM = false;
      }
   }

   if (m_config->verifyP) {
      MRGComponent mrg (m_config->comp[0]->getM(), 1, DECOMP, 0, DECOMP_PRIME,  0);
      if (mrg.maxPeriod (m_config->comp[0]->a)) {
         cout << "Verify maximal period:   true" << endl;
         m_config->maxPeriod = true;
      } else {
         cout << "Verify maximal period:   false" << endl;
         m_config->maxPeriod = false;
      }
      cout << endl;
   }

   PalphaLCG palpha (*m_config);
   const int alpha = m_config->alpha;
   double x;
   int dim;

   if (m_config->calcPalpha == PAL || m_config->calcPalpha == NORMPAL) {

      // SI LE MODULO EST PREMIER  // Si periode maximale
      if (m_config->primeM && m_config->maxPeriod) {
         for (dim = fromDim; dim <= toDim; dim++) {
            switch (alpha) {
            case 2:
               x = palpha.calcPalpha2 (dim);
               break;
            case 4:
               x = palpha.calcPalpha4 (dim);
               break;
            case 6:
               x = palpha.calcPalpha6 (dim);
               break;
            case 8:
               x = palpha.calcPalpha8 (dim);
               break;
            default:
               cerr << " Valeur de alpha invalide" << endl;
               return false;
            }

            if (m_config->calcPalpha == PAL)
               m_merit[dim] = x;
            else {
               m_merit.getMerit(dim) = x;
               if (m_bound->getCst(dim) < 0.0)
                  m_merit[dim] = -1.0;
               else
                  m_merit[dim] = x / m_bound->getCst(dim);
            }
            prepAndDisp (dim);
         }

      } else {
         for (dim = fromDim; dim <= toDim; dim++) {
            switch (alpha) {
            case 2:
               x = palpha.calcPalpha2PerNonMax (dim);
               break;
            case 4:
               x = palpha.calcPalpha4PerNonMax (dim);
               break;
            case 6:
               x = palpha.calcPalpha6PerNonMax (dim);
               break;
            case 8:
               x = palpha.calcPalpha8PerNonMax (dim);
               break;
            default:
               cerr << " Valeur de alpha invalide" << endl;
               return false;
            }
            if (m_config->calcPalpha == PAL)
               m_merit[dim] = x;
            else {
               m_merit.getMerit(dim) = x;
               if (m_bound->getCst(dim) < 0.0)
                  m_merit[dim] = -1.0;
               else
                  m_merit[dim] = x / m_bound->getCst(dim);
            }
            prepAndDisp (dim);
         }
      }

   } else if (m_config->calcPalpha == BAL) {
      for (dim = fromDim; dim <= toDim; dim++) {
         x = m_bound->getCst(dim);
         m_merit[dim] = x;
         prepAndDisp (dim);
      }

   }

   return true;
}

}
