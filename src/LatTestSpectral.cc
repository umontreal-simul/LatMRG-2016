#include "latmrg/LatTestSpectral.h"
#include "latmrg/IntLattice.h"
#include "latcommon/Const.h"
#include "latmrg/Merit.h"
#include "latcommon/Util.h"
#include "latmrg/PolyPE.h"
#include "latmrg/LatticeTest.h"
#include "latcommon/Reducer.h"
#include <cmath>


using namespace std;
using namespace NTL;
using namespace LatCommon;

namespace LatMRG
{

//===========================================================================

LatTestSpectral::LatTestSpectral (const Normalizer * normal,
              LatMRG::IntLattice * lat): LatticeTest (lat)
{
   m_criter = SPECTRAL;
   m_normalizer = normal;
   int dim = lat->getMaxDim();
   m_boundL2.SetLength (1 + dim);
   m_S2toL2 = new double[1 + dim];
   SetZero (m_S2toL2, dim);
}


//===========================================================================

LatTestSpectral::~LatTestSpectral ()
{
   m_boundL2.kill();
   delete [] m_S2toL2;
}


//===========================================================================

void LatTestSpectral::initLowerBoundL2 (int dim1, int dim2)
{
   if (m_lat->getNorm () == L2NORM) {
      for (int i = dim1; i <= dim2; i++) {
         m_S2toL2[i] = m_normalizer->getGamma (i);
         m_S2toL2[i] *= exp2 (m_lat->getLgVolDual2 (i) / i);
      }

   } else if (m_lat->getNorm () == L1NORM) {
      if (m_dualF) {
         for (int i = dim1; i <= dim2; i++)
            m_S2toL2[i] = trunc(exp2 ((m_lat->getLgVolDual2 (i) / 2.0
                          + m_normalizer->getGamma(i)) / i));
      } else {
         // Je ne suis pas sûr que ce soit correct pour le primal
         for (int i = dim1; i <= dim2; i++)
            m_S2toL2[i] = exp2 ((m_lat->getLgVolDual2 (i) / 2.0
                          + m_normalizer->getGamma(i)) / i);
      }

   } else {
      for (int i = dim1; i <= dim2; i++) {
         m_S2toL2[i] = 1.0;
      }
   }
}


//===========================================================================

void LatTestSpectral::setLowerBoundL2 (double S2)
{
   setLowerBoundL2(S2, 0);
}

void LatTestSpectral::setLowerBoundL2 (double S2, const double* weights)
{
   int i;
   if (S2 <= 0.0) {
      for (i = m_fromDim; i <= m_toDim; i++)
         m_boundL2[i] = 0;
      return;
   }

   if (m_lat->getNorm () == L2NORM) {
      if (!m_dualF) {
         // On a multiplié les coordonnées du primal par m pour calculer
         // uniquement avec des entiers: meilleure précision.
         for (i = m_fromDim; i <= m_toDim; i++)
            conv(m_boundL2[i], m_lat->getM2 ()*S2*m_S2toL2[i]);

      } else {
         for (i = m_fromDim; i <= m_toDim; i++)
            conv(m_boundL2[i], S2*m_S2toL2[i]);
      }

   } else {
      if (!m_dualF) {
         for (i = m_fromDim; i <= m_toDim; i++) {
            conv(m_boundL2[i], S2*m_S2toL2[i]);
            // In Reducer, we compare the square length of vectors
            m_boundL2[i] = m_boundL2[i]*m_boundL2[i];
            // On a multiplié les coordonnées du primal par m pour calculer
            // uniquement avec des entiers: meilleure précision.
            m_boundL2[i] = m_boundL2[i]*m_lat->getM2 ();
         }

      } else {
         for (i = m_fromDim; i <= m_toDim; i++) {
            conv(m_boundL2[i], S2*m_S2toL2[i]);
            // In Reducer, we compare the square length of vectors
            m_boundL2[i] = m_boundL2[i]*m_boundL2[i];
         }
      }
   }

   // Merit values will be scaled by the inverse weight after normalization.
   // It is the scaled values that need to be compared, but the Reducer doesn't know
   // about scaling. So we scale the bounds below by the weight, which ensures the
   // ordering of the scaled values is preserved.
   // FIXME: should it be the squared weight for L1NORM?
   if (weights) {
      for (i = m_fromDim; i <= m_toDim; i++)
         m_boundL2[i] *= weights[i];
   }
}


//===========================================================================

bool LatTestSpectral::test (int fromDim, int toDim, double minVal[])
{
   return test (fromDim, toDim, minVal, 0);
}

// weights is the array of the weights of all projections as follows:
//    weights[i] is the weight for projection [1, ..., i]
// weights == 0 means unit weight for all projections
bool LatTestSpectral::test (int fromDim, int toDim, double minVal[], const double* weights)
{
   init ();

   m_fromDim = fromDim;
   m_toDim = toDim;
   resetFromDim (m_lat->getOrder (), fromDim);

   Reducer red (*m_lat);
   if (m_normalizer->getNorm () != m_lat->getNorm ()) {
      cout << "SpectralTest: conflict between NormType and Normalizer" <<
      endl;
      exit (EXIT_FAILURE);
   }

   double mr;
   conv (mr, m_lat->getM ());
   double temp;
   NScal te;
   double lgvv1 = 0.0;
   while (m_lat->getDim () < fromDim)
      m_lat->incDim ();

   if (m_S2toL2[fromDim] <= 0.0)
      initLowerBoundL2 (fromDim, toDim);
   setLowerBoundL2 (minVal[toDim], weights);   // same S2 for all dim
   red.setBoundL2 (m_boundL2, fromDim, toDim);

   while (true) {
      if (m_dualF)
         m_lat->dualize ();
      int dim = m_lat->getDim ();
      if (red.shortestVector (m_lat->getNorm ())) {
         // Calcul de D2. Pour Norm # L2NORM, suppose que VV est a jour.
         if (m_lat->getNorm () == L2NORM) {
            m_lat->getPrimalBasis ().updateScalL2Norm (1);
            conv (temp, m_lat->getPrimalBasis ().getVecNorm (1));
            if (!m_dualF) {
               conv(te, temp);
               te = te / m_lat->getM2 ();
               conv(temp, te);
            }

         } else {
            conv (temp, red.getMinLength ());
            if (!m_dualF)
               temp = temp / mr;
         }
         if (3 == m_detailF) {
            if (m_dualF) {
               dispatchBaseUpdate (m_lat->getDualBasis());
               dispatchBaseUpdate (m_lat->getPrimalBasis());
            } else {
               dispatchBaseUpdate (m_lat->getPrimalBasis());
               dispatchBaseUpdate (m_lat->getDualBasis());
            }
         }
         // weight factor:
         // require a higher merit for more important projections by dividing
         // their merit by the weight of the projection
         double weight = weights ? weights[dim] : 1.0;

         m_merit.getMerit (dim) = temp / weight;

         if (dim <= Normalizer::MAX_DIM) { // Calcul de S2.
            if (m_lat->getNorm () == L2NORM) {
               // Calcul du log(base 2) de ||V1||^2.
               lgvv1 = Lg (temp);
               if ((m_S2toL2[dim] <= 0.0) || (0 != std::isinf(m_S2toL2[dim]))) {
                   m_merit[dim] = exp2(lgvv1 - m_lat->getLgVolDual2 (dim)/dim)
                                  / m_normalizer->getGamma (dim);
               } else {
                  m_merit[dim] = temp / m_S2toL2[dim];
               }

            } else if (m_lat->getNorm () == L1NORM) {
               if ((m_S2toL2[dim] <= 0.0) || (std::isinf(m_S2toL2[dim]))) {
                  double tmp = exp2 ((m_lat->getLgVolDual2 (dim) / 2.0
                               + m_normalizer->getGamma(dim)) / dim);
                  // Je ne suis pas sûr que c'est correct pour le primal
                  m_merit[dim] = temp / trunc(tmp);
               } else {
                  m_merit[dim] = temp / m_S2toL2[dim];
               }

            } else
               m_merit[dim] = temp;

            m_merit[dim] /= weight;

            // Si on sait deja que ce gen. ne pourra etre retenu,
            // on le rejette tout de suite et on arrete le test.
            if ((m_maxAllDimFlag && m_merit[dim] < minVal[toDim])
                  || m_merit[dim] < minVal[dim]) {
           //    m_merit[dim] = 0.0;
               return false;
            }
         }
         else
            m_merit[dim] /= weight;

         prepAndDisp (dim);
         if (m_dualF)
            m_lat->dualize ();

      } else {
         m_merit[dim] = -1.0;
         return false;
      }

      if (dim == toDim)
         break;
      m_lat->incDim ();
   }

   return true;
}


//===========================================================================

void LatTestSpectral::init ()
{
   const int N = 3;
   string header[N];
   if (m_lat->getNorm () == L2NORM) {
      if (m_invertF)
         header[0] = "d_t";
      else
         header[0] = "l_t";
      header[1] = "S_t";
      header[2] = "Cumul CPU t(sec)";
      dispatchTestInit ("SPECTRAL", header, N);

   } else {
      if (m_invertF)
         header[0] = "d_t";
      else
         header[0] = "N_t";
//      header[1] = "N_t^*";
      header[1] = "S_t";
      header[2] = "Cumul CPU t(sec)";
      dispatchTestInit ("SPECTRAL", header, N);
   }
   timer.init ();
}


//===========================================================================

void LatTestSpectral::prepAndDisp (int dim)
{
   const int N = 3;
   double results[N];

   if (2 == m_detailF)
      dispatchBaseUpdate (m_lat->getPrimalBasis());
   else if (1 == m_detailF)
      dispatchBaseUpdate (m_lat->getPrimalBasis(), 1);

   if (m_lat->getNorm () == L2NORM) {
      results[0] = sqrt (m_merit.getMerit (dim));    // L_t
      if (m_invertF)
         results[0] = 1.0 / results[0];   // d_t = 1/L_t
      results[1] = sqrt (m_merit[dim]);
      results[2] = timer.val (Chrono::SEC);
      dispatchResultUpdate (results, N);

   } else {   // L1NORM
      results[0] = m_merit.getMerit (dim);
      if (m_invertF)
         results[0] = 1.0 / results[0];   // d_t = 1/N_t

#if 0
      // Calcule la valeur de N_t^*
      double x;
      conv (x, m_lat->getM());
      x = pow (x, (double) m_lat->getOrder());
      double y;
      if (m_dualF) {
         y = floor (pow (Factorial (dim) * x, 1.0 / dim));
      } else {
         x = 1.0/x;
         y = pow (Factorial (dim) * x, 1.0 / dim);
      }
      results[1] = y;
#endif

      results[1] = m_merit[dim];
      results[2] = timer.val (Chrono::SEC);
      dispatchResultUpdate (results, N);
   }
}

}
