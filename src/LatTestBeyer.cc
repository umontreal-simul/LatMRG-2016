#include "latmrg/LatTestBeyer.h"
#include "latmrg/IntLattice.h"
#include "latcommon/Reducer.h"
#include "latmrg/Merit.h"
#include "latcommon/Const.h"
#include "latcommon/Util.h"
#include "latmrg/PolyPE.h"
#include "latmrg/LatticeTest.h"

using namespace std;
using namespace NTL;
using namespace LatCommon;

namespace LatMRG
{

//===========================================================================

LatTestBeyer::LatTestBeyer (LatMRG::IntLattice * lat): LatticeTest (lat)
{
   m_criter = BEYER;
}

//===========================================================================

bool LatTestBeyer::test (int fromDim, int toDim, double minVal[])
{
   init ();

   m_fromDim = fromDim;
   m_toDim = toDim;
   Reducer red (*m_lat);

   resetFromDim (m_lat->getOrder (), fromDim);
   while (m_lat->getDim () < fromDim)
      m_lat->incDim ();

   m_lat->dualize ();
   red.preRedDieter (0);
   m_lat->dualize ();

   while (true) {
      if (m_dualF)
         m_lat->dualize ();
      bool success = red.reductMinkowski (0);
      int dim = m_lat->getDim ();
      if (success) {
         m_lat->getPrimalBasis ().updateScalL2Norm (1);
         m_lat->getPrimalBasis ().updateScalL2Norm (dim);

         double x1, x2;        // si VV[1] et VV[dim] sont tres
         // grands, il faudrait envisager de changer x1 et x2 en xdouble.
         conv (x1, m_lat->getPrimalBasis ().getVecNorm (1));
         conv (x2, m_lat->getPrimalBasis ().getVecNorm (dim));
         m_merit[dim] = x1 / x2;

         // Si on sait deja que ce gen. ne pourra etre retenu,
         // on le rejette tout de suite et on arrete le test.
         if ((m_maxAllDimFlag && (m_merit[dim] < minVal[toDim]))
             || (m_merit[dim] < minVal[dim])) {
            m_merit[dim] = 0.0;
            return false;
         }

         prepAndDisp (dim);
         if (m_dualF)
            m_lat->dualize ();

      } else {
         m_merit[dim] = 0.0;
         return false;
      }

      if (dim == toDim)
         break;
      m_lat->incDim ();
   }

   return true;
}


//===========================================================================

void LatTestBeyer::init ()
{
   const int N = 2;
   string header[N];
   header[0] = "q_t";
   header[1] = "Cumul CPU t(sec)";
   dispatchTestInit ("Beyer", header, N);
   timer.init();
}


//===========================================================================

void LatTestBeyer::prepAndDisp (int dim)
{
   const int N = 2;
   double results[N];
   results[0] = sqrt (m_merit[dim]);
   results[1] = timer.val (Chrono::SEC);
   dispatchResultUpdate (results, N);
}

//===========================================================================

}
