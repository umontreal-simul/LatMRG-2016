#include "latmrg/TestProjections.h"
#include "latmrg/LatticeTest.h"
#include "latcommon/Util.h"
#include "latmrg/Writer.h"
#include "latmrg/WriterRes.h"
#include "latcommon/UniformWeights.h"
#include "latmrg/ProjIteratorSuccCoords.h"
#include "latmrg/ProjIteratorNonSuccCoords.h"

#include <cstring>
#include <sstream>
#include "latcommon/CoordinateSets.h"

using namespace std;
using namespace LatCommon;


namespace
{
const char* espaceM = " \t\t   ";
const char* espaceG = "                                                   ";

const LatCommon::UniformWeights unitWeights(1.0);

const string formatIndices (const Coordinates & ens)
{
   ostringstream os;
   os << "  ";
   Coordinates::const_iterator it = ens.begin ();
   os << *it;
   ++it;
   while (it != ens.end ()) {
      os << "," << *it;
      ++it;
   }
   os << espaceM;
   return os.str ();
}


void printMerit (bool racF, double merit, LatMRG::Writer * rw)
{
   if (racF) {
      rw->writeString (espaceG);
      rw->writeString ("min merit = ");
      rw->writeDouble (LatCommon::mysqrt (merit));
   } else {
      rw->writeString (espaceG);
      rw->writeString ("min merit = ");
      rw->writeDouble (merit);
   }
   rw->newLine ();
   rw->newLine ();
}


void printMerLen (bool invF, bool racF, double len, double merit,
                  LatMRG::Writer * rw)
{
   double x = len;
   if (racF)
      x = LatCommon::mysqrt (len);
   if (invF)
      x = 1.0 / x;

   if (racF) {
      rw->writeDouble (x);
      rw->writeString (espaceM);
      rw->writeDouble (LatCommon::mysqrt (merit));
   } else {
      rw->writeDouble (x);
      rw->writeString (espaceM);
      rw->writeDouble (merit);
   }
   rw->newLine ();
}

}                                 // namespace


//===========================================================================

namespace LatMRG
{

//===========================================================================

TestProjections::TestProjections (LatMRG::IntLattice * master, LatMRG::IntLattice * lattice,
                                  LatticeTest * test, int td[], int d)
{
   if (d <= 0)
      MyExit(1, "   TestProjections:   d <= 0");
   for (int i = 2; i <= d; i++) {
      if (td[i] >= 100)
         MyExit(1, "   TestProjections:   td[i] >= 100");
   }
   m_d = d;
   m_td = new int[d + 1];
   memcpy (m_td, td, (d + 1) * sizeof (int));
   int maxDim = 0;
   for (int i = 0; i <= d; i++)
      maxDim = max (maxDim, m_td[i]);
   m_weightsTemp = new double[maxDim + 1];
   for (int i = 0; i <= maxDim; i++)
      m_weightsTemp[i] = 1.0;
   m_master = master;
   m_lattice = lattice;
   m_test = test;
   m_dualF = test->getDualFlag ();
   m_printF = true;
   m_invertF = test->getInvertFlag ();
   if ((SPECTRAL == test->getCriterion() && L2NORM == master->getNorm ())
         || (BEYER == test->getCriterion()))
      m_racF = true;
   else
      m_racF = false;
   m_writer = new WriterRes (&cout);
   m_wrFlag = true;
}


//===========================================================================

TestProjections::~TestProjections ()
{
   delete[] m_td;
   delete[] m_weightsTemp;
   if (m_wrFlag) 
      delete m_writer;
}


//===========================================================================

void TestProjections::setOutput (Writer * rw)
{
   delete m_writer;
   m_wrFlag = false;
   m_writer = rw;
}


//===========================================================================

void TestProjections::setDualFlag (bool dual)
{
   m_dualF = dual;
}


//===========================================================================

void TestProjections::setPrintF (bool flag)
{
   m_printF = flag;
}


//===========================================================================

void TestProjections::build (const Coordinates & proj)
{
   m_master->buildProjection (m_lattice, proj);
}

//===========================================================================

int TestProjections::calcNumProjections (bool stationary, bool forceLast)
{
   // Les projections sur dimensions successives
   m_numproj = m_td[1] - m_td[0] + 1;

   if (!stationary) {
      for (int order = 2; order <= m_d; order++) {
         int maxCoord = min (m_td[order], m_master->getMaxDim ());
         if (maxCoord <= order)
            continue;
         for (ProjIteratorSuccCoords projit(2, maxCoord, order, order, stationary, forceLast); projit; ++projit)
            m_numproj++;
      }
   }

   // test for non-successive coordinates
   // loop over projection orders
   for (int order = 2; order <= m_d; order++) {
      int maxCoord = min (m_td[order], m_master->getMaxDim ());
      if (maxCoord <= order)
         continue;
      for (ProjIteratorNonSuccCoords projit(1, maxCoord, order, order, stationary, forceLast); projit; ++projit)
         m_numproj++;
   }

   //! if (false == stationary) {
   //!    // Les autres projections sur dimensions successives
   //!    for (int i = 2; i <= m_d; i++) {
   //!       int dim = min (m_td[i], m_master->getMaxDim ());
   //!       initSuccDims (i, dim, forceLast);
   //!       do {
   //!          m_numproj++;
   //!       } while (nextSuccDims (stationary, forceLast, i, dim));
   //!    }
   //! }

   //! // Les projections sur dimensions non successives
   //! for (int i = 2; i <= m_d; i++) {
   //!    int dim = min (m_td[i], m_master->getMaxDim ());
   //!    if (initNonSuccDims (i, dim, forceLast)) {
   //!       do {
   //!          m_numproj++;
   //!       } while (nextNonSuccDims (stationary, forceLast, i, dim));
   //!    }
   //! }
   return m_numproj;
}


//=========================================================================

double TestProjections::run (bool stationary, bool forceLast, double minVal[])
{
   return run(stationary, forceLast, minVal, unitWeights);
}


double TestProjections::run (bool stationary, bool forceLast, double minVal[], const Weights& weights)
{
   // we assume that m_td[1] >= m_td[i] for all i
   int maxDim = m_td[1];
   double merit = 1.0e100;

   //! cout << "minVal: " << minVal[0];
   //! for (int j = 1; j < maxDim; j++)
   //!    cout << "," << minVal[j];
   //! cout << endl;

   // Le test pour les projections sur dimensions successives
   m_lattice->buildBasis (m_td[0] - 1);
   int minDim = m_td[0];
   m_test->setDualFlag (m_dualF);
   m_test->setMaxAllDimFlag (true);

   // set the temporary weights for successive dimensions
   for (ProjIteratorSuccCoords projit(1, maxDim, 1, maxDim, true, false); projit; ++projit) {
      if ((int)projit->size() >= minDim)
         m_weightsTemp[projit->size()] = weights.getWeight(*projit);
      else
         m_weightsTemp[projit->size()] = 1;
   }
   m_test->test (minDim, maxDim, minVal, m_weightsTemp);

   // ATTENTION: si le test s'est terminé prématurément parce que le réseau
   // est mauvais, les valeurs de mérites ci-après sont n'importe quoi.
   merit = m_test->getMerit().getST (minDim, maxDim);
   m_numproj = maxDim - minDim + 1;

   if (m_printF) {
      //cout << "------------------------------------------" << endl;
      //cout << " a = " << m_lattice->toStringCoef () << endl;
      m_writer->writeString (" ");
      m_writer->getStream() << weights;
      m_writer->newLine ();
      m_writer->writeString (
          m_test->getMerit ().toString(minDim, maxDim, m_racF, m_invertF));
      printMerit (m_racF, merit, m_writer);
   }
   if (merit < minVal[maxDim] || m_d <= 1)
      return merit;

   if (!stationary) {
      // test for successive coordinates
      // loop over projection orders
      for (int order = 2; order <= m_d; order++) {

         int maxCoord = min (m_td[order], m_master->getMaxDim ());
         // if maxCoord < order, there are no projections to consider
         // if maxCoord == order, the merit has already been computed above
         if (maxCoord <= order)
            continue;

         ProjIteratorSuccCoords projit(2, maxCoord, order, order, stationary, forceLast);
         merit = std::min(merit, run(projit, minVal, weights));

         if (merit < minVal[order])
            return merit;

         if (m_printF)
            printMerit (m_racF, merit, m_writer);
      }
   }

   // test for non-successive coordinates
   // loop over projection orders
   for (int order = 2; order <= m_d; order++) {

      int maxCoord = min (m_td[order], m_master->getMaxDim ());
      // if maxCoord <= order, there are no projections with non-successive indices to consider
      if (maxCoord <= order)
         continue;

      ProjIteratorNonSuccCoords projit(1, maxCoord, order, order, stationary, forceLast);
      merit = std::min(merit, run(projit, minVal, weights));

      if (merit < minVal[order])
         return merit;

      if (m_printF)
         printMerit (m_racF, merit, m_writer);
   }

   m_test->getMerit ().setWorstMerit (merit);
   return merit;
}


//===========================================================================

double TestProjections::run (ProjIterator& projit, double minVal[], const Weights& weights)
{
   double minMerit = 1.0e100;

   while (projit) {
      int dim = (int)projit->size();
      int maxCoord = min (m_td[dim], m_master->getMaxDim ());
      if (maxCoord <= dim) // if equal, already computed
         continue;

      if (m_printF)
         m_writer->writeString (formatIndices(*projit));

      m_numproj++;

      m_weightsTemp[dim] = weights.getWeight(*projit);
      m_master->buildProjection (m_lattice, *projit);
      m_test->test (dim, dim, minVal, m_weightsTemp);
      double len = m_test->getMerit ().getMerit (dim);
      double curMerit = m_test->getMerit ().getNormVal (dim);
      if (m_printF)
         printMerLen (m_invertF, m_racF, len, curMerit, m_writer);
      if (curMerit < minVal[dim])
         return curMerit;
      minMerit = std::min (minMerit, curMerit);
      ++projit;
   }
   return minMerit;
}

//===========================================================================

}
