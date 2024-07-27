#include "latmrg/SpectralMeritClassic.h"
#include "latcommon/Types.h"
#include "latmrg/Rank1Lattice.h"
#include "latmrg/LatTestSpectral.h"
#include "latmrg/TestProjections.h"
#include "latcommon/UniformWeights.h"
#include "latcommon/Util.h"

#include <stdexcept>
#include <sstream>


using namespace std;
using namespace LatCommon;

namespace LatMRG
{

const LatCommon::UniformWeights unitWeights(1.0);

//===========================================================================

SpectralMeritClassic::SpectralMeritClassic (long n, int dim)
{
   m_n = n;
   m_dim = dim;

   m_weights = &unitWeights;
   m_normalizer = 0;
   m_coordLimits = 0;
   m_debug = false;
   m_stationary = false;
   m_meritLowerLimits = new double[m_dim + 1];
   resetLimits();
}


//===========================================================================

SpectralMeritClassic::~SpectralMeritClassic()
{
   delete[] m_meritLowerLimits;
}

//===========================================================================

void SpectralMeritClassic::setCoordLimits(int* coordLimits)
{
   m_coordLimits = coordLimits;
}

//===========================================================================

void SpectralMeritClassic::setWeights(const Weights* weights)
{
   m_weights = weights;
}

//===========================================================================

void SpectralMeritClassic::setNormalizer(Normalizer* normalizer)
{
   m_normalizer = normalizer;
}

//===========================================================================

void SpectralMeritClassic::setProjectionStationary(bool stationary)
{
   m_stationary = stationary;
}

//===========================================================================

void SpectralMeritClassic::setDebug(bool debug)
{
   m_debug = debug;
}

//===========================================================================

void SpectralMeritClassic::resetLimits()
{
   for (int i = 0; i <= m_dim; i++)
      m_meritLowerLimits[i] = 0;
   m_last_dim = 0;
}

//===========================================================================

const string SpectralMeritClassic::toString () const
{
   //string S(FigureOfMerit::toString ());
   ostringstream os;
   os << " coord. limits = [";
   for (int i = 0; i <= m_dim; i++)
      os << " " << m_coordLimits[i];
   os << " ]" << endl;
   os << " stationary = " << boolalpha << m_stationary << endl;
   return "SpectralMeritClassic:\n" + os.str ();
}


//===========================================================================

double SpectralMeritClassic::compute (long z[], int dim)
{
   if (!m_coordLimits)
      throw std::invalid_argument("SpectralMeritClassic: compute() called with coordinate limits unset");
   if (!m_normalizer)
      throw std::invalid_argument("SpectralMeritClassic: compute() called with normalizer unset");

   // The previous best merit values do not have sense if we consider a different dimension.
   // Different dimensions mean different problems.
   if (dim != m_last_dim)
      resetLimits();
   m_last_dim = dim;

   MVect big_z;
   big_z.SetLength (dim + 1);
   for (int i = 0; i <= dim; i++)
      big_z[i] = z[i];

   MScal big_n;
   conv (big_n, m_n);

   // in case of CBC search, we must stop at dimension dim
   int coordLimitsReduced[dim + 1];
   for (int i = 0; i <= dim; i++)
      coordLimitsReduced[i] = std::min (dim, m_coordLimits[i]);

   Rank1Lattice master (big_n, big_z, dim, L2NORM);
   Rank1Lattice work (master);
   master.buildBasis (dim);
   LatTestSpectral test (m_normalizer, &work);
   TestProjections proj (&master, &work, &test, coordLimitsReduced, dim);
   proj.setDualFlag (true);
   proj.setPrintF (m_debug);

   // FIXME: propagate lower limits
   if (m_debug) {
      cout << "-------------------------------------------------------\n";
      cout << " z = " << LatCommon::toString (z, 1, dim) << endl;
      cout << " lower limits:" << LatCommon::toString(m_meritLowerLimits, 1, dim) << endl;
   }

   // FIXME: check that this works with CBC as well

   // FIXME: is this really useful anywhere?
   bool alwaysWithLastCoord = false;
   // compute the weighted normalized length of the shortest vector in the dual
   double merit = proj.run (m_stationary, alwaysWithLastCoord, m_meritLowerLimits, *m_weights);

   // Update lower limits.
   // Same for all dimensions: the dimension doesn't matter, all that counts is if we do better than
   // the current best merit.
   if (merit > m_meritLowerLimits[dim])
      for (int i = 0; i <= m_dim; i++)
         m_meritLowerLimits[i] = merit;

   // we want to maximize merit, so we minimize -merit:
   return merit < 0 ? 1e100 : -merit;
}


}
