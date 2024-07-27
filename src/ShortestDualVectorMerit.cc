#include "latmrg/ShortestDualVectorMerit.h"
#include "latcommon/Normalizer.h"

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <cfloat>

using namespace std;
using namespace LatCommon;

namespace LatMRG
{

//===========================================================================

ShortestDualVectorMerit::ShortestDualVectorMerit (const Normalizer& normalizer, double power)
      : HyperplaneDistanceMerit(normalizer, power)
{}


//===========================================================================

const std::string ShortestDualVectorMerit::toString() const
{
   return "minus the normalized length of the shortest dual vector";
}


//===========================================================================

double ShortestDualVectorMerit::compute (const std::vector<long>& a, int n, const std::vector<int>& projection) const
{
   return -1.0 / HyperplaneDistanceMerit::compute(a, n, projection);
}

//===========================================================================

}
