#ifndef SHORTEST_DUAL_VECTOR_MERIT_H
#define SHORTEST_DUAL_VECTOR_MERIT_H
#include "HyperplaneDistanceMerit.h"
#include "latcommon/Normalizer.h"


namespace LatMRG {

/**
 * This class implements the computation of the shortest dual vector in a
 * lattice. Note that since the longer the length of the shortest dual
 * vector, the better the merit, the `compute` function returns a negative
 * value that should be minimized.
 *
 */
class ShortestDualVectorMerit : public HyperplaneDistanceMerit {
public:

   /**
    * Constructor.
    */
   ShortestDualVectorMerit (const LatCommon::Normalizer& normalizer, double power);

   /**
    * Destructor.
    */
   virtual ~ShortestDualVectorMerit() {}

   /**
    * Returns a string representation of the figure of merit.
    */
   virtual const std::string toString() const;

   /**
    * Returns minus the inverse value returned by
    * `HyperplaneDistanceMerit::compute()`.
    */
   virtual double compute (const std::vector<long>& a, int n,
                           const std::vector<int>& projection) const;
};

}
#endif
