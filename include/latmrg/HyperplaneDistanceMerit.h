#ifndef HYPERPLANE_DISTANCE_MERIT_H
#define HYPERPLANE_DISTANCE_MERIT_H
#include "ProjectionMerit.h"
#include "latcommon/Normalizer.h"


namespace LatMRG {

/**
 * This class implements the computation of the shortest dual vector in a
 * lattice.
 *
 * <div class="LatSoft-bigskip"></div>
 */
class HyperplaneDistanceMerit : public ProjectionMerit {
public:

/**
 * Constructor. The real number `power` is the exponent to which the
 * contribution for each contribution is raised.
 */
HyperplaneDistanceMerit (const LatCommon::Normalizer& normalizer, double power);

   /**
    * Destructor.
    */
   virtual ~HyperplaneDistanceMerit() {}

   /**
    * Returns a string representation of the figure of merit.
    */
   virtual const std::string toString() const;

   /**
    * See `ProjectionMerit`.
    */
   virtual double compute (const std::vector<long>& a, int n,
                           const std::vector<int>& projection) const;
protected:

   const LatCommon::Normalizer& m_normalizer;
   const double m_power;
};

}
#endif
