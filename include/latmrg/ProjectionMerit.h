#ifndef PROJECTION_MERIT_H
#define PROJECTION_MERIT_H
#include "latcommon/Coordinates.h"
#include <string>
#include <vector>


namespace LatMRG {

/**
 * This abstract class is the base class of figures of merit for particular
 * projections. The best `FigureOfMerit` is usually the result of comparing
 * many <tt>ProjectionMerit</tt>’s.
 *
 * <div class="LatSoft-bigskip"></div>
 */
class ProjectionMerit {
public:

/**
 * Constructor.
 */
ProjectionMerit() {}

   /**
    * Destructor.
    */
   virtual ~ProjectionMerit() {}

   /**
    * Returns a string representation of the figure of merit.
    */
   virtual const std::string toString() const = 0;

   /**
    * Computes the value of the figure of merit for the given projection
    * of a lattice with generating vector \f$a_j/n =\f$ `a[j]/n`, for \f$j
    * = 1, 2, …, d\f$. Element `a[0]` is unused. It is assumed that the
    * value will be minimized, so a smaller value means a better merit.
    * The projection is specified by the coordinate indices contained in
    * `projection`, including `projection[0]`.
    */
   virtual double compute (const std::vector<long>& a, int n,
                           const LatCommon::Coordinates& projection) const = 0;
};

}
#endif
