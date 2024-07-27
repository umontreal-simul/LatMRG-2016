#ifndef FIGURE_OF_MERIT_H
#define FIGURE_OF_MERIT_H
#include <string>
#include <vector>


namespace LatMRG {

/**
 * This abstract class is the base class of all figures of merit that allow
 * choosing the best lattices or best generators according to some criterion.
 *
 * <div class="LatSoft-bigskip"></div>
 */
class FigureOfMerit {
public:

/**
 * Constructor.
 */
FigureOfMerit() {}

   /**
    * Destructor.
    */
   virtual ~FigureOfMerit() {}

   /**
    * Returns a string representation of the figure of merit.
    */
   virtual const std::string toString() const = 0;

   /**
    * Computes the value of the figure of merit for a lattice with
    * generating vector \f$a_j/n =\f$ `a[j]/n`, for \f$j = 1, 2, …, d\f$.
    * Element `a[0]` is unused. It is assumed that the value will be
    * minimized, so a smaller value means a better merit. If
    * `abortThreshold` is given and non-zero, computation is aborted
    * whenever one becomes certain that the final value will be larger
    * than `abortThreshold`.
    */
   virtual double compute (const std::vector<long>& a, int n,
                           double abortThreshold = 0) const = 0;
};

}
#endif