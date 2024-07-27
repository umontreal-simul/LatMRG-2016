#ifndef SEEKER_H
#define SEEKER_H
#include "latcommon/Types.h"


namespace LatMRG {

/**
 * Classes inheriting from this base class implement different constraints
 * that can be set on the coefficients (or multipliers) \f$a_i\f$ of a
 * recursive generator of order \f$k\f$, of the form
 * \f[
 *   x_n = (a_1x_{n-1} + a_2x_{n-2} + \cdots+ a_kx_{n-k}) \bmod m.
 * \f]
 * This base class itself applies no constraint on the coefficients.
 *
 */
class Seeker {
public:

/**
 * Constructor.
 */
Seeker() {}

   /**
    * Destructor. Does nothing for now.
    */
   virtual ~Seeker() {}
   virtual void init () { }
};

}
#endif
