#ifndef COEFFICIENT_H
#define COEFFICIENT_H
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
class Coefficient {
public:

/**
 * Constructor.
 */
Coefficient() {}

   /**
    * Destructor. Does nothing for now.
    */
   virtual ~Coefficient() {}

   /**
    * Sets the next coefficient `A[i]` \f$ = a_i= q\f$ to be tested. Index
    * \f$i\f$ is not changed for this base class, but may be changed for
    * subclasses depending on the constraints applied on the \f$a_i\f$’s.
    */
   virtual void set (const MScal & q, MVect & A, int & i) { A[i] = q; }
};

}
#endif
