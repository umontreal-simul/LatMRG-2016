#ifndef MRGCOMPONENTFACTORY_H
#define MRGCOMPONENTFACTORY_H
#include "latmrg/MRGComponent.h"
#include "latcommon/Types.h"
#include "latcommon/Util.h"


namespace LatMRG {

/**
 * This class is used to create <tt>MRGComponent</tt>s from other types of
 * recurrences [for instance, multiply-with-carry (MWC)].
 *
 */
class MRGComponentFactory {
public:

/**
 * Creates a `MRGComponent` from an multiply-with-carry type of recurrence.
 * The MWC recurrence has the form
 * \f{align*}{
 *    x_n 
 *    & 
 *    = 
 *    (a_1 x_{n-1} + \cdots+ a_k x_{n-k} + c_{n-1})d \mbox{ mod } b, 
 *  \\ 
 *   c_n 
 *    & 
 *    = 
 *  \lfloor(a_0 x_n + a_1 x_{n-1} + \cdots+ a_k x_{n-k} + c_{n-1}) / b \rfloor
 * \f}
 * where \f$b\f$ is a positive integer, \f$a_0,…,a_k\f$ are arbitrary
 * integers such that \f$a_0\f$ is relatively prime to \f$b\f$, and \f$d\f$
 * is the multiplicative inverse of \f$-a_0\f$ modulo \f$b\f$. The MRG
 * derived from such a MWC is defined by
 * \f[
 *   m = \sum^k_{i=0} a_i b^i
 * \f]
 * where \f$a\f$ is the inverse of \f$b\f$ in arithmetic modulo \f$m\f$.
 */
static MRGComponent * fromMWC (const MScal & b, const MVect & a, int k);
};

}
#endif
