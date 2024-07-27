#ifndef MRGLATTICEFACTORY_H
#define MRGLATTICEFACTORY_H
#include "latcommon/Const.h"
#include "latcommon/Types.h"
#include "latmrg/MRGLattice.h"
#include "latmrg/MRGComponent.h"


namespace LatMRG {

/**
 * This class is used to create <tt>MRGLattice</tt>’s from other types of
 * recurrences than a MRG or from a combination of several MRG lattices.
 *
 */
class MRGLatticeFactory {
public:

/**
 * Creates a `MRGLattice` from the combination of \f$J\f$
 * <tt>MRGComponent</tt>’s, with maximal dimension `maxDim`, lacunary indices
 * `Lac`, lattice type `lat` and vector norm `norm`. `Lac` can be `NULL` if
 * no lacunary indices are to be used. The combined MRG is calculated as
 * described in \cite rLEC96b&thinsp;.
 */
static MRGLattice * fromCombMRG (MRGComponent **comp, int J, int maxDim,
                                    BVect * Lac, LatCommon::LatticeType lat,
                                    LatCommon::NormType norm);

   /**
    * Creates a `MRGLattice` from a multiply-with-carry (MWC) recurrence,
    * with maximal dimension `maxDim`, lacunary indices `Lac`, lattice
    * type `lat` and vector norm `norm`. The MWC recurrence is
    * \f{align*}{
    *    x_n 
    *    & 
    *    = 
    *    (a_1 x_{n-1} + \cdots+ a_k x_{n-k} + c_{n-1}) d \mbox{ mod } b, 
    *  \\ 
    *   c_n 
    *    & 
    *    = 
    *  \lfloor(a_0 x_n + a_1 x_{n-1} + \cdots+ a_k x_{n-k} + c_{n-1}) / b \rfloor
    * \f}
    * where \f$b\f$ is a positive integer, \f$a_0,…,a_k\f$ are arbitrary
    * integers such that \f$a_0\f$ is relatively prime to \f$b\f$, and
    * \f$d\f$ is the multiplicative inverse of \f$-a_0\f$ mod \f$b\f$. The
    * MRG derived from such a MWC is defined by \f$ m = \sum^k_{l=0} a_l
    * b^l \f$ and \f$a\f$ is the inverse of \f$b\f$ in arithmetic modulo
    * \f$m\f$.
    */
   static MRGLattice * fromMWC (const MVect & a, const MScal & b, int maxDim,
                                int k, BVect *Lac, LatCommon::LatticeType lat,
                                LatCommon::NormType norm);

   /**
    * Same as above, but with no lacunary indices.
    */
   static MRGLattice * fromMWC (const MVect & a, const MScal & b, int maxDim,
                                int k, LatCommon::LatticeType lat,
                                LatCommon::NormType norm);
};

}
#endif
