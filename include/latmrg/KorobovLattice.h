#ifndef LATCOMMON__KOROBOVLATTICE_H
#define LATCOMMON__KOROBOVLATTICE_H
#include "latcommon/Types.h"
#include "latcommon/Const.h"
#include "latmrg/IntLattice.h"


namespace LatMRG {

/**
 * This class implements lattice bases built from a Korobov lattice rule. For
 * a given \f$a\f$, a Korobov lattice basis is formed as follows:
 * \f[
 * \mathbf{b_1} = (1, a, a^2, …, a^{d-1}),\quad
 * \mathbf{b_2} = (0, n, 0, …, 0),\quad…,\quad
 * \mathbf{b_d} = (0, …, 0, n).
 * \f]
 * 
 * \remark **Pierre:** Reprogrammer \c incDim de façon efficace comme dans
 * \c MRGLattice
 * 
 */
class KorobovLattice: public IntLattice {
public:

   /**
    * Constructs a Korobov lattice with \f$n\f$ points, maximal dimension
    * `maxDim` using the norm `norm`.
    */
   KorobovLattice (const MScal & n, const MScal & a, int maxDim,
                      LatCommon::NormType norm = LatCommon::L2NORM);

   /**
    * Constructor. Same as above, except the lattice is formed as follow:
    * \f[
    * \mathbf{b_1} = (a^t, a^{t+1}, a^{t+2}, …, a^{t+d-1}),\qquad
    * \mathbf{b_2} = (0, n, 0, …, 0),\qquad…,\qquad
    * \mathbf{b_d} = (0, …, 0, n).
    * \f]
    */
   KorobovLattice (const MScal & n, const MScal & a, int dim, int t,
                   LatCommon::NormType norm = LatCommon::L2NORM);

   /**
    * Copy constructor.
    */
   KorobovLattice (const KorobovLattice & Lat);

   /**
    * Assigns `Lat` to this object.
    */
   KorobovLattice & operator= (const KorobovLattice & Lat);

   /**
    * Destructor.
    */
   ~KorobovLattice();

   /**
    * Returns the multiplier \f$a\f$ as a string.
    */
   std::string toStringCoef() const;

   /**
    * Builds the basis in dimension \f$d\f$.
    */
   void buildBasis (int d);

   /**
    * Increments the dimension of the basis by 1.
    */
   void incDim();

   /**
    * Increments the dimension of the basis by 1 by rebuilding the basis
    * from scratch. This is very slow. It can be used for verification of
    * the fast `incDim` method above.
    */
   void incDimSlow();
protected:

   /**
    * The multiplier of the Korobov lattice rule.
    */
   MScal m_a;

   /**
    * The shift applied to the lattice rule.
    */
   int m_shift;

   /**
    * Initialization.
    */
   void init();
};

}
#endif
