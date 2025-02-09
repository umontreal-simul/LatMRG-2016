#ifndef RANK1LATTICE_H
#define RANK1LATTICE_H
#include "latcommon/Types.h"
#include "latcommon/Const.h"
#include "latmrg/IntLattice.h"


namespace LatMRG {

/**
 * This class implements a general rank 1 lattice basis. For the values
 * \f$a_1, a_2, …, a_d\f$ given, the \f$d\f$-dimensional lattice basis is
 * formed as:
 * \f[
 * \mathbf{b_1} = (a_1, a_2, …, a_d),\quad\mathbf{b_2} = (0, n, 0, …, 0),\quad…, \quad\mathbf{b_d} = (0, …, 0, n)
 * \f]
 * Without loss of generality, one may choose \f$a_1 = 1\f$.
 *
 * \warning There is some code duplication with LatCommon::Rank1Lattice.  This
 * should be fixed in the future.
 *
 */
class Rank1Lattice: public IntLattice {
public:

   /**
    * Constructor. \f$d\f$ represents the number of multipliers in the array
    * `a`.
    */
   Rank1Lattice (const MScal & n, const MVect & a, int d,
                    LatCommon::NormType norm = LatCommon::L2NORM);

   /**
    * Copy constructor.
    */
   Rank1Lattice (const Rank1Lattice & Lat);

   /**
    * Assigns `Lat` to this object.
    */
   Rank1Lattice & operator= (const Rank1Lattice & Lat);

   /**
    * Destructor.
    */
   ~Rank1Lattice();

   /**
    * Returns the vector of multipliers \f$a\f$ as a string.
    */
   std::string toStringCoef() const;

   /**
    * Builds the basis in dimension \f$d\f$.
    */
   void buildBasis (int d);

   /**
    * Increases the dimension by 1.
    */
   void incDim ();
protected:

   /**
    * Initializes the rank 1 lattice.
    */
   void init();

   /**
    * The multipliers of the rank 1 lattice rule.
    */
   MVect m_a;
};

}
#endif
