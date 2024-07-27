#ifndef LATTESTSPECTRAL_H
#define LATTESTSPECTRAL_H
#include "latmrg/IntLattice.h"
#include "LatticeTest.h"
#include "latcommon/Reducer.h"
#include "latcommon/Normalizer.h"
#include "latcommon/Weights.h"
#include "latcommon/Const.h"


namespace LatMRG {

/**
 * This class implements the *spectral* test. It implements the abstract
 * class `LatticeTest`. The figure of merit for this test is the length of
 * the shortest vector in the *primal* or in the *dual* lattice computed with
 * different norms. For the standard spectral test, the figure of merit is
 * based on the length of the shortest non-zero vector in the *dual* lattice,
 * using the \f${\mathcal{L}}_2\f$ norm to compute the length of vectors, and
 * the inverse of this length gives the maximal distance between successive
 * hyperplanes covering all the points in the *primal* lattice. If one
 * computes the length of the shortest non-zero vector in the *dual* lattice
 * using the \f${\mathcal{L}}_1\f$ norm instead, one obtains the minimal
 * number of hyperplanes covering all the points of the *primal* lattice.
 *
 * The main program is obtained by compiling the `LatMain.cc` file (see the
 * description of module `LatMain`, on page (FIXME: page#), for how to run a
 * program).
 *
 */
class LatTestSpectral : public LatticeTest {
public:

   /**
    * Constructor. The *spectral* test will be applied to the lattice `lat`
    * using normalizer `normal` to normalize the figure of merit.
    */
   LatTestSpectral (const LatCommon::Normalizer * normal,
                       IntLattice * lat);

   /**
    * Destructor.
    */
   ~LatTestSpectral ();

   /**
    * Applies the spectral test for dimensions varying from `fromDim` to
    * `toDim`. Whenever the normalized value of the merit is smaller than
    * `minVal` for any dimension, the method returns `false` immediately.
    * The method returns `false` if the test was interrupted for any
    * reason before completion, and it returns `true` upon success. The
    * results of the last test are kept in <tt>m_merit</tt>.
    */
   bool test (int fromDim, int toDim, double minVal[]);

   /**
    * Similar to `test` above, but with the weights `weights`. `weights`
    * is the array of the weights of all projections defined as follows:
    *
    * `weights[0]` is the weight of projections [1, ... ,
    * <tt>fromDim</tt>];<br><tt>weights[1]</tt> is the weight of
    * projections [1, ... , `fromDim` + 1];<br>
    * \f$\vdots\f$ <br><tt>weights</tt>[<tt>toDim - fromDim</tt>] is the
    * weight of projections [1, ... , <tt>toDim</tt>].
    *
    * If `weights` = 0, it means unit weight for all projections.
    */
   bool test (int fromDim, int toDim, double minVal[], const double* weights);

   /**
    * Sets the lower bound on the square length of the shortest vector in
    * each dimension, based on the spectral value `S2`.
    */
   void setLowerBoundL2 (double S2);

   /**
    * Similar to `setLowerBoundL2` above, but with the weights `weights`.
    */
   void setLowerBoundL2 (double S2, const double* weights);

   /**
    * Returns the normalizer used in this test.
    */
   const LatCommon::Normalizer* getNormalizer() { return m_normalizer; }

private:

   /**
    * The normalizer used to normalize the figure of merit.
    */
   const LatCommon::Normalizer* m_normalizer;

   /**
    * The lower bound on the square length of the shortest vector in each
    * dimension. As soon as a vector of length smaller than this bound is
    * found, the search for the shortest vector in this lattice is stopped
    * and the lattice is rejected.
    */
   NVect m_boundL2;

   /**
    * Initializes the constants <tt>m_S2toL2</tt> below, necessary to
    * compute the lower bounds in `setLowerBoundL2`, for all dimensions
    * \f$d\f$ such that `dim1` \f$\le d \le\f$ `dim2`. This function must
    * be called only after the lattice has been built (after a call to
    * `buildBasis` or more specifically `initStates`, since the constants
    * depend on the initialization in <tt>initStates</tt>).
    */
   void initLowerBoundL2 (int dim1, int dim2);

   /**
    * These precomputed constants allows the calculation of the square
    * length of a lattice vector \f$\ell_2\f$ from a value of the merit
    * \f$S_2\f$ for each dimension \f$i\f$, i.e. \f$\ell_2[i] =
    * S_2[i]*\f$<tt>m_S2toL2[</tt>\f$i\f$<tt>]</tt>.
    */
   double *m_S2toL2;

   /**
    * Prepares and dispatches the results for dimension `dim` to all
    * observers attached to this test.
    */
   void prepAndDisp (int dim);

   /**
    * Sends the initialization message to all observers attached to this
    * test.
    */
   void init ();
};

}
#endif
