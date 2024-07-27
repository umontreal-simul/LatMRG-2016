#ifndef LATTICETEST_H
#define LATTICETEST_H
#include "latcommon/Types.h"
#include "latcommon/Util.h"
#include "latcommon/Const.h"
#include "latcommon/Base.h"
#include "latmrg/IntLattice.h"
#include "latmrg/Merit.h"
#include "latmrg/Subject.h"
#include "latmrg/LatticeTestObserver.h"
#include "latcommon/Weights.h"
#include "latmrg/Chrono.h"
#include <string>
#include <list>


namespace LatMRG {

/**
 * This class is the base class that lattice tests must implement. These
 * tests are applied on lattices to assess their structural properties and
 * their qualities with respect to different criteria. Included are
 * well-known tests such as the *spectral* test, the *Beyer* test, the
 * \f$P_{\alpha}\f$ test. The corresponding figures of merit for the lattice
 * are the length of the shortest vector in the *primal* or in the *dual*
 * lattice computed with different norms, the Beyer quotient, or the
 * \f$P_{\alpha}\f$ criterion. For the standard spectral test, the figure of
 * merit is based on the length of the shortest non-zero vector in the *dual*
 * lattice, using the \f${\mathcal{L}}_2\f$ norm to compute the length of
 * vectors, and the inverse of this length gives the maximal distance between
 * successive hyperplanes covering all the points in the *primal* lattice. If
 * one computes the length of the shortest non-zero vector in the *dual*
 * lattice using the \f${\mathcal{L}}_1\f$ norm, one obtains the minimal
 * number of hyperplanes covering all the points of the *primal* lattice.
 *
 */
class LatticeTest: public Subject<LatticeTestObserver *> {
public:

   /**
    * Constructor. The test will be applied on `lattice`.
    */
   explicit LatticeTest (LatMRG::IntLattice * lattice);

   /**
    * Destructor.
    */
   virtual ~LatticeTest ();

   /**
    * If `dualF` is `true`, the tests will be applied on the *dual*
    * lattice; if it is `false`, the tests will be applied on the *primal*
    * lattice.
    */
   void setDualFlag (bool dualF);

   /**
    * Gets the value of the <tt>m_dualF</tt> flag.
    */
   bool getDualFlag () const { return m_dualF; }

   /**
    * Sets the value of the <tt>m_invertF</tt> flag. If `invertF` is
    * `true`, the inverse of the length of the shortest vector will be
    * printed in the results. Otherwise, the length itself will be
    * printed. By default, the value is set to <tt>false</tt>.
    */
   void setInvertFlag (bool invertF);

   /**
    * Gets the value of the <tt>m_invertF</tt> flag.
    */
   bool getInvertFlag () const { return m_invertF; }

   /**
    * When `detail` \f$>0\f$, this flag indicates to print extra detailed
    * results. Default value: 0.
    */
   void setDetailFlag (int detail);

   /**
    * Gets the value of the <tt>m_detailF</tt> flag.
    */
   int getDetailFlag () const { return m_detailF; }

   /**
    * If `maxAllDimF` is `true`, the merit will be maximized in all
    * dimensions.
    */
   void setMaxAllDimFlag (bool maxAllDimF);

   /**
    * Sets the maximum number of nodes in the branch-and-bound tree to
    * `maxNodesBB`. Default value: \f$10^7\f$.
    */
   void setMaxNodesBB (long maxNodesBB);

   /**
    * Gets the results of the last applied test.
    */
   Merit & getMerit () { return m_merit; }

   /**
    * Returns the criterion used for the test.
    */
   LatCommon::CriterionType getCriterion() const { return m_criter; }

   /**
    * Gets the lattice on which the test is applied.
    */
   LatMRG::IntLattice * getLattice () const { return m_lat; }

   /**
    * Returns the lowest dimension on which the test is performed.
    */
   int getMinDim () const { return m_fromDim; }

   /**
    * Returns the highest dimension on which the test is performed.
    */
   int getMaxDim () const { return m_toDim; }

   /**
    * Ensures that `fromDim` is larger than `order`.
    */
   void resetFromDim (int order, int & fromDim);

   /**
    * Starts the test from dimension `minDim` to dimension `maxDim`.
    * Whenever the normalized value of the merit is smaller than `minVal`
    * for any dimension, the method returns `false` immediately. The
    * method returns `false` if the test was interrupted for any reason
    * before completion, and it returns `true` upon success. The results
    * of the test are kept in <tt>m_merit</tt>.
    */
   virtual bool test (int minDim, int maxDim, double minVal[]) = 0;

   /**
    * Similar to `test` above, but with the weights `weights`.
    */
   virtual bool test (int minDim, int maxDim, double minVal[],
                      const double* weights);
protected:

   /**
    * Timer which measures CPU time as test goes on.
    */
   Chrono timer;

   /**
    * Contains the results of the last test.
    */
   Merit m_merit;

   /**
    * The criterion used to evaluate the figure of merit.
    */
   LatCommon::CriterionType m_criter;

   /**
    * If `true`, the test is applied on the dual lattice, otherwise on the
    * primal lattice.
    */
   bool m_dualF;

   /**
    * If `true`, the inverse of the length of the shortest vector is
    * printed in the results. Otherwise, the length is printed.
    */
   bool m_invertF;

   /**
    * When <tt>m_detailF</tt> \f$>0\f$, this flag indicates to print extra
    * detailed results. Default value: 0.
    */
   int m_detailF;

   /**
    * If `true`, the merit is to be maximized in all dimensions.
    */
   bool m_maxAllDimFlag;

   /**
    * The maximum number of nodes in the branch-and-bound tree.
    */
   long m_maxNodesBB;

   /**
    * The lattice on which the test is applied.
    */
   LatMRG::IntLattice* m_lat;

   /**
    * The dimension parameters.
    */
   int m_fromDim, m_toDim;

   /**
    * The list that contains the observers for the lattice test.
    */
   typedef std::list<LatticeTestObserver*> ob_list;

   /**
    * Dispatches a `baseUpdate` signal to all observers.
    */
   void dispatchBaseUpdate (LatCommon::Base &);

   /**
    * Dispatches a `baseUpdate(V, i)` signal to all observers. Only base
    * vector \f$i\f$ will be sent.
    */
   void dispatchBaseUpdate (LatCommon::Base & V, int i);

   /**
    * Dispatches a `resultUpdate` signal to all observers.
    */
   void dispatchResultUpdate (double[], int);

   /**
    * Dispatches a `testInit` signal to all observers.
    */
   void dispatchTestInit (const std::string &, std::string[], int);

   /**
    * Dispatches a `testCompleted` signal to all observers.
    */
   void dispatchTestCompleted ();

   /**
    * Dispatches a `testFailed` signal to all observers.
    */
   void dispatchTestFailed (int);
};

}     // namespace LatMRG
#endif
