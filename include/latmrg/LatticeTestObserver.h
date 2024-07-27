#ifndef	LATTICETESTOBSERVER_H
#define LATTICETESTOBSERVER_H
#include "latcommon/Base.h"
#include <string>


namespace LatMRG {

/**
 * Interface that classes must implement in order to receive information and
 * results from a lattice test. See `LatticeTest`.
 *
 */
class LatticeTestObserver {
public:

/**
 * Destructor.
 */
virtual ~LatticeTestObserver() {}

   /**
    * Called when the base is incremented in a lattice test. `base` is a
    * copy of the base used in the test.
    */
   virtual void baseUpdate (LatCommon::Base & base) = 0;

   /**
    * Called when the base is incremented in a lattice test. `V[i]` is a
    * copy of basis vector \f$i\f$ used in the test.
    */
   virtual void baseUpdate (LatCommon::Base & V, int i) = 0;

   /**
    * Called when new results have been calculated for one dimension. The
    * results are placed in the array `results` of size `n`.
    */
   virtual void resultUpdate (double results[], int n) = 0;

   /**
    * Called when a test is initiated. `test` contains the name of the
    * test being performed, `headers` is the array of names of the values
    * calculated by this test and `n` is the size of the array `headers`.
    */
   virtual void testInit (const std::string & test, std::string headers[],
                          int n) = 0;

   /**
    * Called when the test has terminated successfully.
    */
   virtual void testCompleted() = 0;

   /**
    * Called when the test has terminated but failed. `dim` indicates in
    * which dimension the test failed.
    */
   virtual void testFailed (int dim) = 0;
};

}
#endif
