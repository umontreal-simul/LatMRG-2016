#ifndef LATTESTBEYER_H
#define LATTESTBEYER_H
#include "latmrg/IntLattice.h"
#include "LatticeTest.h"


namespace LatMRG {

/**
 * This class implements the *Beyer* test. It implements the abstract class
 * `LatticeTest`. The figure of merit for this test is the Beyer quotient.
 * The main program is obtained by compiling the `LatMain.cc` file (see the
 * description of module `LatMain`, on page (FIXME: page#), for how to run a
 * program).
 * \remark **Richard:** Dixit Pierre: il est douteux que cette classe devrait
 * exister: ce devrait être une méthode de Lattice ou quelque chose du genre.
 * Idem pour les autres `LatTest*`
 *
 */
class LatTestBeyer : public LatticeTest {
public:

   /**
    * Constructor. The *Beyer* test will be applied on lattice `lat`.
    */
   LatTestBeyer (LatMRG::IntLattice * lat);

   /**
    * Destructor.
    */
   ~LatTestBeyer () {};

   /**
    * Applies the *Beyer* test for dimensions varying from
    * <tt>fromDim</tt> to `toDim`. Whenever the normalized value of the
    * merit is smaller than `minVal` for any dimension, the method returns
    * `false` immediately. The method returns `false` if the test was
    * interrupted for any reason before completion, and it returns `true`
    * upon success. The results of the last test are kept in `merit`.
    */
   bool test (int fromDim, int toDim, double minVal[]);

private:

   /**
    * Prepares and dispatches the results for dimension `dim` to all observers
    * attached to this test.
    */
   void prepAndDisp (int dim);

   /**
    * Sends the initialization message to all observers attached to this
    * test.
    */
   void init ();
};

}
#endif // LATTESTBEYER_H
