#ifndef LATTESTPALPHA_H
#define LATTESTPALPHA_H
#include "LatConfig.h"
#include "LatticeTest.h"
#include "latcommon/Normalizer.h"


namespace LatMRG {

/**
 * This class applies the \f$P_{\alpha}\f$ test by implementing the abstract
 * class `LatticeTest`. The figure of merit is the
 * \f$P_{\alpha}\f$ criterion (see class `Palpha` on page (FIXME: page#)).
 * The main program is obtained by compiling the `LatMain.cc` file, and the
 * executable file is called `latLLDD`, since the test is implemented only
 * for a small number of points (\f$m < 2^{31}\f$). See the description of
 * the `LatMain` program on page (FIXME: page#).
 *
 */
class LatTestPalpha : public LatticeTest {
public:

   /**
    * Constructor. The \f$P_{\alpha}\f$ test will be applied on lattice whose
    * parameters are in `config`. The `bounds` \f$B_{\alpha}(s)\f$ may be used
    * to normalize the \f$P_{\alpha}(s)\f$ values.
    */
   LatTestPalpha (LatCommon::Normalizer * bounds, LatMRG::IntLattice * lat);

   /**
    * Destructor.
    */
   ~LatTestPalpha () {};
   void setConfig (LatConfig * config) {m_config = config;};

   /**
    * Applies the \f$P_{\alpha}\f$ test for dimensions varying from
    * <tt>fromDim</tt> to `toDim`. Whenever the normalized value of the
    * merit is smaller than `minVal` for any dimension, the method returns
    * `false` immediately. The method returns `false` if the test was
    * interrupted for any reason before completion, and it returns `true`
    * upon success. The results of the last test are kept in `merit`.
    */
   bool test (int fromDim, int toDim, double minVal[]);

private:

   /**
    * Contains the parameters of the test.
    */
   LatConfig *m_config;

   /**
    * The \f$B_{\alpha}\f$ bounds used to normalize the
    * \f$P_{\alpha}\f$.
    */
   LatCommon::Normalizer *m_bound;

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
