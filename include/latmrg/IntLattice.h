#ifndef LATMRG__INTLATTICE_H
#define LATMRG__INTLATTICE_H
#include "latcommon/IntLattice.h"
#include "latmrg/MRGComponent.h"

namespace LatMRG {

/**
 * \copydoc LatCommon::IntLattice
 *
 */
class IntLattice : public LatCommon::IntLattice {
public:

   /**
    * \copydoc LatCommon::IntLattice::IntLattice(const MScal&, int, int, NormType)
    */
   IntLattice (const MScal & m, int k, int maxDim, LatCommon::NormType norm = LatCommon::L2NORM) :
      LatCommon::IntLattice(m, k, maxDim, norm) {}

   /**
    * \copydoc LatCommon::IntLattice::IntLattice(const IntLattice&)
    */
   IntLattice (const IntLattice & Lat) :
      LatCommon::IntLattice(Lat) {}

   /**
    * Destructor.
    */
   virtual ~IntLattice ();

#ifdef WITH_NTL
   /**
    * The components of the lattice when it is built out of more than one
    * component. When there is only one component, it is unused as the
    * parameters are the same as above.
    */
   std::vector<MRGComponent *> comp;
#endif


protected:

   /**
    * \copydoc LatCommon::IntLattice::kill()
    */
   virtual void kill ();

};

}
#endif
