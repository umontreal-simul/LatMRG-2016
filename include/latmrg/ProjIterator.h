#ifndef PROJ_ITERATOR_H
#define PROJ_ITERATOR_H
#include "latcommon/Weights.h"
#include "latcommon/CoordinateSets.h"


namespace LatMRG {

/**
 * This abstract class is the basis for different kinds of projection
 * iterators used to walk through sets of projections.
 *
 */
class ProjIterator {
public:

   /**
    * Destructor.
    */
   virtual ~ProjIterator() {}

   /**
    * Resets the iterator to the first projection.
    */
   virtual void reset() = 0;
   virtual const LatCommon::Coordinates * operator->() const {
   return &(operator*());
    }

   /**
    * Returns the current projection.
    */
   virtual const LatCommon::Coordinates & operator*() const = 0;

   /**
    * Advances to the next projection.
    */
   virtual ProjIterator & operator++() = 0;
   bool atEnd() const { return operator bool(); }

   /**
    * Returns `true` if the current projection is valid, or `false` if all
    * projections have already been visited.
    */
   virtual operator bool() const = 0;
};

}
#endif
