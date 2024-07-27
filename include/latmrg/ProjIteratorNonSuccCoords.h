#ifndef PROJ_ITERATOR_NON_SUCC_COORDS_H
#define PROJ_ITERATOR_NON_SUCC_COORDS_H
#include "ProjIteratorDefault.h"


namespace LatMRG {

/**
 * This projection iterator walks through all projections composed of
 * non-successive coordinate indices.
 *
 */
class ProjIteratorNonSuccCoords : public ProjIteratorDefault {
public:

   /**
    * Constructor. Creates a projection iterator for all projections of orders
    * `minOrder` to `maxOrder` with coordinate indices ranging from `minCoord`
    * to `maxCoord`.
    */
   ProjIteratorNonSuccCoords (size_t minCoord, size_t maxCoord, size_t minOrder,
                                 size_t maxOrder);

   /**
    * Constructor. Creates a projection iterator for all projections of
    * orders `minOrder` to `maxOrder` with coordinate indices ranging from
    * `minCoord` to `maxCoord`. Also forces the first coordinate to
    * `minCoord` if `forceMinCoord` is `true`, and forces the last coordinate
    * to `maxCoord` if `forceMaxCoord` is `true`.
    */
   ProjIteratorNonSuccCoords (size_t minCoord, size_t maxCoord, size_t minOrder,
                              size_t maxOrder, bool forceMinCoord, bool forceMaxCoord);

   /**
    * Resets the iterator to the first projection of order `order`.
    */
   virtual void resetAtOrder (size_t order);

   /**
    * Advances to the next projection.
    */
   virtual ProjIterator & operator++();
protected:

   static bool successive (const LatCommon::Coordinates & indices);
};

}
#endif
