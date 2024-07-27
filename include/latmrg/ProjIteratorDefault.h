#ifndef PROJ_ITERATOR_DEFAULT_H
#define PROJ_ITERATOR_DEFAULT_H
#include "ProjIterator.h"


namespace LatMRG {

/**
 * This projection iterator walks through all projections of a given range of
 * orders and coordinate indices.
 *
 */
class ProjIteratorDefault : public ProjIterator {
public:

   /**
    * Constructor. Creates a projection iterator for all projections of orders
    * `minOrder` to `maxOrder` with coordinate indices ranging from `minCoord`
    * to `maxCoord`.
    */
   ProjIteratorDefault (size_t minCoord, size_t maxCoord, size_t minOrder, size_t maxOrder);

   /**
    * Constructor. Creates a projection iterator for all projections of
    * orders `minOrder` to `maxOrder` with coordinate indices ranging from
    * `minCoord` to `maxCoord`. Also forces the first coordinate to
    * `minCoord` if `forceMinCoord` is `true`, and forces the last coordinate
    * to `maxCoord` if `forceMaxCoord` is `true`.
    */
   ProjIteratorDefault (size_t minCoord, size_t maxCoord, size_t minOrder, size_t maxOrder,
                        bool forceMinCoord, bool forceMaxCoord);

   /**
    * Resets the iterator to the first projection of minimal order.
    */
   virtual void reset();

   /**
    * Resets the iterator to the first projection of order `order`.
    */
   virtual void resetAtOrder (size_t order);

   /**
    * Returns the current projection.
    */
   virtual const LatCommon::Coordinates & operator*() const;

   /**
    * Advance to the next projection.
    */
   virtual ProjIterator & operator++();

   /**
    * Returns `true` if the current projection is valid, or `false` if all
    * projections have already been visited.
    */
   virtual operator bool() const { return !m_atEnd; }

protected:

   size_t m_order;

   size_t m_minCoord;
   size_t m_maxCoord;

   size_t m_minOrder;
   size_t m_maxOrder;

   bool m_atEnd;

   bool m_forceMinCoord;     // force first coordinate to the minimum index
   bool m_forceMaxCoord;     // force last coordinate to the maximum index

   LatCommon::Coordinates m_projection;
};

}
#endif
