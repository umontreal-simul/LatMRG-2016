#include "latmrg/ProjIteratorDefault.h"
#include <stdexcept>

using namespace LatCommon;


namespace LatMRG
{

ProjIteratorDefault::ProjIteratorDefault (size_t minCoord,
      size_t maxCoord, size_t minOrder, size_t maxOrder)
      : m_minCoord (minCoord), m_maxCoord (maxCoord), m_minOrder (minOrder),
      m_maxOrder (maxOrder), m_forceMinCoord (false), m_forceMaxCoord (false)
{
   reset ();
}


// ===========================================================================

ProjIteratorDefault::ProjIteratorDefault (size_t minCoord,
      size_t maxCoord, size_t minOrder, size_t maxOrder,
      bool forceMinCoord, bool forceMaxCoord)
      : m_minCoord (minCoord), m_maxCoord (maxCoord), m_minOrder (minOrder),
      m_maxOrder (maxOrder), m_forceMinCoord (forceMinCoord),
      m_forceMaxCoord (forceMaxCoord)
{
   if (forceMinCoord && forceMaxCoord && minOrder < 2)
      throw std::invalid_argument
      ("minOrder must be > 2 to use both forceMinCoord and forceMaxCoord");
   reset ();
}


// ===========================================================================

void ProjIteratorDefault::reset ()
{
   resetAtOrder (m_minOrder);
}


// ===========================================================================

void ProjIteratorDefault::resetAtOrder (size_t order)
{
   m_order = order;
   m_projection.clear ();

   // insert max coordinate if required
   if (m_forceMaxCoord && m_order >= 1)
      m_projection.insert (m_maxCoord);

   // insert other coordinates, starting from min coordinate
   size_t val = m_minCoord;
   while (m_projection.size () < m_order)
      m_projection.insert (val++);

   // check if any projections at all are available
   m_atEnd = m_order > m_maxOrder
             || m_order > (m_maxCoord - m_minCoord + 1);
}


// ===========================================================================

const Coordinates & ProjIteratorDefault::operator* () const
{
   if (m_atEnd)
      throw std::out_of_range ("ProjIteratorDefault: already at end");
   return m_projection;
}


// ===========================================================================

ProjIterator & ProjIteratorDefault::operator++ ()
{
   // Increments the projection indices as if they were a sequence of digits,
   // with the additional constraint that the digits must be in increasing
   // order.

   size_t maxval = m_maxCoord; // max index value

   // iterator before the min element
   Coordinates::reverse_iterator rend = m_projection.rend ();

   // iterator on the max element
   Coordinates::reverse_iterator rbegin = m_projection.rbegin ();

   if (m_forceMinCoord && !m_projection.empty ())
      // leave the minimum coordinate unchanged
      --rend;

   if (m_forceMaxCoord && !m_projection.empty ()) {
      // leave the maximum coordinate unchanged
      ++rbegin;
      --maxval;
   }

   Coordinates::reverse_iterator rit = rbegin;

   // find the maximum index that has not reached its maximum value yet
   while (rit != rend && *rit >= maxval) {
      ++rit;
      --maxval;
   }

   if (rit == rend) {
      // all indices have reached their maximum value

      if (m_order < m_maxOrder)
         // go to next order
         resetAtOrder (m_order + 1);
      else
         // all projections exhausted
         m_atEnd = true;

      return *this;
   }
   // make note of the highest index that is below its maximum value (to be
   // incremented)
   size_t high = *rit;

   // erase all indices that are at their maximum value plus the one to be
   // incremented
   m_projection.erase (--rit.base (), rbegin.base ());

   // insert new values
   while (m_projection.size () < m_order)
      m_projection.insert (++high);

   return *this;
}


// ===========================================================================
}
