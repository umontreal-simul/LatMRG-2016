#include "latmrg/ProjIteratorSuccCoords.h"

#include <stdexcept>

namespace LatMRG
{

ProjIteratorSuccCoords::ProjIteratorSuccCoords (size_t minCoord,
      size_t maxCoord, size_t minOrder, size_t maxOrder)
      : ProjIteratorDefault (minCoord, maxCoord, minOrder, maxOrder)
{
   resetAtOrder (minOrder);
}


//===========================================================================

ProjIteratorSuccCoords::ProjIteratorSuccCoords (size_t minCoord,
      size_t maxCoord, size_t minOrder, size_t maxOrder,
      bool forceMinCoord, bool forceMaxCoord)
      : ProjIteratorDefault (minCoord, maxCoord, minOrder, maxOrder,
                             forceMinCoord, forceMaxCoord)
{
   if (forceMinCoord && forceMaxCoord)
      throw std::invalid_argument (
          "forceMinCoord and forceMaxCoord cannot be both true for successive coordinates");
   resetAtOrder (minOrder);
}


//===========================================================================

void ProjIteratorSuccCoords::resetAtOrder (size_t order)
{
   ProjIteratorDefault::resetAtOrder (order);
   if (m_forceMaxCoord) {
      m_projection.clear ();
      size_t val = m_maxCoord;
      while (m_projection.size () < m_order)
         m_projection.insert (val--);
   }
}


//===========================================================================

ProjIterator & ProjIteratorSuccCoords::operator++ ()
{
   // Increments all projection indices

   if (!m_forceMinCoord && !m_forceMaxCoord && m_order >= 1) {
      // increment all indices
      size_t val = *m_projection.begin ();
      m_projection.clear ();
      while (m_projection.size () < m_order)
         m_projection.insert (++val);
      if (*m_projection.rbegin () <= m_maxCoord)
         return *this;
   }
   // if we reach this point, we have walked all the projections for the
   // current order
   if (m_order < m_maxOrder)
      resetAtOrder (m_order + 1);
   else
      m_atEnd = true;

   return *this;
}


//===========================================================================
}
