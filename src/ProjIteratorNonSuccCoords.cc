#include "latmrg/ProjIteratorNonSuccCoords.h"
#include <stdexcept>

using namespace LatCommon;


namespace LatMRG
{

ProjIteratorNonSuccCoords::ProjIteratorNonSuccCoords (size_t minCoord,
      size_t maxCoord, size_t minOrder, size_t maxOrder)
      : ProjIteratorDefault (minCoord, maxCoord, minOrder, maxOrder)
{
   if (minOrder < 2)
      throw std::invalid_argument ("minOrder must be >= 2");
   resetAtOrder (minOrder);
}


//===========================================================================

ProjIteratorNonSuccCoords::ProjIteratorNonSuccCoords (size_t minCoord,
      size_t maxCoord, size_t minOrder, size_t maxOrder,
      bool forceMinCoord, bool forceMaxCoord)
      : ProjIteratorDefault (minCoord, maxCoord, minOrder, maxOrder,
                             forceMinCoord, forceMaxCoord)
{
   if (minOrder < 2)
      throw std::invalid_argument ("minOrder must be >= 2");
   resetAtOrder (minOrder);
}


//===========================================================================

void ProjIteratorNonSuccCoords::resetAtOrder (size_t order)
{
   ProjIteratorDefault::resetAtOrder (order);
   while (successive (m_projection) && !m_atEnd)
      operator++ ();
}


//===========================================================================

ProjIterator & ProjIteratorNonSuccCoords::operator++ ()
{
   ProjIteratorDefault::operator++ ();
   while (successive (m_projection) && !m_atEnd)
      operator++ ();
   return *this;
}

//===========================================================================

bool ProjIteratorNonSuccCoords::successive (const Coordinates & proj)
{
   Coordinates::iterator it1 = proj.begin ();
   Coordinates::iterator it2 = it1;
   while (++it2 != proj.end ()) {
      if (*it2 != *it1 + 1)
         return false;
      ++it1;
   }
   return true;
}


//===========================================================================
}
