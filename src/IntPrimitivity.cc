#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <cassert>

#include "latmrg/IntPrimitivity.h"
#include "latcommon/Util.h"
#include "NTL/ZZ.h"

using namespace LatCommon;

namespace LatMRG
{

IntPrimitivity::IntPrimitivity () : m_e(1)
{
   m_p = 0;
   m_m = 0;
}


//===========================================================================

void IntPrimitivity::setpe (const MScal & p, long e)
{
   m_p = p;
   m_e = e;
   m_m = power (m_p, m_e);
}

//===========================================================================

IntPrimitivity::IntPrimitivity (const IntFactorization & f,
                                const MScal & p, long e) : m_f(f)
{
   setpe(p, e);
}

//===========================================================================

std::string IntPrimitivity::toString () const
{
   std::ostringstream out;
   out << "p = " << m_p << std::endl;
   out << "e = " << m_e << std::endl;
   out << "m = " << m_m << std::endl;
   out << "\nf is the factorization of  " << m_f.toString () << std::endl;
   return out.str ();
}

//===========================================================================

bool IntPrimitivity::isPrimitiveElement (const MScal & a) const
throw(std::range_error)
{
   if (0 == m_p)
      throw std::range_error("IntPrimitivity::isPrimitiveElement:   p = 0");
   if (0 == a)
      return false;

   MScal t1, t2;
   t1 = a;
   if (t1 < 0)
      t1 += m_m;

   const std::vector<MScal> invList = m_f.getInvFactorList();
//   assert (!(invList.empty ()));
   std::vector<MScal>::const_iterator it = invList.begin();
   while (it != invList.end()) {
      t2 = PowerMod (t1, *it, m_m);
      if (t2 == 1)
         return false;
      ++it;
   }
   return true;
}


//===========================================================================

bool IntPrimitivity::isPrimitiveElement (const MVect & V, int k) const
throw(std::range_error)
{
   MScal a;
   if (k & 1)
      a = V[k];
   else
      a = -V[k];
   return isPrimitiveElement (a);
}


//===========================================================================

}
