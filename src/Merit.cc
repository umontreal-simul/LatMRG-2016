#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cfloat>

#include "latmrg/Merit.h"
#include "latcommon/Util.h"

using namespace std;
using namespace LatCommon;

namespace LatMRG
{

Merit::Merit (int maxDim)
{
   m_val.reserve(1 + maxDim);
   m_normVal.reserve(1 + maxDim);
   m_worstMerit = DBL_MAX;   // Arbitrary large value
   m_dimWorst = 0;
   set (-1.0);
}


//=========================================================================

Merit::Merit (const Merit & mer) :
      m_worstMerit (mer.m_worstMerit), m_dimWorst (mer.m_dimWorst)
{
   m_val.clear();
   m_normVal.clear();
   int dim = mer.getDim();
   m_val.reserve (dim + 1);
   m_normVal.reserve (dim + 1);
   set (-1.0);
   /*
   for (int i = 0; i <= dim; ++i) {
      m_val[i] = mer.m_val[i];
      m_normVal[i] = mer.m_normVal[i];
   }
   */
}


//=========================================================================

void Merit::set (double x)
{
   for (unsigned int i = 0; i < m_val.capacity(); ++i) {
      m_val[i] = x;
      m_normVal[i] = x;
   }
}


//=========================================================================

Merit::~Merit ()
{
   m_val.clear();
   m_normVal.clear();
}


//=========================================================================

Merit & Merit::operator= (const Merit & mer)
{
   if (this != &mer) {
      m_worstMerit = mer.m_worstMerit;
      m_dimWorst = mer.m_dimWorst;
      int dim = mer.getDim();
      if (getDim() < dim) {
         m_val.clear();
         m_normVal.clear();
         m_val.reserve (dim + 1);
         m_normVal.reserve (dim + 1);
      }
      set (-1.0);
      /*
      for (int i = 0; i <= dim; ++i) {
         m_val[i] = mer.m_val[i];
         m_normVal[i] = mer.m_normVal[i];
      }
      */
   }
   return *this;
}


//=======================================================================

double Merit::getST (int fromDim, int T, int & dimWorst)
{
   double min = 1.0e100;
   m_dimWorst = -1;
   for (int i = fromDim; i <= T; i++) {
      if (m_normVal[i] < min) {
         min = m_normVal[i];
         m_dimWorst = i;
      }
   }
   dimWorst = m_dimWorst;
   m_worstMerit = min;
   return min;
}


//=======================================================================

string Merit::toString (int from, int to, bool rac, bool invert) const
{
   double x, y;
   ostringstream os;
   const char *espace = "             ";

   if (rac) {
      for (int i = from; i <= to; ++i) {
         x = mysqrt(m_val[i]);
         y = mysqrt(m_normVal[i]);
         if (invert)
            x = 1.0/x;
         os << setw(5) << right << i << espace << setw(12) << x
            << espace << setw(10) << left << y << "\n";
      }

   } else {
      for (int i = from; i <= to; ++i) {
         x = m_val[i];
         if (invert)
            x = 1.0/x;
         os << setw(5) << right << i << espace << setw(12) << x
            << espace << setw(10) << left << m_normVal[i] << "\n";
      }
   }
//   os << "\n";
   return os.str ();
}


//=======================================================================

}   // namespace LatMRG
