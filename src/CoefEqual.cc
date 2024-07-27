#include "latmrg/CoefEqual.h"


namespace LatMRG {

CoefEqual::CoefEqual(int *I, int s): Coefficient(), m_I(I), m_s(s)
{
   I[0] = 0;
}


//===========================================================================

CoefEqual::~CoefEqual()
{
}


//===========================================================================

void CoefEqual::set (const MScal & q, MVect & A, int & i)
{
   int r = m_s;
   while ((r > 0) && (i <= m_I[r]))
      --r;
   int j = m_I[r];
   for (r = i; r > j; r--)
      A[r] = q;
   i = j + 1;
}

}
