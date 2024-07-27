#include "latmrg/CoefZero.h"


namespace LatMRG {

CoefZero::CoefZero(int *I, int s): Coefficient(), m_I(I), m_s(s)
{
   I[0] = 0;
}


//===========================================================================

CoefZero::~CoefZero()
{
}


//===========================================================================

void CoefZero::set (const MScal & q, MVect & A, int & i)
{
   A[i] = q;
   int r = m_s;
   while ((r > 0) && (i <= m_I[r]))
      --r;
   i = m_I[r] + 1;
}

}
