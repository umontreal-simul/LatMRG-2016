#include "latmrg/IntLattice.h"

#include "NTL/quad_float.h"
#include "NTL/RR.h"
using namespace NTL;

using namespace std;

namespace LatMRG
{


//=========================================================================

void IntLattice::kill ()
{
   LatCommon::IntLattice::kill();
   if (!comp.empty()) {
      for (int s = 0; s < (int) comp.size(); s++)
         delete comp[s];
      comp.clear();
   }
}


//=========================================================================

IntLattice::~IntLattice ()
{
   kill ();
}

}
