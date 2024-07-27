/*
 MRGComponentFactory.cc for ISO C++
 version 1.00

 authors: Hicham Wahbi
          Frédérik Rozon
          Richard Simard
*/

#include "latmrg/MRGComponentFactory.h"
#include "latcommon/Util.h"

using namespace LatCommon;

namespace LatMRG
{

MRGComponent* MRGComponentFactory::fromMWC(const MScal & b, const MVect & a, int k)
{
   MScal _m;
   MScal _b;
   MVect _a;

   MScal d, e, f, g;

   conv(_m, 0);
   conv(_b, 1);
   CreateVect(_a, 1);

   //Calcul de m
   for (int i = 0; i <= k; i++) {
      _m += a[i] * _b;
      _b *= b;
   }

   //Calcul de a
   Euclide(b, _m, _a[1], d, e, f, g);

   if (g < 0)
      _a[1] = -_a[1];

   MRGComponent* lat = new MRGComponent(_m, _a, 1);
   DeleteVect(_a);
   return lat;
}

}
