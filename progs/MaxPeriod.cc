#include "latcommon/Types.h"
#include "latcommon/Const.h"
#include "latcommon/Util.h"
#include "latmrg/MRGComponent.h"
#include "latmrg/ParamReader.h"

#include <string>

using namespace NTL;
using namespace std;
using namespace LatMRG;
using namespace LatCommon;



//===========================================================================

namespace
{

const int MAX_CHARS = 64;

typedef struct {
   GenType typ;
   MScal m;
   int k;
   DecompType decom1;
   string filem1;
   DecompType decor;
   string filer;
   MVect a;
} ParamType;


//===========================================================================

void read (ParamReader & reader, ParamType & par)
{
   reader.getLines ();
   int ln = 0;
   reader.readGenType (par.typ, ++ln, 1);
   long b, e, c;
   reader.readNumber3 (par.m, b, e, c, ++ln, 1);
   reader.readInt (par.k, ++ln, 1);
   reader.readDecompType (par.decom1, ++ln, 1);
   if (par.decom1 == DECOMP_WRITE || par.decom1 == DECOMP_READ) {
      par.filem1.reserve (MAX_CHARS);
      reader.readString (par.filem1, ln, 2);
   }
   reader.readDecompType (par.decor, ++ln, 1);
   if (par.decor == DECOMP_WRITE || par.decor == DECOMP_READ) {
      par.filer.reserve (MAX_CHARS);
      reader.readString (par.filer, ln, 2);
   }
   par.a.SetLength (1 + par.k);
   for (int i = 1; i <= par.k; i++)
      reader.readMScal(par.a[i], ++ln, 1);
   par.a[0] = 0;
}

}  // namespace


//===========================================================================

int main(int argc, char** argv)
{
   if (argc != 2) {
      cout << "Usage: " << argv[0] << " <data file>" << endl;
      return 1;
   }
   string fname(argv[1]);
   fname += ".dat";
   ParamReader reader (fname);
   ParamType par;
   read (reader, par);

   MRGComponent mrg (par.m, par.k, par.decom1, par.filem1.c_str(),
                                   par.decor,  par.filer.c_str());

   cout << "   \nThe generator with" << endl;
   cout << "   m = " << par.m << endl;
   cout << "   k = " << par.k << endl;
   cout << "   a = ";
   cout <<  toString(par.a, 1, par.k);

   if (mrg.maxPeriod (par.a))
      cout << "      HAS maximal period." << endl;
   else
      cout << "      does NOT have maximal period." << endl;

   return 0;
}
