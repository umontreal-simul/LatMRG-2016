#include "latmrg/Primes.h"
#include "latmrg/ParamReader.h"

#include <string>

using namespace std;
using namespace LatMRG;



//===========================================================================

namespace
{

typedef struct {
   int k;
   int e;
   long c1;
   long c2;
   bool safe;
   bool fac;
} ParamType;


//===========================================================================

void read (ParamReader & reader, ParamType & par)
{
   reader.getLines ();
   int ln = 0;
   reader.readInt (par.k, ++ln, 1);
   reader.readInt (par.e, ++ln, 1);
   reader.readLong (par.c1, ++ln, 1);
   reader.readLong (par.c2, ++ln, 1);
   reader.readBool (par.safe, ++ln, 1);
   reader.readBool (par.fac, ++ln, 1);
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

   string oname(argv[1]);
   oname += ".res";

   ofstream fout (oname.c_str());
   Primes primes;
   primes.find (par.k, par.e, par.c1, par.c2, par.safe, par.fac, fout);

   return 0;
}
