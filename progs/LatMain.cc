#include <sys/types.h>
#include <sys/stat.h>
#include <iostream>
#include "latmrg/LatTestAll.h"


using namespace std;
using namespace LatMRG;


//==========================================================================

int main (int argc, char **argv)
{
   if (argc < 2) {
      cerr << "\n*** Usage:\n   "
           << argv[0] << " data_file1 data_file2 ...." << endl
           << "or\n   "
           << argv[0] << " dir1 dir2 ...." << endl
           << endl;
      return -1;
   }

   struct stat buf;    // properties of a file or directory
   LatTestAll testall;
   int status = 0;

   for (int j = 1; j < argc; j++) {
      // Do the test for each data file or directory on the command line
      stat(argv[j], &buf);
      if (0 != S_ISDIR(buf.st_mode))         // directory
         status |= testall.doTestDir (argv[j]);
      else {
         string dataname(argv[j]);
         dataname.append(".dat");
         stat(dataname.c_str(), &buf);
         if (0 != S_ISREG(buf.st_mode))    // data file
            status |= testall.doTest (argv[j]);
      }
   }

   return status;
}
