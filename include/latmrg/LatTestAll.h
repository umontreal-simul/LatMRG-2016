#ifndef LATTESTALL_H
#define LATTESTALL_H
#include "Writer.h"
#include "LatConfig.h"


namespace LatMRG {

/**
 * This class is just an auxiliary class that allows launching the spectral,
 * Beyer, or Palpha test on several different generators successively, each
 * one associated with its own data file. One can launch these tests on all
 * the data files (with extension <tt>".dat"</tt>) in a directory also. The
 * test and generator parameters in the data files must be as described in
 * program `LatMain` (see page (FIXME: page#) of this guide). In fact, the
 * `LatMain` program simply calls the methods of this class.
 *
 */
class LatTestAll {
public:

   /**
    * Reads the parameters of the test and of the generator in input text file
    * `datafile`; then do the test. The data file must always have the extension
    * `".dat"`, but must be given as argument here *without extension*. For
    * example, if the data file is named `mrg.dat`, then the method must be
    * called as `doTest("mrg")`. Returns 0 if the test completed successfully;
    * returns a negative integer if there was an error.
    */
   int doTest (const char *datafile);

   /**
    * Applies the method `doTest` to all the files with extension `".dat"`
    * in directory named `dirname`. Returns 0 if all the tests completed
    * successfully; returns a non-zero integer if there was an error.
    */
   int doTestDir (const char *dirname);

private:

   /**
    * Returns a `Writer` created from the input file `infile` and the given
    * `OutputType`.
    */
   Writer* createWriter (const char *infile, LatCommon::OutputType ot);
};

}
#endif
