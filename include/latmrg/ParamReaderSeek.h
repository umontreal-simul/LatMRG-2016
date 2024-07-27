#ifndef PARAMREADERSEEK_H
#define PARAMREADERSEEK_H
#include "ParamReader.h"
#include "SeekConfig.h"

#include <string>


namespace LatMRG {

/**
 * This class is used to read a configuration file for the executable
 * programs `seek*`, created from the `SeekMain` main program. The format of
 * the configuration file is described in this guide for the program
 * `SeekMain` on page (FIXME: page#).
 *
 */
class ParamReaderSeek : public ParamReader {
public:

/**
 * Default Constructor.
 */
ParamReaderSeek();

   /**
    * Constructor. `fname` is the name of the configuration file.
    */
   ParamReaderSeek (std::string fname);

   /**
    * Destructor.
    */
   ~ParamReaderSeek();

   /**
    * Variables initialization.
    */
   void init();

   /**
    * Reads the configuration parameters and puts them in an instance of
    * <tt>SeekConfig</tt>.
    */
   void read (SeekConfig & config);
};

}
#endif
