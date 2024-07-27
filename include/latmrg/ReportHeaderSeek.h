#ifndef	REPORTHEADERSEEK_H
#define	REPORTHEADERSEEK_H
#include "ReportHeader.h"
#include "SeekConfig.h"
#include "Writer.h"


namespace LatMRG {

/**
 * This class is an implementation of the `ReportHeader` abstract class for
 * the programs `seek*`. It prints the configuration of the **** UNFINISHED
 * ****
 *
 */
class ReportHeaderSeek : public ReportHeader {
public:

/**
 * Constructor. `writer` is the writing engine used to write the report
 * header. `config` is the configuration ********* of the lattice test, which
 * is populated from an instance of `ParamReaderSeek`.
 */
ReportHeaderSeek (Writer *writer, SeekConfig *config);

   /**
    * Does the actual writing of the report header with the `Writer`
    * passed to the constructor.
    */
   void printHeader();
private:

/**
 * Pointer to the configuration of the seek program.
 */
SeekConfig* m_config;
};

}
#endif
