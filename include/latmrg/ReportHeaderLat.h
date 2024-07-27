#ifndef	REPORTHEADERLAT_H
#define	REPORTHEADERLAT_H
#include "ReportHeader.h"
#include "LatConfig.h"
#include "latmrg/IntLattice.h"
#include "Writer.h"


namespace LatMRG {

/**
 * This class is an implementation of the `ReportHeader` abstract class for
 * the programs `lat*`. It prints the configuration of the test launched, the
 * MRG or MWC components, and the combined MRG if applicable.
 *
 */
class ReportHeaderLat : public ReportHeader {
public:

   /**
    * Constructor. `writer` is the writing engine used to write the report
    * header. `config` is the configuration of the lattice test to be performed,
    * which is populated from an instance of `ParamReaderLat`. `lattice` is the
    * final MRG lattice on which the lattice test will be performed. It can be
    * the result of a combination of MRG components, a MWC transformed into a
    * MRG, and so on.
    */
   ReportHeaderLat (Writer *writer, LatConfig *config, IntLattice *lattice);

   /**
    * Does the actual writing of the report header with the `Writer`
    * passed to the constructor.
    */
   void printHeader();
private:

   /**
    * Pointer to the configuration of the lattice test.
    */
   LatConfig* m_config;

   /**
    * Pointer to the final MRG lattice on which the lattice test will be
    * performed.
    */
    IntLattice* m_lattice;
};

}
#endif
