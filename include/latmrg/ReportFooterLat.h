#ifndef REPORTFOOTERLAT_H
#define REPORTFOOTERLAT_H
#include "ReportFooter.h"
#include "Writer.h"
#include "LatticeTest.h"


namespace LatMRG {

/**
 * This class is an implementation of the `ReportFooter` abstract class for
 * the program `lat*`.
 *
 */
class ReportFooterLat : public ReportFooter {
public:

/**
 * Constructor. `writer` is the writing engine used to write the report
 * footer. `test` is the lattice test thas was performed and for which the
 * results are to be written.
 */
ReportFooterLat (Writer *, LatticeTest * test = 0);

   /**
    * `test` is the lattice test thas was performed and for which the
    * results are to be written.
    */
   void setLatticeTest (LatticeTest * test) { m_test = test; }

   /**
    * Defined in abstract class `ReportFooter`.
    */
   void printFooter();
private:

/**
 * Pointer to the final lattice test which was performed.
 */
LatticeTest * m_test;
};

}
#endif
