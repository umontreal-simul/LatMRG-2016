#ifndef REPORTLAT_H
#define REPORTLAT_H
#include "LatticeTestObserver.h"
#include "ReportHeader.h"
#include "ReportFooter.h"
#include "DoubleFormatter.h"
#include "FormatterImpl.h"
#include "Writer.h"
#include "LatConfig.h"

#include <string>


namespace LatMRG {

/**
 * This class formats and prints the actual report for a lattice test for the
 * program `lat*`. It implements the interface `LatticeTestObserver` to be
 * able to receive information and results from the lattice test.
 *
 */
class ReportLat : public LatticeTestObserver {
public:

/**
 * Constructor.
 */
ReportLat (Writer* writer, LatConfig* config, ReportHeader* header,
              ReportFooter* footer);

   /**
    * Destructor.
    */
   ~ReportLat();

   /**
    * Prints the header using the `ReportHeader` passed to the
    * constructor.
    */
   void printHeader();

   /**
    * Prints the footer using the `ReportFooter` passed to the
    * constructor.
    */
   void printFooter();

   /**
    * Prints the table of results obtained from the successive calls of
    * `resultUpdate`. If more than one tests are performed, the results
    * will be concatenated in the same table.
    */
   void printTable();

   /**
    * Defined in interface `LatticeTestObserver`. Prints the base directly
    * in the report.
    */
   void baseUpdate (LatCommon::Base &);

   /**
    * Defined in interface `LatticeTestObserver`. Prints basis vector
    * \f$V[i]\f$ directly in the report.
    */
   void baseUpdate (LatCommon::Base & V, int i);

   /**
    * Defined in interface `LatticeTestObserver`. The results are stacked
    * in the internal table and will be printed upon a call to
    * `printTable`.
    */
   void resultUpdate (double[], int);

   /**
    * Defined in interface `LatticeTestObserver`. The columns in the
    * internal table are set up to be able to receive results from calls
    * to `resultUpdate`.
    */
   void testInit (const std::string &, std::string[], int);

   /**
    * Defined in interface `LatticeTestObserver`. Indicates a successful
    * test.
    */
   void testCompleted();

   /**
    * Defined in interface `LatticeTestObserver`. Indicates a failed test.
    * An error message is printed in the report.
    */
   void testFailed (int);

   /**
    * Returns the writing engine used in this class
    */
   Writer * getWriter() { return m_writer; }
private:

/**
 * Writing engine used to print the report.
 */
Writer * m_writer;

   /**
    * Report header used to print the report.
    */
   ReportHeader * m_header;

   /**
    * Report footer used to print the report.
    */
   ReportFooter * m_footer;

   /**
    * Indicates the index of the first column in the results table to
    * insert results for the current test.
    */
   int m_base_col;

   /**
    * Indicates the index of the first column in the results table to
    * insert results for the next test.
    */
   int m_next_base;

   /**
    * The internal table that contains the results of lattice test(s).
    */
   Table m_results;

   /**
    * The configuration of the lattice test(s).
    */
   LatConfig * m_config;

   /**
    * Formatter used to format the results in the table when writing the
    * report.
    */
   DoubleFormatter m_dFormat;

   /**
    * Formatter used to format the dimensions (first column) in the table
    * when writing the report.
    */
   FormatterImpl<int> m_iFormat;
};

}
#endif
