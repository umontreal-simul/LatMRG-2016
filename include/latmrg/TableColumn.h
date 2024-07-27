#ifndef TABLECOLUMN_H
#define TABLECOLUMN_H
#include "Formatter.h"
#include <string>


namespace LatMRG {

/**
 * This abstract class is the representation of a table column for class
 * `Table`. A column is made of a string header and a variable number of
 * values. Values in one column must be of the same type. A `Formatter` must
 * be provided to format the values in the column into a string.
 *
 */
class TableColumn {
public:

/**
 * Constructor. `f` is used to transform the column’s values into
 * <tt>string</tt>’s. `h` is the header of this column.
 */
TableColumn (Formatter *f, std::string h) { m_format = f; m_header = h; }

   /**
    * Destructor.
    */
   virtual ~TableColumn() {}

   /**
    * Returns the formatted value of cell `pos` using the formatter passed
    * to the constructor. If `pos` is negative, the column’s header is
    * returned.
    */
   virtual std::string getFormattedValue (int pos) = 0;

   /**
    * Returns the void pointer of the value of cell at position `pos`.
    */
   virtual void* getValue (unsigned int pos) = 0;

   /**
    * Sets the cell at position `pos` to `value`.
    */
   virtual void setValue (void* value, unsigned int pos) = 0;

   /**
    * Inserts `value` at position `pos` in the column.
    */
   virtual void insertValue (void* value, unsigned int pos) = 0;

   /**
    * Appends `value` at the end of the column.
    */
   virtual void addValue (void* value) = 0;

   /**
    * Removes the cell at position `pos`.
    */
   virtual void removeValue (unsigned int pos) = 0;

   /**
    * Returns the effective size of the column.
    */
   virtual int size() = 0;

   /**
    * Returns this column’s header.
    */
   virtual std::string getHeader() const { return m_header; }

   /**
    * Sets this column’s header to `s`.
    */
   virtual void setHeader (std::string s) { m_header = s; }
protected:

/**
 * This column’s header.
 */
std::string m_header;

   /**
    * The formatter used to format the values of this column.
    */
   Formatter* m_format;
};

}
#endif
