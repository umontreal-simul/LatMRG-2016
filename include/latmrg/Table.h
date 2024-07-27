#ifndef TABLE_H
#define TABLE_H
#include "TableColumn.h"
#include <vector>


namespace LatMRG {

/**
 * This class implements a table of values. It is made of an arbitrary number
 * of <tt>TableColumn</tt>’s, each one containing one type of data. The table
 * can be printed in several formats using a class derived from `Writer`.
 *
 */
class Table {
public:

/**
 * Constructor.
 */
Table();

   /**
    * Same as above, except that enough memory will be reserved for
    * <tt>size TableColumn</tt>’s.
    */
   Table (int size);

   /**
    * Destructor.
    */
   virtual ~Table();

   /**
    * Appends a `TableColumn` at the end of the table.
    */
   void add (TableColumn * column);

   /**
    * Inserts a `TableColumn` in the table at position `pos`.
    */
   void insert (TableColumn * column, unsigned int pos);

   /**
    * Removes the column in the table at position `pos`.
    */
   void remove (unsigned int pos);

   /**
    * NOT YET IMPLEMENTED.
    */
   void sort (unsigned int pos);

   /**
    * Returns the height of the highest column in the table.
    */
   int getHeight() const;

   /**
    * Return the numbers of columns in the table.
    */
   int size() const { return m_columns.size(); }

   /**
    * Returns the column at position `pos`.
    */
   TableColumn * operator[] (int pos);
protected:

/**
 * The vector containing the columns for the table.
 */
std::vector<TableColumn *> m_columns;

   /**
    * The iterator for the column’s vector.
    */
   std::vector<TableColumn *>::iterator m_iter;
};

}
#endif
