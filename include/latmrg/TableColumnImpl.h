/**
 * This class implements the `TableColumn` interface. A `TableColumnImpl` can
 * hold any type of data, given that the `Formatter` passed to its
 * constructor correctly formats the same type of data.
 *
 */

#ifndef TABLECOLUMNIMPL_H
#define TABLECOLUMNIMPL_H
#include "TableColumn.h"
#include "Formatter.h"

#include <stdexcept>
#include <string>
#include <vector>


namespace LatMRG {

/**
 * Constructor. `format` must be a formatter for type `T` and `header` is the
 * name of the column that will be printed.
 */
template <typename T>
class TableColumnImpl : public TableColumn {
public:

   TableColumnImpl (Formatter *format, std::string header);

   /**
    * Same as above, except that the vector of data will be initialized at
    * size `size`.
    */
   TableColumnImpl (Formatter *format, std::string header, unsigned int size);

   /**
    * Destructor.
    */
   ~TableColumnImpl() {}

   /**
    * Defined in `TableColumn`.
    */
   std::string getFormattedValue (int pos);

   /**
    * Defined in `TableColumn`.
    */
   void * getValue (unsigned int pos);

   /**
    * Defined in `TableColumn`.
    */
   void setValue (void *value, unsigned int pos);

   /**
    * Defined in `TableColumn`.
    */
   void insertValue (void *value, unsigned int pos);

   /**
    * Defined in `TableColumn`.
    */
   void addValue (void *value);

   /**
    * Defined in `TableColumn`.
    */
   void removeValue (unsigned int pos);

   /**
    * Defined in `TableColumn`.
    */
   int size();
protected:

/**
 * Type definition for the vector of values in the column.
 */
typedef std::vector<T> vect;

   /**
    * Internal vector which contains the values of the column.
    */
   vect m_data;
};


template <typename T>
TableColumnImpl<T>::TableColumnImpl(Formatter *f, std::string s) : TableColumn(f,s)
{
    m_data.reserve(5);
}

template <typename T>
TableColumnImpl<T>::TableColumnImpl(Formatter *f, std::string s, unsigned int size):
    TableColumn(f,s)
{
    m_data.reserve(size);
}

template <typename T>
std::string TableColumnImpl<T>::getFormattedValue(int i)
{
 // Returns the effective size of the column.
    if(i < 0){
        return m_header;
    }

    void* value = static_cast<void*>(&m_data.at(i));

    return m_format->format(value);
}

template <typename T>
void* TableColumnImpl<T>::getValue(unsigned int i)
{
    T* value;
    void *ret;
    value = &m_data.at(i);

    ret = static_cast<void*>(value);

    return ret;
}

template <typename T>
void TableColumnImpl<T>::insertValue(void* value, unsigned int i)
{
    if(i > m_data.size()){
        throw std::out_of_range("TableColumnImpl");
    }

    //vect::iterator iter;

    T val;

    val = *static_cast<T*>(value);

    m_data.push_back(val);

    //iter = data.begin();
    //iter += i;

    //data.insert(iter, val);
}

template <typename T>
void TableColumnImpl<T>::addValue(void* value)
{
    T val;

    val = *static_cast<T*>(value);

    m_data.push_back(val);
}

template <typename T>
void TableColumnImpl<T>::removeValue(unsigned int i)
{
    if(i >= m_data.size()){
        throw std::out_of_range("TableColumnImpl");
    }

    /*vect::iterator iter;

    iter = m_data.begin();
    iter += i;

    m_data.erase(iter);*/
}

template <typename T>
void TableColumnImpl<T>::setValue(void* value, unsigned int i)
{
    if(i > m_data.size()){
        throw std::out_of_range("TableColumnImpl");
    }

    T val;

    val = *static_cast<T*>(value);

    m_data.at(i) = val;
}

template <typename T>
int TableColumnImpl<T>::size()
{
    return m_data.size();
}
}
#endif
