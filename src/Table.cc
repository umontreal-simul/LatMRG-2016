/*
 Table.cc for ISO C++
 version 1.00
 modified 15/06/06

 Authors: Hicham Wahbi
 Frédérik Rozon
 Richard Simard
*/

#include <string>
#include <vector>
#include <stdexcept>

#include "latmrg/Table.h"

namespace LatMRG
{

//=========================================================================

Table::Table()
{
   m_columns.reserve(5);
}


Table::Table(int size)
{
   m_columns.reserve(size);
}

Table::~Table()
{
   for (unsigned int i = 0; i < m_columns.size(); i++)
      delete m_columns.at(i);
   m_columns.clear();
}

void Table::add(TableColumn* column)
{
   m_columns.push_back(column);
}

void Table::insert(TableColumn* column, unsigned int col)
{
   if (col >= m_columns.size())
      throw std::out_of_range("Table.insert");
   m_iter = m_columns.begin();
   m_iter += col;
   m_columns.insert(m_iter, column);
}

void Table::remove(unsigned int col)
{
   if (col >= m_columns.size())
      throw std::out_of_range("Table.remove");
   m_iter = m_columns.begin();
   m_iter += col;
   m_columns.erase(m_iter);
}

void Table::sort(unsigned int col)
{
   col = 0;
   //To be done
}

int Table::getHeight() const
{
   int h = 0;
   TableColumn* col;

   for (unsigned int i = 0; i < m_columns.size(); i++) {
      col = m_columns[i];
      h = col->size() > h ? col->size() : h;
   }

   return h;
}

TableColumn * Table::operator[](int i)
{
   if (i < 0 || (unsigned) i >= m_columns.size())
      throw std::out_of_range("Table[]");
   return m_columns[i];
}

}
