/**
 * This template class is used as a base class for any class that needs to
 * have observers attached or detached from it. It was created to avoid
 * recoding the attaching and detaching code of observers more than once.
 *
 */

#ifndef SUBJECT_H
#define SUBJECT_H
#include <list>


namespace LatMRG {

/**
 * Constructor.
 */
template <typename T>
class Subject {
public:

   Subject() {}

   /**
    * Destructor.
    */
   virtual ~Subject();

   /**
    * Attaches the observer `observer` to this object.
    */
   void attach (T observer);

   /**
    * Detaches the observer `observer` to this object.
    */
   void detach (T observer);
protected:

   typedef std::list<T> list_ob;

/**
 * The list that contains the observers.
 */
list_ob m_observers;
};

template <typename T>
void Subject<T>::attach (T observer)
{
   m_observers.push_back (observer);
   m_observers.unique ();
}

template <typename T>
void Subject<T>::detach (T observer)
{
   m_observers.remove (observer);
}

template <typename T>
Subject<T>::~Subject()
{
   m_observers.clear();
}
}
#endif
