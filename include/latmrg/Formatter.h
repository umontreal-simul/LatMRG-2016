#ifndef FORMATTER_H
#define FORMATTER_H
#include <string>


namespace LatMRG {

/**
 * This class is an interface that must implemented to format values in a
 * `TableColumn` which composes a `Table`.
 *
 */
class Formatter {
public:

/**
 * Destructor. Does nothing for now.
 */
virtual ~Formatter() {}

   /**
    * Method that must implemented to format `value` into a string.
    */
   virtual std::string format (void *value) = 0;
};

}
#endif
