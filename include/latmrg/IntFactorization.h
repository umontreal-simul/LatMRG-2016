#ifndef INTFACTORIZATION_H
#define INTFACTORIZATION_H
#include <vector>
#include <list>
#include <string>
#include <stdexcept>
 
#include "latcommon/Types.h"
#include "latcommon/Const.h"
#include "latcommon/IntFactor.h"


namespace LatMRG {

/**
 * Given any natural integer \f$n\f$, there is a unique decomposition in
 * prime factors of the form
 * \f[
 *   n = p_1^{\nu_1}p_2^{\nu_2}\cdots p_s^{\nu_s}
 * \f]
 * where \f$p_i\f$ is a prime integer with \f$\nu_i\f$ its multiplicity, and
 * where the factors are sorted in increasing order. In the case of very
 * large integers, it may not be possible to find all the prime factors
 * within a reasonable amount of time. In that case, a similar decomposition
 * to the above may be used with some of the factors composite.
 *
 * The class `IntFactorization` implements the decomposition of integers in
 * factors, preferably prime (see class <tt>IntFactor</tt>). It contains
 * functions to factorize an integer in prime factors, to sort and print the
 * list of its factors. \anchor REF__IntFactorization_IntFactorization_class
 * Integers are factorized by calling the MIRACL software
 * \cite iSCO03a&thinsp;, which uses many different methods in succession to
 * effect the factorization.
 *
 */
class IntFactorization {  
public:

   /**
    * Integer \f$x\f$ is the number whose "prime" factors decomposition is kept
    * in this object.
    */
   explicit IntFactorization (const MScal & x);

   /**
    * Integer `number` and its prime factors will be read from file
    * `fname` by calling method `read` below. If no argument is given,
    * `number` is initialized to 0.
    */
   explicit IntFactorization (const char *fname = 0);

   /**
    * Destructor.
    */
   ~IntFactorization ();

   /**
    * Copy constructor.
    */
   IntFactorization (const IntFactorization & f);

   /**
    * Assignment operator.
    */
   IntFactorization & operator= (const IntFactorization & f);

   /**
    * Empties the lists of primes factors and set this number to 0.
    */
   void clear ();

   /**
    * Reads the list of (possibly) prime factors of a number from file
    * `f`. The first line contains the number itself. The following lines
    * contain one factor per line: the factor (first field) with its
    * multiplicity (second field) and its status (third field). The status
    * field is written as `P` if the factor is known to be prime, `Q` if
    * the factor is probably prime, `C` if the factor is composite, and
    * `U` if its status is unknown or unimportant. For example, given the
    * number \f$120= 2^3*3*5\f$, then the file must be
    *
    * <center>
    *  <div class="LatSoft-fbox"> <div class="LatSoft-parbox">  120
    * <br>2\qquad3 \qquad P <br>3\qquad1 \qquad P <br>5\qquad1 \qquad P <br>
    * </div> </div>
    * </center>
    */
   void read (const char *f) throw (std::invalid_argument);

   /**
    * Adds factor \f$p\f$ with multiplicity `mult` and prime status `st`
    * to this object.
    */
   void addFactor (const MScal & p, int mult = 1,
                   LatCommon::PrimeType st = LatCommon::UNKNOWN);

   /**
    * Replaces repeated equal factors in `factorList` by one factor with
    * its multiplicity. Also sorts the factors.
    */
   void unique ();

   /**
    * Tries to find all the prime factors of this number.
    */
   void factorize ();

   /**
    * Checks that the number is equal to the product of its factors.
    * Returns `true` if it is, otherwise `false`.
    */
   bool checkProduct () const;

   /**
    * Given the list of prime factors \f$p\f$ of `number`, fills the list
    * of inverse factors with the values <tt>number</tt>/\f$p\f$.
    */
   void calcInvFactors ();

   /**
    * Returns the value of this number.
    */
   MScal getNumber () const { return m_number; }

   /**
    * Returns a non-mutable list of the factors.
    */
   const std::list<LatCommon::IntFactor> & getFactorList () const { return m_factorList; }

   /**
    * Returns a non-mutable list of the inverse factors.
    */
   const std::vector<MScal> & getInvFactorList () const { return m_invFactorList; }

   /**
    * Sets the value of this number to \f$x\f$.
    */
   void setNumber (const MScal & x) { m_number = x; }

   /**
    * Returns the status of this number.
    */
   LatCommon::PrimeType getStatus () const { return m_status; }

   /**
    * Sets the status of this number to \f$s\f$.
    */
   void setStatus (LatCommon::PrimeType s) { m_status = s; }

   /**
    * Returns the list of (possibly) prime factors of this object in the
    * same format as described in method `read` above.
    */
   std::string toString () const;
private:

/**
 * The number whose "prime" factor decomposition is kept in this object.
 */
MScal m_number;

   /**
    * The status of this number, i.e. whether it is prime, composite, ...
    */
   LatCommon::PrimeType m_status;

   /**
    * The list of the "prime" factors in the decomposition of `number`.
    */
   std::list<LatCommon::IntFactor> m_factorList;

   /**
    * Given the list of prime factors \f$p\f$ of `number`, `invFactorList`
    * contains all the sorted values <tt>number</tt>/\f$p\f$ (indexing
    * starts at 0). However, one must have called the function
    * `calcInvFactors` beforehand. For example, if `number` = 24, its
    * prime factors decomposition is \f$24 = 2^3\cdot3\f$, and
    * `invfactorList` = \f$[8, 12]\f$.
    */
   std::vector<MScal> m_invFactorList;

   /**
    * Nested class used to sort the prime factors of a number in
    * increasing order. Returns `true` if the factor of `f1` is smaller
    * than the factor of `f2`; otherwise returns `false`.
    */
   class CompFactor {
   public:
      bool operator() (const LatCommon::IntFactor & f1, const LatCommon::IntFactor & f2) {
         return f1.getFactor() < f2.getFactor(); }
   };

   /**
    * Sorts the list of factors in increasing order, the smallest factor
    * first.
    */
   void sort () { CompFactor comp; m_factorList.sort (comp); }
};

}
#endif
