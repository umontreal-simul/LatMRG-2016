#ifndef	MRGCOMPONENT_H
#define	MRGCOMPONENT_H
#include "latcommon/Types.h"
#include "latcommon/Const.h"
#include "IntFactorization.h"
#include "Modulus.h"
#include <string>


namespace LatMRG {

/**
 * This class is used to implement a MRG component in a combined MRG. It
 * exists in order to avoid creating numerous relatively heavy `MRGLattice`
 * objects to represent MRG components. Each MRG component is defined by a
 * modulus \f$m\f$, an order \f$k\f$ and a vector of multipliers \f$a\f$,
 * where \f$a[i]\f$ represents \f$a_i\f$. This MRG satisfies the recurrence
 * \f[
 *   x_n = (a_1 x_{n-1} + \cdots+ a_k x_{n-k}) \mbox{ mod } m
 * \f]
 */
class MRGComponent {
public:

/**
 * Constructor with modulus \f$m\f$, vector \f$a\f$ and order \f$k\f$.
 */
MRGComponent (const MScal & m, const MVect & a, int k);

   /**
    * Constructor with modulus \f$m=b^e + c\f$, vector \f$a\f$ and order
    * \f$k\f$.
    */
   MRGComponent (long b, long e, long c, const MVect & a, int k);

   /**
    * Constructor with modulus \f$m\f$ and order \f$k\f$. Arguments
    * `decom1` and `decor` refer to the prime factor decomposition of
    * \f$m-1\f$ and \f$r=(m^k-1)/(m-1)\f$, respectively. If `decor` equals
    * `DECOMP`, the constructor will factorize \f$r\f$. If `decor` equals
    * <tt>DECOMP_WRITE</tt>, the constructor will factorize \f$r\f$ and
    * write the prime factors to file `filer`. If `decor` equals
    * <tt>DECOMP_READ</tt>, the constructor will read the factors of
    * \f$r\f$ from file `filer`. If `decor` equals <tt>DECOMP_PRIME</tt>,
    * \f$r\f$ is assumed to be prime. Similar considerations apply to
    * `decom1` and `filem1` with respect to \f$m-1\f$.
    */
   MRGComponent (const MScal & m, int k, LatCommon::DecompType decom1,
                 const char *filem1,     LatCommon::DecompType decor,
                 const char *filer);

   /**
    * Constructor similar to the above, except that the modulus of
    * congruence \f$m\f$ is inside the object `modul`.
    */
   MRGComponent (Modulus & modul, int k, LatCommon::DecompType decom1,
                 const char *filem1,     LatCommon::DecompType decor,
                 const char *filer);

   /**
    * Destructor.
    */
   ~MRGComponent();

   /**
    * Copy constructor;
    */
   MRGComponent (const MRGComponent & comp);

   /**
    * Assignment operator.
    */
   MRGComponent & operator= (const MRGComponent & comp);

   /**
    * Sets the multipliers of the recurrence to \f$A\f$.
    */
   void setA (const MVect & A);

   /**
    * Returns `true` if coefficients \f$A\f$ give a MRG with maximal
    * period; returns `false` otherwise.
    */
   bool maxPeriod (const MVect & A);

   /**
    * Returns `true` if coefficients \f$A\f$ give a MRG with maximal
    * period; returns `false` otherwise. This method supposes that
    * condition 1 is `true` and tests only conditions 2 and 3. See method
    * `isPrimitive` of class `PolyPE` on page (FIXME: page#) of this
    * guide.
    */
   bool maxPeriod23 (const MVect & A);

   /**
    * The prime factor decomposition of \f$m-1\f$.
    */
   IntFactorization ifm1;

   /**
    * The prime factor decomposition of \f$r=(m^k-1)/(m-1)\f$, where
    * \f$k\f$ is the order of the recurrence.
    */
   IntFactorization ifr;

   /**
    * The modulus \f$m\f$ of the recurrence.
    */
   Modulus module;

   /**
    * Returns the value of the modulus \f$m\f$ of the recurrence.
    */
   MScal getM() { return module.m; }

   /**
    * The order \f$k\f$ of the recurrence.
    */
   int k;

   /**
    * The multipliers \f$a_i\f$ of the recurrence, \f$i = 1, …, k\f$.
    */
   MVect a;

   /**
    * The length of the period \f$\rho\f$ for this MRG. For now, the
    * recurrence is assumed to correspond to a primitive polynomial and
    * `rho` is calculated as
    * \f[
    * \rho= m^k - 1
    * \f]
    * This value is calculated by `MRGLatticeFactory` and stored here for
    * simplicity.
    */
   MScal rho;

   /**
    * Value needed for the calculation of the multipliers of a combined
    * MRG. It is defined by
    * \f[
    *   n_j = (m/m_j)^{-1} \mbox{ mod } m_j \qquad\mbox{ for } j = 1,…,J,
    * \f]
    * where \f$n_j = \f$ `nj`, \f$m_j\f$ is this object’s modulus `m`,
    * \f$m\f$ is the calculated modulus for the combined MRG (see class
    * <tt>MRGLatticeFactory</tt>), and \f$(m/m_j)^{-1} \mbox{ mod } m_j\f$
    * is the inverse of \f$m/m_j\f$ modulo \f$m_j\f$. This value is
    * calculated by `MRGLatticeFactory` and stored here for simplicity.
    */
   MScal nj;

   /**
    * Contains the starting state of the component for the case when the
    * lattice type is `ORBIT`. It is made of \f$k\f$ numbers.
    */
   MVect orbitSeed;

   /**
    * Returns this object as a string.
    */
   std::string toString ();
private:

/**
 * Does the same as the constructor above with similar arguments.
 */
void init (const MScal & m, int k, LatCommon::DecompType decom1,
              const char *filem1,     LatCommon::DecompType decor,
              const char *filer);
};

}
#endif
