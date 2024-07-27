#ifndef MODULUS_H
#define MODULUS_H
#include "latcommon/Types.h"


namespace LatMRG {

/**
 * This class keeps parameters closely associated with a modulus of
 * congruence. Using it, it will not be necessary to recalculate the square
 * roots of large integers, which are used repeatedly in searches for good
 * generators.
 *
 */
class Modulus {
public:

   Modulus ();

/**
 * Constructor with modulus of congruence \f$m\f$.
 */
Modulus (const MScal & m);

   /**
    * Constructor with value \f$m =b^e + c\f$. Restrictions: \f$b>1\f$ and
    * \f$e > 0\f$.
    */
   Modulus (long b, long e, long c);

   /**
    * Destructor.
    */
   virtual ~Modulus ();

   /**
    * Initializes with value \f$m\f$. Computes `mRac` and `mRacNeg`.
    */
   void init (const MScal & m);

   /**
    * Initializes with value \f$m =b^e + c\f$. Restrictions: \f$b>1\f$ and
    * \f$e > 0\f$. Computes `mRac` and `mRacNeg`.
    */
   void init (long b, long e, long c);

   /**
    * Reduces the modulus \f$m\f$ and sets the variable `mRed` to the
    * reduced modulus. The modulus must have the form \f$m=p^e\f$. The
    * multiplier of the LCG is \f$a\f$.
    */
   void reduceM (const MScal & a);

   /**
    * Assumes that \f$m\f$ is a power of a prime \f$p=b\f$, the order \f$k
    * = 1\f$, and the recurrence is homogeneous. Returns `true` iff the
    * maximal period conditions are satisfied.
    */
   bool perMaxPowPrime (const MScal & a);

   /**
    * Value \f$m\f$ of the modulus.
    */
   MScal m;

   /**
    * Reduced value of the modulus. Computed by `reduceM`.
    */
   MScal mRed;

   /**
    * This flag is `true` when \f$m\f$ is prime, otherwise `false`.
    */
   bool primeF;

   /**
    * When this flag is `true`, the value of \f$m\f$ is built out of the
    * three numbers \f$b\f$, \f$e\f$ and \f$c\f$ as described below;
    * otherwise, the flag is set `false`.
    */
   bool threeF;
long b;
   long e;

/**
 * When `threeF` is `true`, then \f$m\f$ is given in the form \f$m = b^e +
 * c\f$; otherwise, \f$b\f$, \f$e\f$ and \f$c\f$ are undefined.
 */
long c;

   /**
    * \f$\sqrt{\lfloor m \rfloor}\f$.
    */
   MScal mRac;

   /**
    * \f$-\sqrt{\lfloor m \rfloor}\f$.
    */
   MScal mRacNeg;
private:

/**
 * The constant \f$b - 1\f$.
 */
MScal bm1;

   /**
    * The constant \f$b^2\f$.
    */
   MScal b2;

   /**
    * Work variables.
    */
   MScal Y, Eight, Four;
};

}
#endif
