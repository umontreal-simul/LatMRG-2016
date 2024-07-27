#ifndef COEFAPPROXFACT_H
#define COEFAPPROXFACT_H
#include "Coefficient.h"
#include "Zone.h"


namespace LatMRG {

/**
 * This class chooses the coefficients \f$a_i\f$ of a recursive generator of
 * order \f$k\f$ such that the coefficients satisfy the condition \f$|a_i|
 * (m\mod|a_i|) < m\f$, called "approximate factoring", for each \f$i\f$.
 * MRGs are often easier to implement under this condition
 * \cite rLEC90a&thinsp;.
 *
 */
class CoefApproxFact: public Coefficient {
public:

   /**
    * Constructor. \f$m\f$ is the modulus of congruence of the generator.
    */
   CoefApproxFact (Zone *Z, MScal & m);

   /**
    * Destructor.
    */
   ~CoefApproxFact();

   /**
    * Sets \f$A[i]\f$ as a function of \f$q\f$, depending on the zone and
    * \f$m\f$. \f$i\f$ is unchanged.
    */
   void set (const MScal & q, MVect & A, int & i);
private:

   /**
    * Contains the zone as defined in the constructor. It is a pointer to a zone
    * defined elsewhere.
    */
   Zone *m_Z;

   /**
    * The modulus of congruence \f$m\f$ as defined in the constructor.
    */
   MScal m_m;
};

}
#endif
