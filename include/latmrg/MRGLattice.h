#ifndef MRGLATTICE_H
#define MRGLATTICE_H
#include "latcommon/Types.h"
#include "latcommon/Const.h"
#include "latcommon/Lacunary.h"
#include "latmrg/IntLattice.h"
#include "latmrg/MRGComponent.h"
#include <string>


namespace LatMRG {

/**
 * This class implements lattice basis built from multiple recursive linear
 * congruential generators (MRGs). One must first call the constructor with a
 * given congruence modulus \f$m\f$, a given order \f$k\f$ for the
 * recurrence, and a maximal dimension for the basis. One must then build the
 * lattice basis associated to a vector of multipliers for a given dimension.
 * Each MRG is defined by a vector of multipliers \f$A\f$, where \f$A[i]\f$
 * represents \f$a_i\f$. This MRG satisfies the recurrence
 * \f[
 *   x_n = (a_1 x_{n-1} + \cdots+ a_k x_{n-k}) \mod m.
 * \f]
 */
class MRGLattice: public LatMRG::IntLattice {
public:

/**
 * Constructor with modulus of congruence \f$m\f$, order of the recurrence
 * \f$k\f$, multipliers \f$a\f$, maximal dimension `MaxDim`, and lattice type
 * `Latt`. Vectors and (square) matrices of the basis have maximal dimension
 * `maxDim`, and the indices of vectors and matrices vary from dimension 1 to
 * `maxDim`. The norm to be used for the basis vectors is `norm`.
 */
MRGLattice (const MScal & m, const MVect & a, int maxDim, int k,
               LatCommon::LatticeType latt,
               LatCommon::NormType norm = LatCommon::L2NORM);

   /**
    * As in the constructor above but the basis is built for the lacunary
    * indices `lac`.
    */
   MRGLattice (const MScal & m, const MVect & a, int maxDim, int k, BVect & lac,
               LatCommon::LatticeType latt,
               LatCommon::NormType norm = LatCommon::L2NORM);

   /**
    * Copy constructor. The maximal dimension of the created basis is set
    * equal to <tt>Lat</tt>’s current dimension.
    */
   MRGLattice (const MRGLattice & Lat);

   /**
    * Assigns `Lat` to this object. The maximal dimension of this basis is
    * set equal to <tt>Lat</tt>’s current dimension.
    */
   MRGLattice & operator= (const MRGLattice & Lat);

   /**
    * Destructor.
    */
   ~MRGLattice();

   /**
    * Cleans and releases memory used by this object.
    */
   void kill();

   /**
    * Builds the basis in dimension \f$d\f$.
    */
   virtual void buildBasis (int d);

   /**
    * Increments the dimension of the basis by 1 by calling either
    * `incDimBasis` or `incDimLaBasis`.
    */
   virtual void incDim();

   /**
    * Returns `true` for the case of lacunary indices, returns `false` for
    * non-lacunary indices.
    */
   bool isLacunary() const { return m_lacunaryFlag; }

   /**
    * Returns the \f$j\f$-th lacunary index.
    */
   BScal & getLac (int j);

   /**
    * Sets the lacunary indices for this lattice to `lat`.
    */
   virtual void setLac (const LatCommon::Lacunary & lat);

   /**
    * \name Sets and gets the values of <tt>m_rho</tt> and <tt>m_lossRho</tt>.
    *
    * @{
    */
   MScal getRho() const { return m_rho; }
   MScal getLossRho() const { return m_lossRho; }
   void setRho (const MScal & val) { m_rho = val; }
   void setLossRho (const MScal & val) { m_lossRho = val; }
   /*
    * @}
    */

   /**
    * Returns a non-mutable copy of the multipliers (coefficients) of the
    * MRG.
    */
   const MVect & getCoef() const { return m_aCoef; }

   /**
    * Returns the vector of multipliers \f$A\f$ as a string.
    */
   std::string toStringCoef() const;

protected:

   /**
    * Initializes a square matrix of order \f$k\f$. This initial matrix contains
    * a system of generators for the given group of states.
    */
   void initStates ();

   /**
    * Initializes some of the local variables.
    */
   void init();

   /**
    * Initializes this object when the lattice type is `ORBIT`.
    */
   void initOrbit();
   void insertion (BMat & Sta);
   void lemme2 (BMat & Sta);

   /**
    * For debugging purposes.
    */
   void trace (char* msg, int d);

   /**
    * Increments the basis by 1 in case of non-lacunary indices.
    */
   virtual void incDimBasis ();

   /**
    * Increments the basis by 1 in case of lacunary indices.
    */
   void incDimLaBasis (int);

   /**
    * Builds the basis of the MRG recurrence in case of non-lacunary
    * indices.
    */
   void buildNaBasis (int d);

   /**
    * Builds the basis of the MRG recurrence in case of lacunary indices.
    */
   void buildLaBasis (int d);

   /**
    * \name Used for the calculation of a combined MRG.
    *
    * @{
    */
   MScal m_lossRho;
   MScal m_rho;
   /**
    * @}
    */

   /**
    * The coefficients of the recurrence.
    */
   MVect m_aCoef;

   /**
    * Indicates which lattice or sublattice is analyzed.
    */
   LatCommon::LatticeType m_latType;

   /**
    * Is `true` in the case of lacunary indices, `false` otherwise.
    */
   bool m_lacunaryFlag;

   /**
    * Contains the lacunary indices when `LacunaryFlag` is `true`,
    * otherwise is undefined.
    */
   LatCommon::Lacunary m_lac;


   /**
    * Work variables.
    *
    * @{
    */
    MScal m_t4, m_t5, m_t6, m_t7, m_t8, m_e;
    MVect m_xi;
    /**
     * @}
     */

   /**
    * \f$\clubsuit\f$ To be completed.
    */
   BMat m_sta, m_wSI;

   /**
    * When the flag <tt>m_ip[i]</tt> is `true`, the \f$i\f$-th diagonal
    * element of matrix <tt>m_sta</tt> is non-zero (modulo \f$m\f$) and
    * divides \f$m\f$. Otherwise (when <tt>m_ip[i]</tt> is
    * <tt>false</tt>), the \f$i\f$-th line of matrix <tt>m_sta</tt> is
    * identically 0.
    */
   bool *m_ip;
};

}
#endif
