#ifndef SEARCHERKOROBOV_H
#define SEARCHERKOROBOV_H
#include "FigureOfMerit.h"
#include "Searcher.h"


namespace LatMRG {

/**
 * This class implements searches to find the best *Korobov* rank-1 lattice
 * rules with respect to a given figure of merit. Given a positive integer
 * \f$n\f$, a dimension \f$s\f$ and a multiplier \f$a\f$, the
 * \f$s\f$-dimensional points are defined by
 * \f[
 * \mathbf{u}_i = (i/n)(1, a, a^2, …, a^{s-1}) \bmod1
 * \f]
 * for \f$i=0,…,n-1\f$ or as an equivalent rank-1 lattice
 * \f[
 * \mathbf{u}_i = (i/n)(z_1, z_2, …, z_s) \bmod1
 * \f]
 * where \f$z_s = a^{s-1} \bmod n\f$. 
 *
 * The figure of merit object in the constructor *must* set the number of
 * points \f$n\f$, the maximal dimension \f$s\f$ of the lattice and possibly,
 * the weight factors \f$\gamma_j\f$. Then the search program will examine
 * different lattices with \f$n\f$, \f$s\f$ and \f$\gamma_j\f$ fixed in
 * order to find the best amongst those examined.
 *
 * <div class="LatSoft-bigskip"></div>
 */
class SearcherKorobov : public Searcher {
public:

/**
 * The number of points \f$n\f$, the dimension \f$s\f$, and possibly the
 * \f$s\f$ weight factors \f$\gamma_j\f$ must be given in `merit`. The
 * \f$n\f$ points of the lattice will be generated in the search.
 */
SearcherKorobov (const FigureOfMerit* merit);

   /**
    * Destructor.
    *
    * `double exhaust (int s, long n);`
    *
    * Exhaustive search to find the best Korobov lattice (i.e. with the
    * smallest discrepancy) in dimension \f$s\f$. The search runs
    * exhaustively over all values \f$ a = 1, 2, 3, …, (n-1)\f$, where
    * \f$a\f$ is the generator of the lattice. The first component
    * \f$z_1\f$ of the generator is always set to 1; \f$z_0\f$ is unused.
    * The method returns the best value found for the discrepancy.
    *
    * `double exhaustPrime (int s, long n);`
    *
    * Similar to  {@link #exhaust(int) exhaust(s)}, except that only
    * values of \f$a\f$ *relatively prime* to \f$n\f$ are considered.
    *
    * `double random (int s, long n, int k);`
    *
    * Random search to find the best Korobov lattice (i.e. with the
    * smallest discrepancy) in dimension \f$s\f$. \f$k\f$ random values
    * \f$a\f$ are examined as generators of the lattice. The \f$a\f$ takes
    * values over the integers \f$ a = 1, 2, 3, …, (n-1)\f$. The method
    * returns the best value found for the discrepancy.
    *
    * `double randomPrime (int s, long n, int k);`
    *
    * Similar to  {@link #random(int,int) random(s, k)}, except that only
    * values of \f$a\f$ *relatively prime* to \f$n\f$ are considered.
    */
   virtual ~SearcherKorobov();

   /**
    * Returns the generator \f$a\f$ of the lattice which gave the best
    * value of the discrepancy in the last search.
    */
   long getBestA();
protected:

   double exhaust (int s, long n, bool relPrime);
   double random (int s, long n, int k, bool relPrime);
   void print (long a);
private:

   long m_bestA;           // best generator a for the Korobov lattice
   void calcAs (long n, int s, long a, std::vector<long>& z);
   double exhaustPrimePower2 (int s, long n);
};

}
#endif