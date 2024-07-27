#ifndef SEARCHER_H
#define SEARCHER_H
#include "FigureOfMerit.h"
#include "latcommon/Util.h"
#include <vector>


namespace LatMRG {

/**
 * This class implements searches to find rank-1 lattices that are good with
 * respect to some measure or figure of merit (see \cite vSLO94a&thinsp;).
 * Given a positive integer \f$n\f$ and a \f$s\f$-dimensional integer vector
 * \f$(a_1, a_2,…,a_s)\f$, where \f$1 \le a_j < n\f$ for each \f$j\f$, the
 * points are defined by
 * \f[
 * \mathbf{u}_i = (i/n)(a_1, a_2, …, a_s) \bmod1
 * \f]
 * for \f$i=0,…,n-1\f$. Here we always choose \f$a_1=1\f$.
 *
 * The merit object in the `Searcher` constructor *must* set the number of
 * points \f$n\f$, the maximal dimension \f$s\f$ of the lattice and some
 * weight factors \f$\gamma_j\f$. Then the search program will examine
 * different lattices with the given \f$n\f$, \f$s\f$ and \f$\gamma_j\f$ in
 * order to find the best amongst those examined. Searches may be exhaustive
 * or random.
 * \remark **Richard:** P. voudrait aussi pouvoir faire des recherches sur
 * les projections et les dimensions successives.
 *
 * One may consider only a subset of the possible lattices according to some
 * criterion. Some of these restricted searches are implemented in classes
 * derived from `Searcher`.
 *
 * <div class="LatSoft-bigskip"></div>
 */
class Searcher {
public:

/**
 * The number of points \f$n\f$, the dimension \f$s\f$, and possibly the
 * \f$s\f$ weight factors \f$\gamma_j\f$ must be given in `merit`. The
 * \f$n\f$ points of the lattice will be generated in the search.
 */
Searcher (const FigureOfMerit* merit);

   /**
    * Destructor.
    */
   virtual ~Searcher();

   /**
    * Exhaustive search to find the lattice with the best (the smallest)
    * figure of merit in dimension \f$s\f$. The search runs over all
    * values of the generator \f$a_j = 1, 2, …, (n-1)\f$ and over all
    * dimensions up to \f$s\f$ inclusively. The first component \f$a_1\f$
    * of the generator is always set to 1; \f$a_0\f$ is unused. The method
    * returns the best value of the figure of merit found.
    */
   double exhaust (int s, long n);

   /**
    * Similar to `exhaust(s)`, except that only values of \f$a_j\f$
    * *relatively prime* to \f$n\f$ are considered.
    */
   double exhaustPrime (int s, long n);

   /**
    * Random search to find the lattice with the best (the smallest)
    * figure of merit in dimension \f$s\f$. At most \f$k\f$ random
    * \f$s\f$-dimensional generators \f$\mathbf{a}\f$ are examined. Each
    * random component \f$a_j\f$ takes values over the integers \f$1, 2,
    * …, (n-1)\f$. The first component \f$a_1\f$ is always set to 1;
    * \f$a_0\f$ is unused. The method returns the best value of the figure
    * of merit found.
    */
   double random (int s, long n, int k);

   /**
    * Similar to `random(s, k)`, except that only values of \f$a_j\f$
    * *relatively prime* to \f$n\f$ are considered.
    */
   double randomPrime (int s, long n, int k);

   /**
    * Returns the best value of the figure of merit found in the last
    * search.
    */
   double getBestVal();

   /**
    * Returns the generator of the lattice which gave the best value of
    * the figure of merit in the last search. The components of this
    * generator are returned as `A[j]`\f$ =a_j\f$, for \f$j = 1, 2, …,
    * s\f$. `A[0]` is unused.
    */
   const std::vector<long>& getBestAs();

   /**
    * Initializes the random number generator used in random searches with
    * the starting seed `seed`. If this method is not called, a default
    * seed is used.
    */
   void initRand (unsigned long seed);
protected:

   std::vector<long> m_bestAs;         // best generator vector for lattices
   double m_bestVal;                   // best value of merit found
   const FigureOfMerit* m_merit;
   void print (std::vector<long>& y, int s, double val);  // for debugging
   virtual double exhaust (int s, long n, bool relPrime);
   virtual double random (int s, long n, int k, bool relPrime);
   virtual bool isPrime (long n);
   virtual bool isPower2 (long n);
private:

   void incr (std::vector<long>& y, long n, int s);
   void incrPrime (std::vector<long>& y, long n, int s, bool power2F);
   double randomPower2 (int s, long n, int k, bool relPrime);
};

}
#endif