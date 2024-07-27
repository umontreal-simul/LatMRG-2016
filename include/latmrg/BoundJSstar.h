#ifndef BOUNDJSSTAR_H
#define BOUNDJSSTAR_H
#include "Discrepancy.h"


namespace LatMRG {

/**
 * This class computes bounds on the weighted star discrepancy for lattice
 * rules constructed with the CBC algorithm as proposed by Joe and Sinescu
 * \cite rJOE06a, \cite rSIN08a, \cite rSIN08b&thinsp;. In the case of
 * product weights, for a lattice rule with \f$n\f$ points, this bound,
 * denoted by \f$D^*_{n,\boldsymbol {\gamma}}\f$, is given by
 * \anchor REF__BoundJSstar_eq_sinescu1
 * \f{align}{
 *    e_{n,d}^2(\mathbf{z}) 
 *    & 
 *   =
 *  \frac{1}{n} \sum_{k=0}^{n-1} \prod_{j=1}^d \left( \beta_j + \gamma_j \sideset{}’\sum_{-N/2<h\le N/2}\; \frac{e^{i2\pi hkz_j/n}}{|h|}\right) - \prod_{j=1}^d \beta_j \tag{eq.sinescu1} 
 *  \\ 
 *    D^*_{n,\boldsymbol {\gamma}}(\mathbf{z}) 
 *    & 
 * \le
 *  \frac{e_{n,d}^2(\mathbf{z})}{2} + \prod_{j=1}^d\beta_j-\prod_{j=1}^d\left(\beta_j-\gamma_j/N\right), \nonumber
 * \f}
 * where \f$\beta_j = 1 + \gamma_j\f$ for all \f$j=1,…,d\f$,
 * \f$i=\sqrt{-1}\f$, \f$N = nr\f$, and \f$\mathbf{z} = \{z_1, z_2, …,
 * z_d\}\f$ is the generating vector of the lattice rule.
 * \remark **Richard:** Il faudrait uniformiser la notation: \f$z_j\f$ ou
 * \f$a_j\f$ comme dans le reste de LatMRG. Joe et Sinescu utilise toujours
 * \f$z_j\f$. Pierre croit que cela devrait être \f$a_j\f$ comme dans tous
 * ses articles.
 *
 *  Each component \f$z_j\f$ of the generating vector is taken from the set
 * \f$\{1, 2, …, n-1\}\f$ such that \f$z_j\f$ is relatively prime with
 * \f$n\f$. By \f$\boldsymbol {\gamma}\f$, we denote the weights. The symbol
 * \f$\sum^{\prime}\f$ means that the \f$h=0\f$ term is omitted from the
 * sum.
 *
 * The constructors below precomputes the sum
 * \f[
 * \sideset{}’\sum_{-N/2<h\le N/2}\; \frac{e^{i2\pi hkz_j/n}}{|h|}
 * \f]
 * for all possible values of \f$z_j\f$ and keeps them in an array. This will
 * accelerate the computation of the bounds tremendously, when we do a search
 * over many values of \f$z_j\f$.
 * \remark **Richard:**  Rajouter une méthode `static compute` pour calculer
 * la borne sans précalculer toutes ces sommes de Fourier. On ne veut pas les
 * précalculer dans les cas on a besoin de la borne pour une seule valeur de
 * \f$z_j\f$. Par exemple, on pourrait vouloir calculer la borne pour un
 * \f$z_j\f$ et plusieurs valeurs différentes de \f$n\f$. (Dixit Pierre).
 * Idem pour les autres bornes.
 *
 * <div class="LatSoft-bigskip"></div>
 */
class BoundJSstar: public Discrepancy {
public:

/**
 * Constructor. The bound {@link REF__BoundJSstar_eq_sinescu1 (eq.sinescu1)}
 * will be computed for rank-1 lattice rules with \f$n\f$ points and weights
 * \f$\gamma_j\f$ = `gamma[j]` for \f$j = 1, 2, …, d\f$. Element `gamma[0]`
 * is unused. `prime` indicates whether \f$n\f$ is a prime number
 * (<tt>true</tt>) or not (<tt>false</tt>). In this case, \f$N=n\f$ in eq.
 * {@link REF__BoundJSstar_eq_sinescu1 (eq.sinescu1)}.
 */
BoundJSstar (long n, bool prime, double gamma[], int d);

   /**
    * Constructor for a lattice rule with \f$n\f$ points as above, except
    * that now \f$nr=N\f$.
    */
   BoundJSstar (long n, long r, bool prime, double gamma[], int d);

   /**
    * Destructor.
    */
   ~BoundJSstar();

   /**
    * Returns the (constructor) parameters of this object as a string.
    */
   const std::string toString() const;

   /**
    * Computes and returns the bound {@link REF__BoundJSstar_eq_sinescu1
    * (eq.sinescu1)} in dimension \f$d\f$ for the lattice rule with
    * generating vector \f$z_j =\f$ `z[j]`, for \f$j = 1, 2, …, d\f$.
    * Element `z[0]` is unused.
    */
   double compute (long z[], int d);
protected:

/**
 * Sets the factors \f$\beta_j\f$ = `1 + gamma[j]`, for \f$j = 1, 2, …,
 * d\f$. Element \f$\beta_0\f$ is unused.
 */
void setBeta (double gamma[], int d);
double *m_Beta;

/**
 * Initialization method used to compute the \f$\beta_j\f$ and the Fourier
 * factors \f$C_k\f$.
 */
// Factors beta_j for j = 1, 2,..., d


private:

   void init (double gamma[], int d);
};

}
#endif // BOUNDJSSTAR_H