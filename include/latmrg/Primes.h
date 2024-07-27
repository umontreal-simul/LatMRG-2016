#ifndef PRIMES_H
#define PRIMES_H
#include "latmrg/Chrono.h"
#include "latcommon/Types.h"
#include "latmrg/IntFactorization.h"
#include <fstream>


namespace LatMRG {

/**
 * This class provides methods to search for integers \f$m\f$ that are prime,
 * for which the integer \f$r = (m^k-1)/(m-1)\f$ is also prime for a given
 * \f$k\f$, and possibly for which \f$(m-1)/2\f$ is also prime.
 * \anchor REF__Primes_clas_Primes
 *
 */
class Primes {
public:

   /**
    * Constructor.
    */
   Primes();

   /**
    * Destructor.
    */
   ~Primes();

   /**
    * Finds the \f$s\f$ prime integers \f$m<2^e\f$ that are closest to
    * \f$2^e\f$. If `facto` is `true`, then \f$m-1\f$ is factorized in its
    * prime factors. The results are printed on stream `fout`.
    */
   void find (int e, int s, bool facto, std::ofstream & fout);

   /**
    * Finds the \f$s\f$ integers \f$m<2^e\f$ that are closest to
    * \f$2^e\f$, such that \f$m\f$ and \f$r = (m^k-1)/(m-1)\f$ are prime.
    * If `safe` is `true`, \f$(m-1)/2\f$ is also required to be prime. The
    * results are printed on stream `fout`. If `facto` is `true`, then
    * \f$m-1\f$ is factorized in its prime factors. If \f$k=1\f$, \f$r\f$
    * is considered to be prime.
    */
   void find (int k, int e, int s, bool safe, bool facto, std::ofstream & fout);

   /**
    * Finds all integers \f$m\f$, in \f$2^e + c_1 \le m \le2^e + c_2\f$,
    * such that \f$m\f$ and \f$r = (m^k-1)/(m-1)\f$ are prime. If `safe`
    * is `true`, \f$(m-1)/2\f$ is also required to be prime. The results
    * are printed on stream `fout`. If `facto` is `true`, then \f$m-1\f$
    * is factorized in prime factors. If \f$k=1\f$, \f$r\f$ is considered
    * to be prime.
    */
   void find (int k, int e, long c1, long c2, bool safe, bool facto,
              std::ofstream & fout);
private:

   /**
    * Writes the parameters of the find to the stream `fout`.
    */
   void writeHeader (int k, int e, long c1, long c2, bool safe, bool facto,
                        std::ofstream & fout);

   /**
    * Writes the CPU time of the find to the stream `fout`.
    */
   void writeFooter (std::ofstream & fout);
   void find (int k, int e, int s, const MScal & S1, const MScal & S2, bool safe,
              bool facto, std::ofstream & fout);
   Chrono timer;
   IntFactorization ifac;
   void nextM (MScal & m) {
      m -= 2;
      if (0 == m % 5)
         m -= 2;    }
};

}
#endif
