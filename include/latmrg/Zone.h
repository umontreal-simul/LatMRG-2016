/* Zone.h for ISO C++ */
#ifndef ZONE_H
#define ZONE_H
#include "SeekConfig.h"


namespace LatMRG {

/**
 * This class implements search zones in parameter space for the coefficients
 * of the recurrences defining generators or lattices.
 *
 */
class Zone {
public:

   /**
    * Possible zone number. <tt>ZONE1</tt> is the case where \f$b
    * < -\sqrt{m}\f$, <tt>ZONE3</tt> is the case where \f$b > \sqrt{m}\f$, and
    * <tt>ZONE2</tt> is the case where \f$ -\sqrt{m} \le b \le\sqrt{m}\f$.
    */
   enum ZoneType { ZONE1, ZONE2, ZONE3, NZONES };

   /**
    * Constructor.
    */
   Zone ();

   /**
    * Destructor.
    */
   ~Zone ();

   /**
    * Initializes the research zones for the multiplier \f$a_i\f$ of
    * <tt>comp</tt> which is the \f$s\f$-th component of the combined
    * generator. In the case of a random search, creates also the region
    * for this multiplier.
    */
   void init (const Component & comp, int s, int i);

   /**
    * Computes the lower bound <tt>inf</tt> of the intersection of the
    * zone with the search interval \f$[b,c]\f$. If <tt>z</tt> is equal to
    * <tt>ZONE1</tt> or <tt>ZONE3</tt>, the computed bound is such that
    * \f$\lfloor m/q\rfloor= c\f$ (approximately). If <tt>z = ZONE2</tt>,
    * the computed bound is \f$b\f$. In any case, the lower bound of the
    * zone is always \f$\le b\f$.
    */
   void calcInfBound (const MScal & b, const MScal & c, const Component & comp,
                      int z, MScal & inf);

   /**
    * Computes the upper bound <tt>sup</tt> of the intersection of the
    * zone with the search interval \f$[b,c]\f$. If <tt>z</tt> is equal to
    * <tt>ZONE1</tt> or <tt>ZONE3</tt>, the computed bound is such that
    * \f$\lfloor m/q\rfloor= b\f$ (approximately). If <tt>z = ZONE2</tt>,
    * the computed bound is \f$\min\{c, \sqrt{m}\}\f$.
    */
   void calcSupBound (const MScal & b, const MScal & c, const Component & comp,
                      int z, MScal & sup);

   /**
    * Returns the lower boundary of this region.
    */
   MScal & getInf() { return inf; }

   /**
    * Returns the upper boundary of this region.
    */
   MScal & getSup() { return sup; }
   MScal & getSupMsH() { return supMsH; }

   /**
    * Returns the value of `frac`.
    */
   double getFrac() { return frac; }

   /**
    * Is <tt>true</tt> if and only if <tt>sup - inf \f$\le\f$ H</tt> (or
    * <tt>Hk</tt>).
    */
   bool smallF;

   /**
    * Returns the zone number for this region.
    */
   ZoneType getNo () { return No; }

   /**
    * Selects a random region and initializes its boundaries. The program
    * will search this region exhaustively.
    */
   void chooseBoundaries (const Component & comp, Zone *zone);

   /**
    * Sets the values of the upper boundaries in the three zones based on
    * \f$m_j\f$ for each of the \f$J\f$ components generators. Also sets
    * the seed for the random number generator used in the random search.
    */
   static void initFrontieres (const SeekConfig & config);

   /**
    * Upper boundary of each zone.
    */
   static MMat Frontiere;

   /**
    * Is <tt>true</tt> if at upper boundary of each zone.
    */
   static const bool DivQ[NZONES];

   /**
    * List of zones.
    */
   Zone *nextZone;

   /**
    * Returns this zone as a string.
    */
   std::string toString();

protected:

   /**
    * The zone number where this region is found.
    */
   ZoneType No;

   /**
    * Lower boundary defining the region.
    */
   MScal inf;

   /**
    * Upper boundary defining the region.
    */
   MScal sup;

   /**
    * The fraction of the acceptable values of the multiplier \f$a\f$
    * lying in this zone.
    */
   double frac;

   /**
    * When <tt>small</tt> is <tt>true</tt>, <tt>supMsH</tt> is equal to
    * <tt>sup\f${}+1-{}\f$H</tt> (or <tt>sup\f${}+1-{}\f$Hk</tt>).
    */
   MScal supMsH;

   /**
    * Work variables.
    */
   MScal T1, T2;
};

}
#endif
