#ifndef CHRONO_H
#define CHRONO_H
#include <string>

/**
 * \file Chrono.h
 */

/**
 * On a MS-Windows platform, the MS-Windows function `GetProcessTimes` will
 * be used to measure the CPU time used by programs. On Linux\f$|\f$Unix
 * platforms, if the macro <tt>USE_ANSI_CLOCK</tt> is defined, the timers
 * will call the ANSI C `clock` function. When <tt>USE_ANSI_CLOCK</tt> is
 * left undefined, class `Chrono` gets the CPU time used by a program via an
 * alternate non-ANSI C timer based on the POSIX (The Portable Operating
 * System Interface) function `times`, assuming this function is available.
 * The POSIX standard is described in the IEEE Std 1003.1-2001 document (see
 * The Open Group web site at
 * [http://www.opengroup.org/onlinepubs/007904975/toc.htm](http://www.opengroup.org/onlinepubs/007904975/toc.htm)).
 */
#undef USE_ANSI_CLOCK
namespace LatMRG {

/**
 * This class acts as an interface to the system clock to compute the CPU
 * time used by parts of a program. Even though the ANSI/ISO macro
 * <tt>CLOCKS_PER_SEC = 1000000</tt> is the number of clock ticks per second
 * for the value returned by the ANSI-C standard `clock` function (so this
 * function returns the number of microseconds), on some systems where the
 * 32-bit type `long` is used to measure time, the value returned by `clock`
 * wraps around to negative values after about 36 minutes. On some other
 * systems where time is measured using the 32-bit type `unsigned long`, the
 * clock may wrap around to 0 after about 72 minutes. When the macro
 * <tt>USE_ANSI_CLOCK</tt> is undefined, a non-ANSI-C clock is used. On
 * Linux-Unix systems, it calls the POSIX function `times` to get the CPU
 * time used by a program.   On a Windows platform (when the macro
 * <tt>HAVE_WINDOWS_H</tt> is defined), the Windows function
 * `GetProcessTimes` will be used to measure the CPU time used by programs.
 *
 * Every object `Chrono` acts as an independent *stopwatch*. Several such
 * stopwatchs can run at any given time. An object of type `Chrono` must be
 * declared for each of them. The method `init` resets the stopwatch to zero,
 * `val` returns its current reading, and `write` writes this reading to the
 * current output. The returned value includes part of the execution time of
 * the functions from class `Chrono`. The `TimeFormat` allows one to choose
 * the kind of time units that are used.
 *
 * Below is an example of how the functions may be used. A stopwatch named
 * `timer` is declared and created. After 2.1 seconds of CPU time have been
 * consumed, the stopwatch is read and reset. Then, after an additional 330
 * seconds (or 5.5 minutes) of CPU time the stopwatch is read again, printed
 * to the output and deleted.
 *
 * <tt>Chrono timer; <br> \f$\vdots\f$ (<em>suppose 2.1 CPU seconds are used
 * here</em>.)<br> double t = timer.val (Chrono::SEC); \qquad// Here, t =
 * 2.1 <br>timer.init(); <br> \f$\vdots\f$ (<em>suppose 330 CPU seconds are
 * used here</em>.) <br> t = timer.val (Chrono::MIN); \qquad// Here, t = 5.5
 * <br>timer.write (Chrono::HMS); \qquad// Prints: 00:05:30.00 <br></tt>
 *
 */
class Chrono {
public:

   /**
    * Types of units in which the time on a `Chrono` can be read or printed: in
    * seconds (<tt>SEC</tt>), minutes (<tt>MIN</tt>), hours (<tt>HOUR</tt>),
    * days (<tt>DAYS</tt>), or in the `HH:MM:SS.xx` format, with hours, minutes,
    * seconds and hundreths of a second (<tt>HMS</tt>).
    */
   enum TimeFormat { SEC, MIN, HOURS, DAYS, HMS };

   /**
    * Constructor for a stopwatch; initializes it to zero. One may
    * reinitialize it later by calling `init`.
    */
   Chrono();

   /**
    * Destructor.
    */
   ~Chrono() {}

   /**
    * Initializes this stopwatch to zero.
    */
   void init ();

   /**
    * Returns the CPU time measured by this `Chrono` since the last call
    * to `init()`. The parameter `unit` specifies the time unit.
    * Restriction: `unit = HMS` is not allowed here; it will cause an
    * error.
    */
   double val (TimeFormat unit);

   /**
    * Prints, on standard output, the CPU time measured by this `Chrono`
    * since its last call to `init()`. The parameter `unit` specifies the
    * time unit.
    */
   void write (TimeFormat unit);

   /**
    * Returns as a string the CPU time measured by this `Chrono` since its
    * last call to `init()`. The time format used is `HMS`.
    */
   std::string toString();

   /**
    * Returns `true` if this `Chrono` has reached the time `limit` (in
    * seconds), otherwise returns `false`.
    */
   bool timeOver (double limit);
private:

   unsigned long microsec;         // microseconds
   unsigned long second;           // seconds

   void tick();
};

/**
 * Returns the value of the duration from `timer`.
 */
std::string toString (Chrono& timer);

}
#endif
