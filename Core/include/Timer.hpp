/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Timer Header File                                                */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : Timer.hpp                                 */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef TIMER_H_INCLUDED
#define TIMER_H_INCLUDED

#include <ctime>

class Timer
{

    public:

        Timer(bool start_now = false);
        ~Timer();

        void Start(bool reset = false);
        void Stop();

        unsigned long Elapsed( ) const;

    private:

        std::clock_t start, stop;
        bool running;

};

#endif /* TIMER_H_INCLUDED */

