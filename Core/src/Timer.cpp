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

#include "Timer.hpp"

Timer::Timer( bool start_now )
    : start( 0 ), stop( 0 ), running( false )
{

    /* Constructor */

    if ( start_now )
        Start( true );

} /* End of Timer::Timer */

Timer::~Timer( )
{

    /* Destructor */
     
} /* End of Timer::~Timer */

void Timer::Start( bool reset )
{

    if ( !running ) {
        
        if ( reset )
            start = clock();

        running = true;

    }


} /* End of Timer::Start */

void Timer::Stop( )
{

    if ( running ) {
     
        stop = clock();
        running = false;

    }

} /* End of Timer::Stop */

unsigned long Timer::Elapsed( ) const
{

    clock_t ticks = ( running ? std::clock() : stop ) - start;
    double seconds = (double)ticks / CLOCKS_PER_SEC;
    unsigned long ms = seconds * 1000;

    return ms;

} /* End of Timer::Elapsed */

/* End of Timer.cpp */
