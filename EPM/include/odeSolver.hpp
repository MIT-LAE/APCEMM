/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*                      Early Plume Microphysics                    */
/*                              (EPM)                               */
/*                                                                  */
/* odeSolver Header File                                            */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 1/10/2018                                 */
/* File                 : odeSolver.hpp                             */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef ODESOLVER_H_INCLUDED
#define ODESOLVER_H_INCLUDED

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <boost/numeric/odeint.hpp>

#include "../../Headers/ForwardsDecl.hpp"

namespace EPM
{

    template<class System> class odeSolver;
    class streamingObserver;

}

template<class System> class EPM::odeSolver
{

    public:
    
        odeSolver( System system, Vector_1D &x );
        ~odeSolver();
        odeSolver( const odeSolver &solver );
        odeSolver& operator=( const odeSolver &solver );
        void updateTime( RealDouble t );
        void updateStep( RealDouble dt );

        template<class Observer> UInt integrate( RealDouble start_time, RealDouble end_time, RealDouble dt, Observer observer );
        UInt integrate( RealDouble start_time, RealDouble end_time, RealDouble dt );
        bool done( const RealDouble &x, RealDouble threshold = 0.0 );

    protected:

        const System *system;
        Vector_1D vars;
        RealDouble currentTime;
        RealDouble timeStep;

    private:


};

class EPM::streamingObserver
{

    public:

        streamingObserver( Vector_2D &states, Vector_1D &times, UInt write_every = 100 );
        ~streamingObserver();
        streamingObserver( const streamingObserver &obs );
        streamingObserver& operator=( const streamingObserver &obs );
        void operator()( const Vector_1D &x, double t );

        UInt m_write_every;
        UInt m_count;

    protected:

        Vector_2D m_states;
        Vector_1D m_times;

    private:

};

#endif /* ODESOLVER_H_INCLUDED */
