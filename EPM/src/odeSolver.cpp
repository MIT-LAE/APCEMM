/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*                      Early Plume Microphysics                    */
/*                              (EPM)                               */
/*                                                                  */
/* odeSolver Program File                                           */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 1/10/2018                                 */
/* File                 : odeSolver.cpp                             */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "odeSolver.hpp"

namespace EPM
{
     
    template<class System> odeSolver<System>::odeSolver( System system_, Vector_1D &x_ ):
        system( system_ ),
        vars( x_ ),
        currentTime( 0.0 ),
        timeStep( 0.0 )
    {

        /* Constructor */

    } /* End of odeSolver::odeSolver */

    template<class System> odeSolver<System>::~odeSolver()
    {

        /* Destructor */

    } /* End of odeSolver::~odeSolver */

    template<class System> odeSolver<System>::odeSolver( const odeSolver &solver )
    {

        /* Copy */

        system = solver.system;
        vars = solver.vars;
        currentTime = solver.currentTime;
        timeStep = solver.timeStep;

    } /* End of odeSolver::odeSolver */
    
    template<class System> odeSolver<System>& odeSolver<System>::operator=( const odeSolver &solver )
    {

        if ( &solver == this )
            return *this;

        system = solver.system;
        vars = solver.vars;
        currentTime = solver.currentTime;
        timeStep = solver.timeStep;
        return *this;

    } /* End of odeSolver::operator= */

    template<class System> void odeSolver<System>::updateTime( RealDouble t_ )
    {

        currentTime = t_;

    } /* End of odeSolver::updateTime */
    
    template<class System> void odeSolver<System>::updateStep( RealDouble dt_ )
    {

        timeStep = dt_;

    } /* End of odeSolver::updateStep */

    template<class System> template<class Observer> UInt odeSolver<System>::integrate( RealDouble start_time, RealDouble end_time, RealDouble dt, Observer observer )
    {

        unsigned int nStep;

        if ( end_time > start_time ) {
            nStep = boost::numeric::odeint::integrate( system, vars, start_time, end_time, dt, observer );
        } else {
            std::cout << "\nIn odeSolver::integrate: end_time is smaller than start_time!";
            return 0;
        }

        if ( nStep >= 1 ) {
            updateTime( end_time );
        }

        return nStep;

    } /* End of odeSolver::integrate */
    
    template<class System> UInt odeSolver<System>::integrate( RealDouble start_time, RealDouble end_time, RealDouble dt )
    {

        unsigned int nStep;

        if ( end_time > start_time ) {
            nStep = boost::numeric::odeint::integrate( system, vars, start_time, end_time, dt );
        } else {
            std::cout << "\nIn odeSolver::integrate: end_time is smaller than start_time!";
            return 0;
        }

        if ( nStep >= 1 ) {
            updateTime( end_time );
        }

        return nStep;

    } /* End of odeSolver::integrate */

    template<class System> 
    bool odeSolver<System>::done( const RealDouble &x, RealDouble threshold )
    {

        return ( x <= threshold );

    } /* End of odeSolver::done */



    streamingObserver::streamingObserver( Vector_2D &states, Vector_1D &times, UInt write_every ):
        m_write_every( write_every ),
        m_count( 0 ),
        m_states( states ),
        m_times( times )

    {

        /* Constructor */

    } /* End of streamingObserver::streamingObserver */

    streamingObserver::~streamingObserver()
    {

        /* Destructor */

    } /* End of streamingObserver::~streamingObserver */

    streamingObserver::streamingObserver( const streamingObserver &obs )
    {

        /* Copy */

        m_write_every = obs.m_write_every;
        m_count = obs.m_count;
        m_states = obs.m_states;
        m_times = obs.m_times;

    } /* End of streamingObserver::streamingObserver */

    streamingObserver& streamingObserver::operator=( const streamingObserver &obs )
    {

        if ( &obs == this )
            return *this;

        m_write_every = obs.m_write_every;
        m_count = obs.m_count;
        m_states = obs.m_states;
        m_times = obs.m_times;
        return *this;

    } /* End of streamingObserver::operator= */

    void streamingObserver::operator()( const Vector_1D &x, double t )
    {

        if ( ( m_count % m_write_every ) == 0 )
        {
           m_states.push_back( x );
           m_times.push_back( t );
        }

    } /* End of streamingObserver::operator() */

    template class odeSolver<void ( const Vector_1D&, Vector_1D&, RealDouble )>;

}

/* End of odeSolver.cpp */
