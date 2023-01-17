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
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <boost/range/algorithm.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp>
#include <boost/numeric/odeint/stepper/controlled_runge_kutta.hpp>
#include <boost/numeric/odeint/iterator/adaptive_iterator.hpp>

#include "Util/ForwardDecl.hpp"
#include "Core/Parameters.hpp"
#include "Util/PhysConstant.hpp"
#include "Util/PhysFunction.hpp"

typedef boost::numeric::odeint::runge_kutta_fehlberg78< Vector_1D > error_stepper_type;
typedef boost::numeric::odeint::controlled_runge_kutta< error_stepper_type > controlled_stepper_type;

namespace EPM
{

    template<class System> class odeSolver;
    class streamingObserver;

}

template<class System> class EPM::odeSolver
{

    public:
    
        odeSolver( System system, Vector_1D &x, bool adaptive = 1, bool stop = 0, RealDouble threshold = 0.0 );
        inline odeSolver( System system, Vector_1D &x, bool stop, RealDouble threshold ) {
            odeSolver( system, x, 0, stop, threshold );
        }
        ~odeSolver();
        odeSolver( const odeSolver &solver );
        odeSolver& operator=( const odeSolver &solver );
        void updateTime( RealDouble t );
        void updateStep( RealDouble dt );

        UInt integrate( RealDouble start_time, RealDouble end_time, RealDouble dt, streamingObserver observer );
       
        UInt integrate( RealDouble start_time, RealDouble end_time, RealDouble dt );
       
        void observer( const Vector_1D &x, const RealDouble t );

        void updateThreshold( RealDouble t );

        bool done( const RealDouble &x );
        Vector_1D getState() const;

    protected:

        const System *system;
        Vector_1D vars;
        RealDouble currentTime;
        RealDouble timeStep;
        bool adaptive;
        controlled_stepper_type controlled_stepper;
        bool stop;
        RealDouble threshold;

    private:


};

class EPM::streamingObserver
{

    public:

        streamingObserver( Vector_2D &states, Vector_1D &times, std::vector<UInt> indices, std::string fileName, UInt write_every = 100 );
        ~streamingObserver( );
        streamingObserver& operator=( const streamingObserver &obs );
        void operator()( const Vector_1D &x, double t );
        RealDouble getLastElement() const;
        void print2File( ) const;
        bool checkwatersat( ) const;

        UInt m_write_every;
        UInt m_count = 0;
        //const char* fileName;
        std::string fileName;

        Vector_2D &m_states;
        Vector_1D &m_times;

    protected:

    private:
        
        const std::vector<UInt> m_indices;

};

#endif /* ODESOLVER_H_INCLUDED */
