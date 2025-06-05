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

#ifndef EPM_MODELS_ORIGINAL_ODESOLVER_H_INCLUDED
#define EPM_MODELS_ORIGINAL_ODESOLVER_H_INCLUDED

#include <boost/range/algorithm.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp>
#include <boost/numeric/odeint/stepper/controlled_runge_kutta.hpp>
#include <boost/numeric/odeint/iterator/adaptive_iterator.hpp>
#include "Util/ForwardDecl.hpp"

typedef boost::numeric::odeint::runge_kutta_fehlberg78<Vector_1D> error_stepper_type;
typedef boost::numeric::odeint::controlled_runge_kutta<error_stepper_type> controlled_stepper_type;

namespace EPM::Models::OriginalImpl
{
    class streamingObserver
    {
    public:
        streamingObserver( Vector_2D &states, Vector_1D &times, std::vector<UInt> indices, std::string fileName, UInt write_every = 100 );
        ~streamingObserver( );
        streamingObserver& operator=( const streamingObserver &obs );
        void operator()( const Vector_1D &x, double t );
        double getLastElement() const;
        void print2File( ) const;
        bool checkwatersat( ) const;

        UInt m_write_every;
        UInt m_count = 0;
        //const char* fileName;
        std::string fileName;

        Vector_2D &m_states;
        Vector_1D &m_times;

    private:
        const std::vector<UInt> m_indices;
    };
}

#endif
