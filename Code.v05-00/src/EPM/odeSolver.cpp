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

#include "EPM/odeSolver.hpp"

namespace EPM
{
     
    template<class System> odeSolver<System>::odeSolver( System system_, Vector_1D &x_, bool adapt, bool stop_, double t ):
        system( system_ ),
        vars( x_ ),
        currentTime( 0.0 ),
        timeStep( 0.0 ),
        adaptive( adapt ),
        controlled_stepper(),
        stop( stop_ ),
        threshold( t )
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

    template<class System> void odeSolver<System>::updateTime( double t_ )
    {

        currentTime = t_;

    } /* End of odeSolver::updateTime */
    
    template<class System> void odeSolver<System>::updateStep( double dt_ )
    {

        timeStep = dt_;

    } /* End of odeSolver::updateStep */

    template<class System> UInt odeSolver<System>::integrate( double start_time, double end_time, double dt, streamingObserver observer )
    {

        unsigned int nStep;

        if ( end_time > start_time ) {
            if ( adaptive ) {
                nStep = boost::numeric::odeint::integrate_adaptive( boost::numeric::odeint::make_controlled< error_stepper_type >( EPM_ATOLS, EPM_RTOLS ), system, vars, start_time, end_time, dt, observer );
                //boost::find_if( boost::numeric::odeint::make_adaptive_range( boost::numeric::odeint::make_controlled< error_stepper_type>( EPM_ATOLS, EPM_RTOLS ), system, vars, start_time, end_time, dt, observer), done );
            } else {
                nStep = boost::numeric::odeint::integrate( system, vars, start_time, end_time, dt, observer );
            }
        } else {
            std::cout << "\nIn odeSolver::integrate: end_time is smaller than start_time!";
            nStep = 0;
        }

        if ( nStep >= 1 ) {
            updateTime( end_time );
        }

        return nStep;

    } /* End of odeSolver::integrate */
    
    template<class System> UInt odeSolver<System>::integrate( double start_time, double end_time, double dt )
    {

        unsigned int nStep;

        if ( end_time > start_time ) {
            if ( adaptive ) {
                nStep = boost::numeric::odeint::integrate_adaptive( controlled_stepper, system, vars, start_time, end_time, dt );
            } else {
                nStep = boost::numeric::odeint::integrate( system, vars, start_time, end_time, dt );
            }
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
    void odeSolver<System>::updateThreshold( double t )
    {

        threshold = t;

    } /* End of odeSolver::updateThreshold */

        
    template<class System> 
    bool odeSolver<System>::done( const double &x )
    {

        return ( x <= threshold );

    } /* End of odeSolver::done */

    template<class System>
    Vector_1D odeSolver<System>::getState() const
    {

        return vars;

    } /* End of odeSolver::getState */

    streamingObserver::streamingObserver( Vector_2D &states, Vector_1D &times, std::vector<UInt> indices, std::string filename, UInt write_every ):
        m_write_every( write_every ),
        fileName( filename ),
        m_states( states ),
        m_times( times ),
        m_indices( indices )

    {

        /* Constructor */

    } /* End of streamingObserver::streamingObserver */

    streamingObserver::~streamingObserver( )
    {

        /* Destructor */

    } /* End of streamingObserver::~streamingObserver */

    streamingObserver& streamingObserver::operator=( const streamingObserver &obs )
    {

        if ( &obs == this )
            return *this;

        m_write_every = obs.m_write_every;
        m_count = obs.m_count;
        m_states = obs.m_states;
        m_times = obs.m_times;
        fileName = obs.fileName;
        return *this;

    } /* End of streamingObserver::operator= */

    void streamingObserver::operator()( const Vector_1D &x, double t )
    {

        if ( ( m_count % m_write_every ) == 0 ) {
        
            m_states.push_back( x );
            m_times.push_back( t );

            m_count++;

        }

    } /* End of streamingObserver::operator() */

    double streamingObserver::getLastElement( ) const
    {

        return m_states[m_count][0];

    } /* End of streamingObserver::getLastElement */

    void streamingObserver::print2File( ) const
    {

        std::ofstream file (fileName);

        if ( file.is_open() == 0 ) {

            std::cout << "\nIn streamingObserver::streamingObserver: Couldn't open " << fileName << "!\n";

        }
        else {

            unsigned int counter;
            const char* sep = ", ";
            const unsigned int prec = 6;

            double n_air; 

            /* Variable list: 
             * - Temperature [K]
             * - Water molecular concentration [molecules/cm^3] 
             * - Saturation with respect to ice [-]
             * - Saturation with respect to liquid water [-]
             * - TBC ...
             * - */

            file << std::setw(prec+8) << "Time [s], ";
            file << std::setw(prec+8) << "Tracer [-], ";
            file << std::setw(prec+8) << "Temp. [K], ";
            file << std::setw(prec+8) << "Pres. [Pa], ";
            file << std::setw(prec+8) << "H2O [/cm3], ";
            file << std::setw(prec+8) << "RH_i [-], ";
            file << std::setw(prec+8) << "RH_w [-], ";
            file << std::setw(prec+8) << "SO4 [/cm3], ";
            file << std::setw(prec+8) << "SO4g [/cm3], ";
            file << std::setw(prec+8) << "SO4l [/cm3], ";
            file << std::setw(prec+8) << "SO4s [/cm3], ";
            file << std::setw(prec+8) << "SO4Sat [-], ";
            file << std::setw(prec+8) << "HNO3 [/cm3], ";
            file << std::setw(prec+8) << "HNO3Sat[-], ";
            file << std::setw(prec+8) << "Part[/cm3], ";
            file << std::setw(prec+8) << "Rad[mum], ";
            file << std::setw(prec+8) << "Theta1[-], ";
            file << std::setw(prec+8) << "Theta2[-], ";

            /* New line */
            file << "\n";
            file << std::setfill('-') << std::setw(18*(prec+8)) << "-";

            file << std::setfill(' ');

            for ( counter = 0; counter < m_states.size(); counter++ ) {
                /* New line */
                file << "\n";

                file << std::scientific << std::setprecision(prec) << std::setfill(' ');

                /* Output data */

                /* Print time [s] */
                file << m_times[counter];
                file << sep;

                /* Print tracer dilution ratio [-] */
                file << m_states[counter][m_indices[0]];
                file << sep;

                /* Print temperature [K] */
                file << m_states[counter][m_indices[1]];
                file << sep;
                
                /* Print pressure [Pa] */
                file << m_states[counter][m_indices[2]];
                file << sep;
                
                /* Compute number concentration of air for conversions sake */
                n_air = m_states[counter][m_indices[2]] / (physConst::kB * m_states[counter][m_indices[1]] * 1.0e6) ; 

                /* Print gaseous water molecular concentration [molec/cm^3] */
                file << m_states[counter][m_indices[3]] * n_air ;
                file << sep;

                /* Print rel. humidities [-] */
                file << m_states[counter][m_indices[3]] * m_states[counter][m_indices[2]] / physFunc::pSat_H2Os( m_states[counter][m_indices[1]] );
                file << sep;
                file << m_states[counter][m_indices[3]] * m_states[counter][m_indices[2]] / physFunc::pSat_H2Ol( m_states[counter][m_indices[1]] );
                file << sep;
                
                /* Print gaseous SO4 molecular concentration [molec/cm^3] */
                file << m_states[counter][m_indices[4]] * n_air ;
                file << sep;
                
                /* Print gaseous SO4 gaseous molecular concentration [molec/cm^3] */
                file << m_states[counter][m_indices[6]] * n_air;
                file << sep;

                /* Print gaseous SO4 liquid molecular concentration [molec/cm^3] */
                file << m_states[counter][m_indices[5]] * n_air;
                file << sep;
                
                /* Print gaseous SO4 on part [molec/cm^3] */
                file << m_states[counter][m_indices[7]] * m_states[counter][m_indices[2]] / ( physConst::kB * m_states[counter][m_indices[1]] * 1.0E+06 ) ;
                file << sep;

                /* Print SO4 saturation [-] */
                file << ( m_states[counter][m_indices[5]] + m_states[counter][m_indices[6]] ) * m_states[counter][m_indices[2]] / physFunc::pSat_H2SO4( m_states[counter][m_indices[1]] );
                file << sep;
                
                /* Print gaseous HNO3 molecular concentration [molec/cm^3] */
                file << m_states[counter][m_indices[8]] * m_states[counter][m_indices[2]] / ( physConst::kB * m_states[counter][m_indices[1]] * 1.0E+06 ) ;
                file << sep;
                
                /* Print HNO3 saturation [-] */
                file << m_states[counter][m_indices[8]] * physConst::kB * m_states[counter][m_indices[1]] * 1.0E+06 / physFunc::pSat_HNO3( m_states[counter][m_indices[1]], \
                        m_states[counter][m_indices[2]] * physConst::kB * m_states[counter][m_indices[1]] * 1.0E+06 );
                file << sep;
                
                /* Print particle concentration [#/cm^3] */
                file << m_states[counter][m_indices[9]] * n_air;
                file << sep;
                
                /* Print particle radius [mum] */
                file << m_states[counter][m_indices[10]] * 1.0E+06;
                file << sep;
                
                /* Print soot coverage [-] */
                file << m_states[counter][m_indices[11]];
                file << sep;
                
                /* Print soot coverage [-] */
                file << m_states[counter][m_indices[12]];
                file << sep;

            }

            file << "\n";
            file.close();

        }

    } /* End of streamingObserver::print2File */

    bool streamingObserver::checkwatersat( ) const
    {

        /* Variable definitions */
        std::size_t counter = 0;
        float RHw;
        bool watersup { false };

        /* Loop over until RHw >= 1.0 */
        while ( !watersup && counter < m_states.size() ) {

            /* Check value of RHw and update watersub */
            RHw = m_states[counter][m_indices[3]] * m_states[counter][m_indices[2]] / physFunc::pSat_H2Ol( m_states[counter][m_indices[1]] );
            if ( RHw >= 1.0 ) {
                watersup = true;
            }

            /* Iterate the counter */
            counter += 1;

        }

        return watersup;

    } /* End of streamingObserver::checkwatersat */

    template class odeSolver<void ( const Vector_1D&, Vector_1D&, double )>;
}

/* End of odeSolver.cpp */
