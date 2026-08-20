#include <iostream>
#include "Util/PlumeModelUtils.hpp"

namespace PlumeModelUtils {
    std::vector<double> BuildTime( const double tStart, const double tEnd, \
                                const double DYN_DT )
    {

        unsigned int nT = 0;

        std::vector<double> timeArray;
        double time = tStart;
        double timeStep, nextTimeStep;

        timeStep = 0.0E+00;
        nextTimeStep = 0.0E+00;

        while ( time < tEnd ) {

            /* For debug purposes, uncomment the next line */
            // std::cout << time/3600 << ", " << timeStep << "  ";

            timeArray.push_back( time );
            timeStep = UpdateTime( time, tStart, DYN_DT, nextTimeStep );
            time += std::min( timeStep, std::abs( ( tEnd - time ) ) );
            nT++;
        }

        timeArray.push_back( time );
        return timeArray;

    } /* End of BuildTime */

    double UpdateTime( double time, const double tStart, \
                   const double DYN_DT, double &nextTimeStep )
    {

        const double default_TimeStep = DYN_DT;
        double timeStep;
        const int LVL_TIMESTEP_REFINEMENT = 0;

        if ( nextTimeStep > 0.0E+00 ) {
            timeStep = nextTimeStep;
            nextTimeStep = 0.0E+00;
            return timeStep;
        }

        if ( ( LVL_TIMESTEP_REFINEMENT >= 1 ) && \
        ( ( time - tStart ) < 3600.0 && ( time - tStart ) >= 0.0 ) ) {
            if ( LVL_TIMESTEP_REFINEMENT == 2 ) {
                /* If default_TimeStep = 10 mins, then:
                * 0.0 - .5 - 1.0 - 1.5 - 2.0 - 2.5 - 3.0 - 3.5 - ... - 59.5 - 60.0 */
                timeStep = default_TimeStep / (double) 20.0; 
            } else if ( LVL_TIMESTEP_REFINEMENT == 1 ) {
                /* If default_TimeStep = 10 mins, then:
                * 0.0 - 5.0 - 10.0 - ... - 55.0 - 60.0 */
                timeStep = default_TimeStep / (double) 2.0;
            } else {
                std::cout << "In UpdateTime:: Wrong value for LVL_TIMESTEP_REFINEMENT = ";
                std::cout << LVL_TIMESTEP_REFINEMENT << std::endl;
                exit(-1);
            }
        } else
            timeStep = default_TimeStep;

        return timeStep;

    } /* End of Updatetime */

    void AdvGlobal( const double time, const double T_UPDRAFT, \
                    const double V_UPDRAFT,                    \
                    double &v_x, double &v_y,                  \
                    double &dTrav_x, double &dTrav_y )
    {

        /* AdvGlobal:
        * Computes domain advection parameters
        *
        * INPUTS:
        * (double) time     : current time since simulation started in [s]
        * (double) T_UPDRAFT: updraft timescale in [s]
        * (double) V_UPDRAFT: initial upward velocity in [s]
        *
        * OUTPUTS:
        * (double) v_x: current horizontal velocity at time in [m/s]
        * (double) v_y: current vertical velocity at time in [m/s]
        * (double) dTrav_x: vertical distance traveled since simulation started [m]
        * (double) dTrav_y: vertical distance traveled since simulation started [m]
        */


        v_x = VX;
        dTrav_x = VX * time;

        if ( SYNPROF == 0 ) {

            /* Velocity step */
            if ( time < T_UPDRAFT ) {
                v_y = V_UPDRAFT;
                dTrav_y = v_y * time;
            }
            else
                v_y = 0.0;

        }
        else if ( SYNPROF == 1 ) {

            /* Exponentially decaying velocity */
            v_y = V_UPDRAFT * exp ( - time / T_UPDRAFT );

            dTrav_y = V_UPDRAFT * T_UPDRAFT * ( 1.0 - exp ( - time / T_UPDRAFT ) );

        }
        else {
            std::string const currFile("AdvGlobal.cpp");
            std::cout << "ERROR: In " << currFile << ":SYNPROF set to " << SYNPROF << std::endl;
        }

    } /* End of AdvGlobal */

    void DiffParam( const double time, const double dt, double &d_x, double &d_y, \
                    const double D_X, const double D_Y )
    {

        /* DiffParam:
        * Computes time-averaged diffusion parameters over [time, time + dt]
        *
        * INPUTS:
        * (double) time: current time since simulation started in [s]
        * (double) dt: transport timestep duration in [s]
        *
        * OUTPUTS:
        * (double) d_x: average horizontal diffusion coefficient over [time, time + dt] in [m^2/s]
        * (double) d_y: average vertical diffusion coefficient over [time, time + dt] in [m^2/s]
        */

        if ( dt <= 0.0 ) {
            DiffParam( time, d_x, d_y, D_X, D_Y );
            return;
        }

        const double t1 = std::max(0.0, time);
        const double t2 = t1 + dt;

        if ( DPROF == 0 ) {
            // Piecewise constant enhancement
            const double dt_enh_x = std::max(0.0, std::min(t2, tH0) - std::min(t1, tH0));
            const double dt_amb_x = dt - dt_enh_x;
            d_x = ( 1.13 * D_X * dt_enh_x + D_X * dt_amb_x ) / dt;

            const double dt_enh_y = std::max(0.0, std::min(t2, tV0) - std::min(t1, tV0));
            const double dt_amb_y = dt - dt_enh_y;
            d_y = ( 7.0 * D_Y * dt_enh_y + D_Y * dt_amb_y ) / dt;
        }
        else if ( DPROF == 1 ) {
            // Linear decay
            auto int_linear = [](double t, double D_0, double D_amb, double t_0) {
                if ( t <= 0.0 ) return 0.0;
                if ( t <= t_0 ) {
                    return D_0 * t + 0.5 * (D_amb - D_0) * (t * t) / t_0;
                } else {
                    double int_to_t0 = 0.5 * (D_0 + D_amb) * t_0;
                    return int_to_t0 + D_amb * (t - t_0);
                }
            };
            d_x = ( int_linear(t2, 1.13 * D_X, D_X, tH0) - int_linear(t1, 1.13 * D_X, D_X, tH0) ) / dt;
            d_y = ( int_linear(t2, 7.00 * D_Y, D_Y, tV0) - int_linear(t1, 7.00 * D_Y, D_Y, tV0) ) / dt;
        }
        else if ( DPROF == 2 ) {
            // Exponential decay
            auto avg_exp = [](double t_start, double t_end, double delta_t, double D_0, double D_amb, double t_0) {
                return D_amb + (D_0 - D_amb) * (t_0 / delta_t) * (std::exp(-t_start / t_0) - std::exp(-t_end / t_0));
            };
            d_x = avg_exp(t1, t2, dt, 1.13 * D_X, D_X, tH0);
            d_y = avg_exp(t1, t2, dt, 7.00 * D_Y, D_Y, tV0);
        }
        else {
            std::string const currFile("DiffParam.cpp");
            std::cout << "ERROR: In " << currFile << ": DPROF set to " << DPROF << "\n";
        }

        if ( d_x < 0.0 ) {
            std::cout << "d_x is negative: d_x = " << d_x << " [m^2/s]" << "\n";
            std::cout << "Setting d_x to 0.0" << "\n";
            d_x = 0.0;
        }

        if ( d_y < 0.0 ) {
            std::cout << "d_y is negative: d_y = " << d_y << " [m^2/s]" << "\n";
            std::cout << "Setting d_y to 0.0" << "\n";
            d_y = 0.0;
        }

    } /* End of DiffParam (time, dt) */

    void DiffParam( const double time, double &d_x, double &d_y, \
                    const double D_X, const double D_Y )
    {

        /* DiffParam:
        * Computes instantaneous diffusion parameters at time
        *
        * INPUTS:
        * (double) time: current time since simulation started in [s]
        *
        * OUTPUTS:
        * (double) d_x: current horizontal diffusion coefficient at time in [m^2/s]
        * (double) d_y: current vertical diffusion coefficient at time in [m^2/s]
        */

        if ( DPROF == 0 ) {
            if ( time <= tH0 )
                d_x = 1.13 * D_X;
            else
                d_x = D_X;
            if ( time <= tV0 )
                d_y = 7.0 * D_Y;
            else
                d_y = D_Y;
        }
        else if ( DPROF == 1 ) {
            d_x = std::max( D_X, 1.13 * D_X + (D_X - 1.13 * D_X) * time / tH0 );
            d_y = std::max( D_Y, 7.00 * D_Y + (D_Y - 7.00 * D_Y) * time / tV0 );
        }
        else if ( DPROF == 2 ) {
            d_x = D_X + (1.13 * D_X - D_X) * exp( -time / tH0 );
            d_y = D_Y + (7.00 * D_Y - D_Y) * exp( -time / tV0 );
        }
        else {
            std::string const currFile("DiffParam.cpp");
            std::cout << "ERROR: In " << currFile << ": DPROF set to " << DPROF << "\n";
        }

        if ( d_x < 0.0 ) {
            std::cout << "d_x is negative: d_x = " << d_x << " [m^2/s]" << "\n";
            std::cout << "Setting d_x to 0.0" << "\n";
            d_x = 0.0;
        }

        if ( d_y < 0.0 ) {
            std::cout << "d_y is negative: d_y = " << d_y << " [m^2/s]" << "\n";
            std::cout << "Setting d_y to 0.0" << "\n";
            d_y = 0.0;
        }

    } /* End of DiffParam (time) */

}