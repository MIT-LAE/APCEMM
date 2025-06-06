/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*                      Early Plume Microphysics                    */
/*                              (EPM)                               */
/*                                                                  */
/* Right-hand side Program File                                     */
/*                                                                  */
/* File                 : RHS.cpp                                   */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "AIM/Nucleation.hpp"
#include "Util/PhysConstant.hpp"
#include "Util/PhysFunction.hpp"
#include "EPM/Utils/Physics.hpp"
#include "EPM/Models/Original/Indexes.hpp"
#include "EPM/Models/Original/RHS.hpp"


using physConst::Na, physConst::R, physConst::kB, physConst::PI, physConst::RHO_ICE;

using physFunc::pSat_H2SO4, physFunc::pSat_H2Ol, physFunc::thermalSpeed;

using namespace EPM::Utils::Physics;
using namespace EPM::Models::OriginalImpl;


namespace EPM::Models::OriginalImpl {

    void gas_aerosol_rhs::operator()(
        const Vector_1D &x, Vector_1D &dxdt, const double t) const {
        /* DESCRIPTION:
            * Returns the right hand side of the system of ordinary differential equations:
            * y' = f(t,y)
            * where y belongs to an n-dimensional Euclidean space and t is time.
            * If f doesn't explicitly depend on time, the system is autonomous and the parameter t can be dropped
            * Default t value is 0. */

        /* INPUT:
            * - Vector_1D x    :: vector of variables
            * - Vector_1D dxdt :: vector of rate of change
            *
            * (optional)
            * - double t   :: time expressed in s
            *
            * OUTPUT:
            */

        double entrainRate = entrainmentRate(t);

        /* Compute air number concentration to convert number
        concentration rates to mixing ratio rates */ 
        double n_air = Na * m_pressure_Pa/(R * x[EPM_ind_T] * 1.0e6);

        /* Ensure that particle radius is larger than bare soot radius and smaller than some threshold */
        double m_part_r;
        m_part_r = (x[EPM_ind_ParR] > m_part_r0) ? x[EPM_ind_ParR] : m_part_r0;
        m_part_r = (m_part_r        > 1.00E-04 ) ? 1.00E-04        : m_part_r;

        /* Parameters for binary homogeneous nucleation */
        double nucRate, x_SO4, nMolec;

        /* Compute SO4 mixing ratios */

        /* Liquid and gaseous sulfate concentration */
        double SO4_rl, SO4_g;

        /* SO4_g represents the gaseous molecular concentration of SO4
            * SO4_rl is the "ready to be liquid" molecular concentration of SO4. This is still gaseous SO4 because it is limited by kinetics.
            * Gaseous SO4 is in phase equilibrium. Ensure limitations! */

        SO4_g  = ( pSat_H2SO4( x[EPM_ind_T] ) / ( kB * x[EPM_ind_T] * 1.0E+06 ) > x[EPM_ind_SO4g] * n_air ) ? x[EPM_ind_SO4g] : pSat_H2SO4( x[EPM_ind_T] ) / ( kB * x[EPM_ind_T] * 1.0E+06 * n_air); 

        // SO4_rl = x[EPM_ind_SO4g] * x[EPM_ind_P] / ( kB * x[EPM_ind_T] * 1.0E+06 ) - SO4_g;
        SO4_rl = x[EPM_ind_SO4g] - SO4_g; 
        SO4_g  = ( SO4_g > 0.0 ) ? SO4_g : 0.0;
        SO4_rl = ( SO4_rl > 0.0 ) ? SO4_rl : 0.0;


        double RH_liquid = n_air * x[EPM_ind_H2O] * 1.0e6 * kB * x[EPM_ind_T] / pSat_H2Ol( x[EPM_ind_T] ) ; 
        if ( ( n_air * SO4_rl >= AIM::nThresh( x[EPM_ind_T], RH_liquid))) {

            /* Mole fraction of sulfuric acid */
            x_SO4 = AIM::x_star(x[EPM_ind_T], RH_liquid, n_air*SO4_rl); 

            /* Number of molecules per cluster */
            nMolec = AIM::nTot(      x[EPM_ind_T], x_SO4, RH_liquid, n_air*SO4_rl );

            /* Nucleation rate */
            nucRate = AIM::nuclRate( x[EPM_ind_T], x_SO4, RH_liquid, n_air*SO4_rl );

        } else {

            x_SO4 = 0.0;
            nMolec = 0.0;
            nucRate = 0.0;

        }

    /* Compute sulfate - soot coagulation rate */
        double CoagRate = 0;
        Vector_1D pdf = nPDF_SO4.getPDF();
        Vector_1D binCenters = nPDF_SO4.getBinCenters();
        for ( unsigned int iBin = 0; iBin < nPDF_SO4.getNBin(); iBin++ ) {
            CoagRate += ( KernelSO4Soot[iBin] * pdf[iBin] * PI * binCenters[iBin] * binCenters[iBin] ) * ( 1.0 - x[EPM_ind_the1] - x[EPM_ind_the2] );
            /* Unit check:
                * [cm^3/s] * [#/cm3] * [m^2] = [m^2/s]
                */
        }


        /* Tracer differential equation
            * \frac{d[.]}{dt} = - w_T(t) * [.] */
        dxdt[EPM_ind_Trac] = - entrainRate * x[EPM_ind_Trac];

        /* Temperature differential equation
            * \frac{dT}{dt} = v_z(t) * \frac{dT}{dz} - w_T(t) * ( T - T_amb(t) ) */
        dxdt[EPM_ind_T] = dT_Vortex( t, m_delta_T, 1 ) \
                        - entrainRate * ( x[EPM_ind_T] - ( m_temperature_K + dT_Vortex ( t, m_delta_T ) ) );
        /* Unit check:
            * dT_Vortex(.,.,1): [K/s]
            * dT_Vortex(.,.,0): [K]
            * entrainRate * ...  : [1/s] * ( [K] - ( [K] + [K] ) ) = [K/s] */

        /* Pressure differential equation
            * \frac{dP}{dt} = 0 */
        dxdt[EPM_ind_P] = 0.0;

        /* Gaseous water differential equation 
            * \frac{d[H2O]}{dt} = - w_T(t) * ( [H2O] - [H2O]_amb ) */
        dxdt[EPM_ind_H2O] = - entrainRate * ( x[EPM_ind_H2O] - m_H2O_mixingratio ) \
                            - ( isFreezable( m_part_r, x[EPM_ind_T], x[EPM_ind_H2O]*n_air, m_part_r0 ) \
                                * depositionRate( m_part_r, x[EPM_ind_T], x[EPM_ind_P], x[EPM_ind_H2O]*n_air, m_part_r0, x[EPM_ind_the1] + x[EPM_ind_the2] ) \
                            + condensationRate( m_part_r, x[EPM_ind_T], x[EPM_ind_P], x[EPM_ind_H2O]*n_air, x[EPM_ind_the1] + x[EPM_ind_the2] ) ) \
                                * x[EPM_ind_Part] * Na / MW_H2O \
                            -  nucRate * ( 1.0 - x_SO4 ) * nMolec / n_air;
        /* Unit check:
            * [1/s] * ( [molec/cm3] - [molec/cm^3] ) = [molec/cm^3/s]
            * ([-] * [kg/s] + [kg/s]) * [#/cm^3] * [molec/kg] = [molec/cm^3/s]
            * */

        /* Liquid SO4 differential equation */
        dxdt[EPM_ind_SO4l] = - entrainRate * ( x[EPM_ind_SO4l] - m_SO4l_mixingratio ) \
                                + nucRate * x_SO4 * nMolec / n_air;

        /* Gaseous SO4 differential equation */
        dxdt[EPM_ind_SO4g] = - entrainRate * ( x[EPM_ind_SO4g] - m_SO4g_mixingratio ) \
                                - sticking_SO4 * thermalSpeed( x[EPM_ind_T], MW_H2SO4 / Na ) \
                                * x[EPM_ind_Part] * n_air * 1.0E+06 * PI * m_part_r * m_part_r \
                                * x[EPM_ind_SO4g] * ( 1.0 - x[EPM_ind_the1] - x[EPM_ind_the2] ) \
                                - nucRate * x_SO4 * nMolec / n_air;

        /* On soot SO4 differential equation 
            * \frac{d[SO4_s]}{dt} = alpha * v_th / 4.0 * n_part * 4.0 * \pi * r^2 * ( 1.0 - \theta ) * [SO4] */ 
        dxdt[EPM_ind_SO4s] = - entrainRate * ( x[EPM_ind_SO4s] - 0.0 ) \
                                + sticking_SO4 * thermalSpeed( x[EPM_ind_T], MW_H2SO4 / Na ) \
                                * x[EPM_ind_Part] * n_air * 1.0E+06 * PI * m_part_r * m_part_r \
                                * x[EPM_ind_SO4g] * ( 1.0 - x[EPM_ind_the1] - x[EPM_ind_the2] );

        /* Total SO4 differential equation
            * \frac{d[SO4]}{dt} = - w_T(t) * ( [SO4] - [SO4]_amb ) */
        dxdt[EPM_ind_SO4]  = dxdt[EPM_ind_SO4g] + dxdt[EPM_ind_SO4l] + dxdt[EPM_ind_SO4s];

        /* Gaseous HNO3 differential equation 
            * \frac{d[HNO3]}{dt} = - w_T(t) * ( [HNO3] - [HNO3]_amb ) */
        dxdt[EPM_ind_HNO3] = - entrainRate * ( x[EPM_ind_HNO3] - m_HNO3_mixingratio );

        /* Particle differential equation 
            * \frac{d[part]}{dt} = - w_T(t) * ( [part] - [part]_amb )*/
        dxdt[EPM_ind_Part] = - entrainRate * ( x[EPM_ind_Part] - m_part_mixingratio );

        /* Particle radius differential equation
            * \frac{dr}{dt} = \frac{dm}{dt}/(\rho * 4 * \pi * r^2) */
        dxdt[EPM_ind_ParR] = ( isFreezable( m_part_r, x[EPM_ind_T], n_air*x[EPM_ind_H2O], m_part_r0 ) \
                                * depositionRate( m_part_r, x[EPM_ind_T], x[EPM_ind_P], x[EPM_ind_H2O] * n_air , m_part_r0, x[EPM_ind_the1] + x[EPM_ind_the2] ) \
                                + condensationRate( m_part_r, x[EPM_ind_T], x[EPM_ind_P], x[EPM_ind_H2O] * n_air , x[EPM_ind_the1] + x[EPM_ind_the2] ) ) \
                                / ( RHO_ICE * 4.0 * PI * m_part_r * m_part_r );
        /* Unit check: 
            * [kg/s] / ( [kg/m^3] * [m^2]) = [m/s] */

        /*Diff Eqs for theta1 (fractional soot coverage due to adsorption) and theta2 (due to scavenging)
        CoagRate and d(theta1)/dt and d(theta2)/dt from Karcher (1998)*/
        dxdt[EPM_ind_the1] = sticking_SO4 * thermalSpeed( x[EPM_ind_T], MW_H2SO4 / Na ) * 1.0E+02 * 0.25 * n_air * ( SO4_g + SO4_rl ) / sigma_SO4 * ( 1.0 - x[EPM_ind_the1] - x[EPM_ind_the2] );

        dxdt[EPM_ind_the2] = CoagRate / ( 4.0 * PI * m_part_r0 * m_part_r0 );

    } /* End of gas_aerosol_rhs::operator() */

};
