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

#include "EPM/RHS.hpp"
#include "EPM/Indexes.hpp"

using physConst::Na, physConst::R, physConst::kB, physConst::PI, physConst::RHO_ICE;

using physFunc::pSat_H2SO4, physFunc::pSat_H2Ol, physFunc::pSat_H2Os,
    physFunc::pSat_HNO3, physFunc::thermalSpeed, physFunc::CorrDiffCoef_H2O,
    physFunc::Kelvin, physFunc::LHeatSubl_H2O, physFunc::ThermalCond;

namespace EPM {

    /* Vortex sinking timescales, taken from Unterstrasser et al., 2008 */
    const double t_Vortex_0 = 8.00E+00;
    const double t_Vortex_1 = 1.10E+02;

    /* Dilution timescales for a B747, taken from:
     * B. Kärcher, "A trajectory box model for aircraft exhaust plumes", Journal of Geophysical Research, 1995 */
    const double t_0 = 1.00E-04; /* [s],  */
    const double t_j = 1.00E-02; /* [s],  */
    const double t_1 = 8.00E+00; /* [s], Transition to vortex regime */
    const double t_2 = 6.60E+01; /* [s], Transition to dispersion regime */

    const double m = 2.0;
    const double n = 50.0;
    const double Cv = 3.0;

    /* Engine exit plane characteristics for a B747, taken from:
     * B. Kärcher, "A trajectory box model for aircraft exhaust plumes", Journal of Geophysical Research, 1995 */
    /* Engine exit core area im m^2 */
    const double Ac0 = 0.604;
    /* Engine exit core velocity in m/s */
    const double uc0 = 475.7;
    /* Engine exit core temperature in K */
    /* const double Tc0 = 547.3; */
    /* double Tc0 */
    /* Engine exit bypass area in m^2 */
    /* const double Ab0 = 1.804; */
    /* double Ab0 */


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

    double entrainmentRate(const double t) {
        /* DESCRIPTION:
         * Returns the plume entrainment rate in the plume early stages.
         * Based on Karcher 1995 */

        /* INPUT:
         * - double t :: time since engine exit plane expressed in s
         *
         * OUTPUT:
         * - double :: plume entrainment rate 1/s */

        double dRat = 0;

        if ( ( t >= t_0 ) && ( t <= t_j ) ) {
            if ( t < 1.0E-03 ) {
                /* Obtained from best fit */
                dRat = exp( 12.175448395306601 * log(t) + 83.748342982715030 );
            } else if ( t >= 1.0E-03 ) { 
                dRat = std::max( 0.8 / t_j - 0.8 / ( t_j * t_j ) * ( t - t_j ) + ( 1 - 0.8 / t_j + 0.8 / ( t_j * t_j ) * ( 1.0E-03 - t_j ) ) / ( ( 1.0E-03 - t_j ) * ( 1.0E-03 - t_j ) ) * ( t - t_j ) * ( t - t_j ), 1.0E-10 );
            }
        } else if ( ( t > t_j ) && ( t <= t_1 ) ) {
            dRat = 0.8 / t ;
        } else if ( ( t > t_1 ) && ( t <= t_2 ) ) {
            dRat = ( 0.8 / t_1 ) / ( 1.0 + Cv * ( 0.8 / t_1 ) * ( t - t_1 ) );
        } else if ( ( t > t_2 ) ) {
            dRat = ( 0.8 / t_1 ) / ( 1.0 + Cv * ( 0.8 / t_1 ) * ( t_2 - t_1 ) ) / ( 1.0 + 1.0 / ( 1.0 + m * exp( - ( t - t_2 ) / ( n * t_2 ) ) ) * ( 0.8 / t_1 ) / ( 1.0 + Cv * ( 0.8 / t_1 ) * ( t_2 - t_1 ) ) * ( t - t_2 ) );
        } else {
            dRat = 0;
        }

        return dRat;
    } /* End of entrainmentRate */

    double depositionRate(const double r, const double T, const double P,
                          const double H2O, const double r_0,  const double theta) {
        /* DESCRIPTION:
         * Returns the water deposition rate in kg/s on a single spherical particle */

        /* INPUT:
         * - double r     :: radius in m
         * - double T     :: temperature in K
         * - double P     :: pressure in Pa
         * - double H2O   :: gaseous water concentration in molec/cm^3
         * - double r_0   :: bare soot radius in m
         * - double theta :: particle surface coverage in [-]
         *
         * OUTPUT:
         * - double :: deposition rate */

        /* Capacitance factor for spherical particles is C = r; */

        double p_h2o = H2O * kB * T * 1.0E+06; // Partial pressure of water

        // If supersaturated w.r.t. liquid 
        if ( ( p_h2o / pSat_H2Ol( T ) >= 1.0 ) && ( p_h2o / pSat_H2Ol( T ) < 2.0 ) ) {
            if ( ( r * r - r_0 * r_0 ) / ( r_0 * r_0 ) <= 0.5 ) {
                return 4.0 * PI * r * theta * CorrDiffCoef_H2O( r, T, P ) * \
                       ( p_h2o - pSat_H2Os( T ) * Kelvin( r ) ) / \
                       ( R / MW_H2O * T + LHeatSubl_H2O( T ) * CorrDiffCoef_H2O( r, T, P ) * \
                         pSat_H2Os( T ) * ( LHeatSubl_H2O( T ) * MW_H2O / ( R * T ) - 1.0 ) / \
                           ( ThermalCond( r, T, P ) * T ) );
            } else {
                return 4.0 * PI * r * CorrDiffCoef_H2O( r, T, P ) * \
                       ( p_h2o - pSat_H2Os( T ) * Kelvin( r ) ) / \
                       ( R / MW_H2O * T + LHeatSubl_H2O( T ) * CorrDiffCoef_H2O( r, T, P ) * \
                         pSat_H2Os( T ) * ( LHeatSubl_H2O( T ) * MW_H2O / ( R * T ) - 1.0 ) / \
                           ( ThermalCond( r, T, P ) * T ) );
            }
        } else { // Melting
            return 4.0 * PI * r * CorrDiffCoef_H2O( r, T, P ) * \
                   ( p_h2o  - pSat_H2Os( T ) * Kelvin( r ) ) / \
                   ( R / MW_H2O * T + LHeatSubl_H2O( T ) * CorrDiffCoef_H2O( r, T, P ) * \
                     pSat_H2Os( T ) * ( LHeatSubl_H2O( T ) * MW_H2O / ( R * T ) - 1.0 ) / \
                       ( ThermalCond( r, T, P ) * T ) );
        }
    } /* End of depositionRate */

    double dT_Vortex(const double time, const double delta_T, bool deriv) {
        /* DESCRIPTION:
         * Returns the temperature increase / rate of change from cruise altitude as a
         * function of time */

        /* INPUT:
         * - double time    :: time since engine exit plane expressed in s
         * - double delta_T :: temperature increase in K after vortex sinking in K
         * - bool deriv         :: rate of change or delta?
         *
         * OUTPUT:
         * - double :: temperature increase / rate of change from cruise altitude in K */

        if ( deriv == 0 ) {
            if ( ( time >= t_Vortex_0 ) )
                return std::min( delta_T, delta_T * ( time - t_Vortex_0 ) / ( t_Vortex_1 - t_Vortex_0 ) );
            else
                return 0.0;
        } else {
            if ( ( time >= t_Vortex_0 ) && ( time <= t_Vortex_1 ) )
                return delta_T / ( t_Vortex_1 - t_Vortex_0 );
            else
                return 0.0;
        }

    } /* End of dT_Vortex */

    bool isFreezable(const double r, const double T, const double H2O, const double r0) {
        /* DESCRIPTION:
         * Returns whether or not water can freeze on a particle */

        /* INPUT:
         * - double r   :: radius in m
         * - double T   :: temperature in K
         * - double P   :: pressure in Pa
         * - double H2O :: gaseous water concentration in molec/cm^3
         * - double r0  :: bare soot radius in m
         *
         * OUTPUT:
         * - bool :: isFreezable? */

        /* Returns 1 if:
         * - Current conditions are supersaturated with respect to liquid water (freezing)
         * or
         * - Current radius is greater than bare soot radius (melting) */

        return ( ( H2O * kB * T * 1.0E+06 >= Kelvin( r ) * pSat_H2Ol( T ) ) || ( r > 1.005 * r0 ) );

    } /* End of isFreezable */

    double condensationRate(const double r, const double T, const double P, const double H2O, const double theta) {
        /* DESCRIPTION:
         * Returns the gaseous water condensation rate in kg/s on a single spherical particle */

        /* INPUT:
         * - double r     :: radius in m
         * - double T     :: temperature in K
         * - double P     :: pressure in Pa
         * - double H2O   :: gaseous water concentration in molec/cm^3
         * - double theta :: fractional coating of the particle
         *
         * OUTPUT:
         * - double :: condensation rate */

        double p_h2o = H2O * kB * T * 1.0E+06; // Partial pressure of water

        if ( p_h2o >= Kelvin( r ) * pSat_H2Ol( T ) ) {
            return 4.0 * PI * MW_H2O * CorrDiffCoef_H2O( r, T, P ) * r * theta \
                       * ( H2O / Na * 1.0E+06 \
                         - Kelvin( r ) * pSat_H2Ol( T ) / ( R * T ) );
        } else {
            return 0.0;
        }
    } /* End of condensationRate */

};
