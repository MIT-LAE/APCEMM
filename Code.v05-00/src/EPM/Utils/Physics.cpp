#include "Util/PhysConstant.hpp"
#include "Util/PhysFunction.hpp"
#include "EPM/Utils/Physics.hpp"

using physConst::Na, physConst::R, physConst::kB, physConst::PI;

using  physFunc::pSat_H2Ol, physFunc::pSat_H2Os, physFunc::CorrDiffCoef_H2O,
    physFunc::Kelvin, physFunc::LHeatSubl_H2O, physFunc::ThermalCond;

namespace EPM::Utils::Physics {

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

} // namespace EPM::Utils::Physics
