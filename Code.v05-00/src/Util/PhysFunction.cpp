/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* PhysFunction Program File                                        */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : PhysFunction.cpp                          */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "Util/PhysFunction.hpp"

namespace physFunc
{

    double pSat_H2Ol( const double T )
    {
        
        /* DESCRIPTION:
         * Returns water liquid saturation pressure in Pascal.
         * Source: Sonntag (1994) */

        /* INPUT PARAMETERS:
         * - double T :: temperature expressed in K 
         *
         * OUTPUT PARAMETERS:
         * - double :: H2O liquid saturation pressure in Pascal */

        return 100.0 * exp( - 6096.9385 / T \
                            + 16.635794 \
                            - 0.02711193 * T \
                            + 1.673952E-5 * T * T \
                            + 2.433502 * log( T ) );

    } /* End of pSat_H2Ol */

    double pSat_H2Os( const double T )
    {
        
        /* DESCRIPTION:
         * Returns water solid saturation pressure in Pascal.
         * Source: Sonntag (1990) */

        /* INPUT PARAMETERS:
         * - double T :: temperature expressed in K 
         *
         * OUTPUT PARAMETERS:
         * - double :: H2O solid saturation pressure in Pascal */
        
        return 100.0 * exp( - 6024.5282 / T \
                            + 24.7219 \
                            + 0.010613868 * T \
                            - 1.3198825E-5 * T * T \
                            - 0.49382577 * log( T ) );

    } /* End of pSat_H2Os */

    double dpSat_H2Os( const double T )
    {

        /* DESCRIPTION:
         * Returns the derivative of the water solid saturation pressure with 
         * respect to temperature in Pascal/K. */

        /* INPUT PARAMETERS:
         * - double T :: temperature expressed in K 
         *
         * OUTPUT PARAMETERS:
         * - double :: dPsat/dT(T) in Pa/K */

        /* pSat(T) = C * exp( g )
         * dpSat(T)/dT = pSat(T) * g' */

        return pSat_H2Os( T ) * ( 6024.5282 / ( T * T ) + 0.010613868 - \
                                  2.0 * 1.3198825E-5 * T - 0.49382577 / T );

    } /* End of dpSat_H2Os */

    double pSat_H2SO4( const double T )
    {
        
        /* DESCRIPTION:
         * Returns sulfate liquid saturation pressure in Pascal.
         * Source: Noppel et al. (2002) */

        /* INPUT PARAMETERS:
         * - double T :: temperature expressed in K 
         *
         * OUTPUT PARAMETERS:
         * - double :: H2SO4 liquid saturation pressure in Pascal */
        
        return physConst::ATM * exp(-11.94 \
                            + 10156.0 * ( 1.0 / 360.15 \
                                          - 1.0 / T\
                                          + 0.38/(545.0)*( \
                                                           + 1.0 \
                                                           + log( 360.15 * 1.0 / T ) \
                                                           - 360.15 * 1.0 / T ) ) );

    } /* End of pSat_H2SO4 */

    double pSat_HNO3( const double T , const double PPH2O ) 
    {
        
        /* DESCRIPTION:
         * Returns nitric acid saturation pressure in Pascal.
         * Source: Hanson and Mauersberger (1988) */

        /* INPUT PARAMETERS:
         * - double T     :: temperature expressed in K 
         * - double PPH2O :: H2O partial pressure in Pa
         *
         * OUTPUT PARAMETERS:
         * - double :: HNO3 saturation pressure in Pascal */
        
        return ( physConst::ATM / 760.0 ) * pow( 10.0, ( ( ( - 2.7836 - 0.00088 * T ) * log10( PPH2O * ( 760.0 / physConst::ATM ) ) ) + ( 38.9855 - 11397.0 / T + 0.009179 * T ) ) );  

    } /* End of pSat_HNO3 */

    double rhoAir( const double T, const double P )
    {

        /* DESCRIPTION:
         * Returns the density of air in kg/m^3 */

        /* INPUT PARAMETERS:
         * - double T :: temperature in K
         * - double P :: pressure in Pa
         *
         * OUTPUT PARAMETERS:
         * - double :: Air density in kg/m^3 */

        return P / ( physConst::R_Air * T );

    } /* End of airDens */ 

    double dynVisc( const double T )
    {

        /* DESCRIPTION:
         * Returns the dynamic viscosity of air in kg/m/s */

        /* INPUT PARAMETERS:
         * - double T :: temperature in K
         *
         * OUTPUT PARAMETERS:
         * - double :: Dynamic viscosity in kg/m/s */

        return 1.8325E-05 * ( 416.16 / ( T + 120.0 ) ) * pow( T / 296.16, 1.5 );

    } /* End of dynVisc */

    double kinVisc( const double T, const double P )
    {

        /* DESCRIPTION:
         * Returns the kinematic viscosity of air in m^2/s */

        /* INPUT PARAMETERS:
         * - double T :: temperature in K
         * - double P :: pressure in Pa
         *
         * OUTPUT PARAMETERS:
         * - double :: Kinematic viscosity in m^2/s */

        return dynVisc( T ) / rhoAir( T, P );

    } /* End of kinVisc */

    double thermalSpeed( const double T, const double m )
    {

        /* DESCRIPTION:
         * Returns the thermal speed of a molecule/particle in m/s */

        /* INPUT PARAMETERS:
         * - double T :: temperature in K
         * - double m :: mass of the molecule/particle in kg ( set to the mass
         *               of an air molecule by default )
         *
         * OUTPUT PARAMETERS:
         * - double :: Thermal speed in m/s */

        return sqrt( 8.0 * physConst::kB * T / ( physConst::PI * m ) );

    } /* End of thermalSpeed */

    double lambda( const double T, const double P )
    {

        /* DESCRIPTION:
         * Returns the mean free path of an air molecule in m */

        /* INPUT PARAMETERS:
         * - double T :: temperature in K
         * - double P :: pressure in Pa
         *
         * OUTPUT PARAMETERS:
         * - double :: Mean free path in m */

        return 2.0 * kinVisc( T, P ) / thermalSpeed( T );

    } /* End of lambda */

    double mass_sphere( const double r, const double rho )
    {
    
        /* DESCRIPTION:
         * Returns the mass of a sphere in kg */

        /* INPUT PARAMETERS:
         * - double r   :: radius of the spherical particle in m
         * - double rho :: density of the particle in kg/m^3
         *
         * OUTPUT PARAMETERS:
         * - double :: Mass in kg */

        return 4.0 / double(3.0) * physConst::PI * rho * r * r * r;

    } /* End of mass_sphere */

    double vFall( const double r, const double rho, \
                      const double T, const double P )
    {

        /* DESCRIPTION:
         * Returns the terminal fall speed of a spherical particle in m/s */

        /* INPUT PARAMETERS:
         * - double r   :: radius of the spherical particle in m
         * - double rho :: density of the particle in kg/m^3
         * - double T   :: temperature in K
         * - double P   :: pressure in Pa
         *
         * OUTPUT PARAMETERS:
         * - double :: Terminal fall speed in m/s */

        return 2.0 * rho * physConst::g * r * r / ( 9.0 * dynVisc( T ) ) * slip_flowCorrection( Kn( r, T, P ) ); 

    } /* End of vFall */
    
    double Kn( const double r, const double T, \
                   const double P )
    {

        /* DESCRIPTION:
         * Returns the dimensionless particle Knudsen number */

        /* INPUT PARAMETERS:
         * - double r :: radius of the spherical particle in m
         * - double T :: temperature expressed in K
         * - double P :: pressure expressed in Pa
         *
         * OUTPUT PARAMETERS:
         * - double :: Particle Knudsen number */

        return lambda( T, P ) / r;

    } /* End of Kn */

    double partDiffCoef( const double r, const double T, \
                             const double P )
    {

        /* DESCRIPTION:
         * Returns the particle diffusion coefficient in m^2/s 
         * Source: Jacobson (2005) Equation 15.29*/

        /* INPUT PARAMETERS:
         * - double r :: radius of the spherical particle in m
         * - double T :: temperature expressed in K
         * - double P :: pressure expressed in Pa
         *
         * OUTPUT PARAMETERS:
         * - double :: Particle diffusion coefficient in m^2/s */

        return physConst::kB * T / ( 6.0 * physConst::PI * dynVisc( T ) * r ) * slip_flowCorrection( Kn ( r, T, P ) );

    } /* End of partDiffCoef */

    double slip_flowCorrection( const double Kn )
    {

        /* DESCRIPTION:
         * Returns the dimensionless Cunningham slip-flow correction.
         * Source: Jacobson (2005) page 509 */

        /* INPUT PARAMETERS:
         * - double Kn :: particle Knudsen number
         *
         * OUTPUT PARAMETERS:
         * - double :: slip-flow correction factor */

        static const double A = 1.249;
        static const double B = 0.42;
        static const double C = 0.87;

        return 1 + Kn * ( A + B * exp( - C / Kn ) );

    } /* End of slip_flowCorrection */
    
    double lambda_p( const double r, const double m, \
                         const double T, const double P )
    {

        /* DESCRIPTION:
         * Returns the mean free path of a particle in air in m 
         * Source: Jacobson (2005) Equation 15.34 */

        /* INPUT PARAMETERS:
         * - double r :: particle radius in m
         * - double m :: particle mass in kg
         * - double T :: temperature in K
         * - double P :: pressure in Pa
         *
         * OUTPUT PARAMETERS:
         * - double :: Particle mean free path in m */

        return 8.0 * partDiffCoef( r, T, P ) / ( physConst::PI * thermalSpeed( T, m ) );

    } /* End of lambda_p */

    double delta_p( const double r, const double m, \
                        const double T, const double P )
    {

        /* DESCRIPTION:
         * Returns the mean distance in m from the center of a sphere reached by particles 
         * leaving the sphere's surface and traveling a distance of particle mean free 
         * path lambda_p
         * Source: Jacobson (2005) Equation */

        /* INPUT PARAMETERS:
         * - double r :: particle radius in m
         * - double m :: particle mass in kg
         * - double T :: temperature in K
         * - double P :: pressure in Pa
         *
         * OUTPUT PARAMETERS:
         * - double :: Mean distance from the center reached by particles leaving the surface in m */

        double l = lambda_p( r, m, T, P );

        return ( pow( 2.0 * r + l , 3.0 ) - pow( 4.0 * r * r + l * l , 1.5 ) ) / ( 6.0 * r * l ) - 2.0 * r;

    } /* End of delta_p */

    double Reynolds_p( const double r, const double rho, \
                           const double T, const double P )
    {
        
        /* DESCRIPTION:
         * Returns the particle Reynolds number 
         * Ratio of the inertial force exerted by a particle to the viscous force 
         * exerted by the air.
         * Returns the mean free path of an air molecule in m
         *
         * If the inertial force is due to the particles's fall speed, the dimensionless Reynolds number is:
         * Re = 2 * r * V_t( r, rho ) / kinVisc( T, P ); */

        /* INPUT PARAMETERS:
         * - double r   :: particle radius in m
         * - double rho :: particle density in kg/m^3
         * - double T   :: temperature in K
         * - double P   :: pressure in Pa
         *
         * OUTPUT PARAMETERS:
         * - double :: Particle Reynolds number */

        return 2 * r * vFall( r, rho, T, P ) / kinVisc( T, P );

    } /* End of Reynolds_p */

    double Schmidt_p( const double r, const double T, \
                          const double P )
    {

        /* DESCRIPTION:
         * Returns the particle Schmidt number
         * Ratio of viscous to diffusive forces
         * 
         * The dimensionless Schmidt number is defined as:
         * Sc = kinVisc( T, P ) / partDiffCoef( r, T, P ); */

        /* INPUT PARAMETERS:
         * - double r   :: particle radius in m
         * - double T   :: temperature in K
         * - double P   :: pressure in Pa
         *
         * OUTPUT PARAMETERS:
         * - double :: Particle Schmidt number */

        return kinVisc( T, P )/ partDiffCoef( r, T, P );

    } /* End of Schmidt_p */

    double Stokes_p( const double r_1, const double rho_1, \
                         const double r_2, const double rho_2, \
                         const double T, const double P )
    {

        /* DESCRIPTION:
         * Returns the dimensionless particle Stokes number 
         * 
         * The particle Stokes number is defined as:
         * St = vFall( r_m, rho_m, T, P ) * | vFall( r_M, rho_M, T, P ) - vFall( r_m, rho_m, T, P ) | / ( r_M * g )
         * where r_M = max( r_1, r_2 ) and r_m = min( r_1, r_2 ); */
         
        /* INPUT PARAMETERS:
         * - double r_1   :: particle 1 radius in m
         * - double rho_1 :: particle 1 density in kg/m^3
         * - double r_2   :: particle 2 radius in m
         * - double rho_2 :: particle 2 density in kg/m^3
         * - double T     :: temperature in K
         * - double P     :: pressure in Pa
         *
         * OUTPUT PARAMETERS:
         * - double :: Particle Stokes number */

        if ( r_1 >= r_2 )
            return vFall( r_2, rho_2, T, P ) * std::abs( vFall( r_1, rho_1, T, P ) - vFall( r_2, rho_2, T, P ) ) / ( r_1 * physConst::g );
        else
            return vFall( r_1, rho_1, T, P ) * std::abs( vFall( r_2, rho_2, T, P ) - vFall( r_1, rho_1, T, P ) ) / ( r_2 * physConst::g );

    } /* End of Stokes_p */

    double E_agg( const double r_1, const double rho_1, \
                      const double r_2, const double rho_2, \
                      const double T, const double P )
    {

        /* DESCRIPTION:
         * Returns the aggregation efficiency for liquid particles
         * Data for aggregation efficiency under cirrus-like conditions are
         * not yet available in the literature (Sölch and Kärcher, 2010) 
         *
         * The aggregation efficiency for liquid particles is:
         * E_agg = ( 60.0 * E_V(r_m, rho_m, r_M, rho_M, T, P) + E_A(r_m, rho_m, r_M, rho_M, T, P) * Re(r_M, rho_M, T, P)) /\
         *         ( 60.0 + Re(r_M, rho_M, T, P));
         *
         * where r_M = max(r_1,r_2) and r_m = min(r_1,r_2) and:
         * E_V(r_1, rho_1, r_2, rho_2, T, P) = \
         *              0; ( if Stokes_p(r_1, r_2) <= 1.214 )
         *              pow( 1.0 + 0.75 * log( 2 * Stokes_p(r_1, rho_1, r_2, rho_2, T, P)) / \
         *                                (Stokes_p(r_1, rho_1, r_2, rho_2, T, P) - 1.214), -2.0 );  (else)
         * E_A(r_1, rho_1, r_2, rho_2, T, P) = pow( Stokes_p(r_1, rho_1, r_2, rho_2, T, P), 2.0 ) /\
         *              pow( Stokes_p(r_1, rho_1, r_2, rho_2, T, P) + 0.5, 2.0 );
         * (M.Z. Jacobson, Fundamentals of Atmospheric Modeling, 2005) */
         
        /* INPUT PARAMETERS:
         * - double r_1   :: particle 1 radius in m
         * - double rho_1 :: particle 1 density in kg/m^3
         * - double r_2   :: particle 2 radius in m
         * - double rho_2 :: particle 2 density in kg/m^3
         * - double T     :: temperature in K
         * - double P     :: pressure in Pa
         *
         * OUTPUT PARAMETERS:
         * - double :: Liquid particles aggregation efficiency */

        double s;
        double E_V, E_A;

        if ( r_1 >= r_2 ) {
            /* r_M = r_1 and r_m = r_2 */
            s = Stokes_p( r_2, rho_2, r_1, rho_1, T, P );
            if ( s > 1.214 ) {
                E_V = pow(1 + 0.75*log(2*s)/(s - 1.214), -2);
            } else {
                E_V = 0.0;
            }
            E_A = s * s / ( ( s + 0.5 ) * ( s + 0.5 ) );
            return ( 60.0 * E_V + E_A * Reynolds_p( r_1, rho_1, T, P ) )/( 60.0 + Reynolds_p( r_1, rho_1, T, P ) );
        } else {
            /* r_M = r_2 and r_m = r_1 */
            s = Stokes_p( r_1, rho_1, r_2, rho_2, T, P );
            if ( s > 1.214 ) {
                E_V = pow(1 + 0.75*log(2*s)/(s - 1.214), -2);
            } else {
                E_V = 0.0;
            }
            E_A = s * s / ( ( s + 0.5 ) * ( s + 0.5 ) );
            return ( 60.0 * E_V + E_A * Reynolds_p( r_2, rho_2, T, P ) )/( 60.0 + Reynolds_p( r_2, rho_2, T, P ) );
        }

    } /* End of E_agg */

    double DiffCoef_H2O( const double T, const double P )
    {

        /* DESCRIPTION: 
         * Returns the gas phase diffusion coefficient of water in m^2/s
         * Source: Pruppacher and Klett (2010), Equation 13.3 */

        /* INPUT PARAMETERS:
         * - double r :: particle radius expressed in m
         * - double T :: temperature expressed in K
         * - double P :: pressure expressed in Pa 
         *
         * OUTPUT PARAMETERS:
         * - double :: Water diffusion coefficient */

        return 2.1100E-05 * pow( ( T / physConst::TEMP_REF ), 1.94 ) * physConst::ATM / P;

    } /* End of DiffCoef_H2O */
    
    double DiffCoef_H2SO4( const double T, const double P )
    {

        /* DESCRIPTION: 
         * Returns the gas phase diffusion coefficient of sulfuric acid in m^2/s */

        /* INPUT PARAMETERS:
         * - double r :: particle radius expressed in m
         * - double T :: temperature expressed in K
         * - double P :: pressure expressed in Pa 
         *
         * OUTPUT PARAMETERS:
         * - double :: Sulfuric acid diffusion coefficient 
         *
         * (M.J. Tang, Gas phase diffusion coefficients: inorganic compounds, 
         * Atmospheric Chemistry and Physics, 2014)
         *
         * D(296K) = D(T) * (296.0/T)^(1.75) = 74.0 Torr cm^2 s^-1 */

        return 9.736842E-06 * pow( ( T / 296.0 ), 1.75 ) * physConst::ATM / P;

    } /* End of DiffCoef_H2SO4 */
    
    double DiffCoef_HNO3( const double T, const double P )
    {

        /* DESCRIPTION: 
         * Returns the gas phase diffusion coefficient of nitric acid in m^2/s */

        /* INPUT PARAMETERS:
         * - double r :: particle radius expressed in m
         * - double T :: temperature expressed in K
         * - double P :: pressure expressed in Pa 
         *
         * OUTPUT PARAMETERS:
         * - double :: Nitric acid diffusion coefficient 
         *
         * (M.J. Tang, Gas phase diffusion coefficients: inorganic compounds, 
         * Atmospheric Chemistry and Physics, 2014)
         *
         * D(296K) = D(T) * (296.0/T)^(1.75) = 87.0 Torr cm^2 s^-1 */

        return 1.144736E-05 * pow( ( T / 296.0 ), 1.75 ) * physConst::ATM / P;

    } /* End of DiffCoef_HNO3 */
    
    double CorrDiffCoef_H2O( const double r, const double T, \
                                 const double P )
    {

        /* DESCRIPTION: 
         * Returns the corrected gas phase diffusion coefficient of gaseous water in m^2/s
         * This correction accounts for the transition between the diffusion-limited 
         * ( r >> \lambda ) and the free molecular ( r << \lambda ) growth regime
         * Source: Pruppacher and Klett (2010) Equation 13.14 */

        /* INPUT PARAMETERS:
         * - double r :: particle radius expressed in m
         * - double T :: temperature expressed in K
         * - double P :: pressure expressed in Pa 
         *
         * OUTPUT PARAMETERS:
         * - double :: Corrected water diffusion coefficient */

        /* alpha represents the deposition coefficient for H2O molecules impinging on the 
         * surface. It is experimentally derived */

        static const double alpha = 0.036; //Pruppacher and Klett Table 5.5

        return DiffCoef_H2O( T, P ) \
            / ( r / ( r + lambda( T, P ) ) \
              + 4.0 * DiffCoef_H2O( T, P ) / ( alpha * r * thermalSpeed( T, MW_H2O / physConst::Na ) ) );

    } /* End of CorrDiffCoef_H2O */
    
    double CorrDiffCoef_H2SO4( const double r, const double T, \
                                   const double P )
    {

        /* DESCRIPTION: 
         * Returns the corrected gas phase diffusion coefficient of gaseous sulfuric acid in m^2/s
         * This correction accounts for the transition between the diffusion-limited 
         * ( r >> \lambda ) and the free molecular ( r << \lambda ) growth regime
         * Source: Pruppacher and Klett (2010) Equation 13.14 */

        /* INPUT PARAMETERS:
         * - double r :: particle radius expressed in m
         * - double T :: temperature expressed in K
         * - double P :: pressure expressed in Pa 
         *
         * OUTPUT PARAMETERS:
         * - double :: Corrected sulfuric acid diffusion coefficient */

        /* alpha represents the deposition coefficient for H2SO4 molecules impinging on the 
         * surface. It is experimentally derived */

        static const double alpha = 1.0;

        return DiffCoef_H2SO4( T, P ) \
            / ( r / ( r + lambda( T, P ) ) \
              + 4.0 * DiffCoef_H2SO4( T, P ) / ( alpha * r * thermalSpeed( T, MW_SO4 / physConst::Na ) ) );

    } /* End of CorrDiffCoef_H2SO4 */

    double CorrDiffCoef_HNO3( const double r, const double T, \
                                  const double P )
    {

        /* DESCRIPTION: 
         * Returns the corrected gas phase diffusion coefficient of gaseous nitric acid in m^2/s
         * This correction accounts for the transition between the diffusion-limited 
         * ( r >> \lambda ) and the free molecular ( r << \lambda ) growth regime.
         * Source: Pruppacher and Klett (2010) Equation 13.14 */

        /* INPUT PARAMETERS:
         * - double r :: particle radius expressed in m
         * - double T :: temperature expressed in K
         * - double P :: pressure expressed in Pa 
         *
         * OUTPUT PARAMETERS:
         * - double :: Corrected nitric acid diffusion coefficient */

        /* alpha represents the deposition coefficient for HNO3 molecules impinging on the 
         * surface. It is experimentally derived */

        static const double alpha = 1.0;

        return DiffCoef_HNO3( T, P ) \
            / ( r / ( r + lambda( T, P ) ) \
              + 4.0 * DiffCoef_HNO3( T, P ) / ( alpha * r * thermalSpeed( T, MW_HNO3 / physConst::Na ) ) );

    } /* End of CorrDiffCoef_HNO3 */

    double ThermalCond( const double r, const double T, \
                            const double P )
    {

        /* DESCRIPTION:
         * Returns the corrected thermal conductivity of dry air in J / (m s K)
         * through which heat is transported */

        /* INPUT PARAMETERS:
         * - double r :: particle radius expressed in m
         * - double T :: temperature expressed in K
         * - double P :: pressure expressed in Pa 
         *
         * OUTPUT PARAMETERS:
         * - double :: Thermal conductivity of dry air in J / ( m s K ) 
         *
         * (H.R. Pruppacher and J.D. Klett, Microphysics of Clouds and Precipitation,
         *  Kluwer Academic Publishers, 1997)*/

        /* alpha_T is experimentally derived */

        static const double k_a = 2.50E-02; /* [J / (m s K)] */
        static const double alpha_T = 0.7;

        return k_a \
            / ( r / ( r + 2.16E-07 ) \
              + 4.0 * k_a / ( alpha_T * r * 1000.* physConst::CP_Air * rhoAir( T, P ) * thermalSpeed( T, MW_Air / physConst::Na ) ) );

    } /* End of ThermalCond */

    double LHeatSubl_H2O( const double T )
    {

        /* DESCRIPTION:
         * Returns the latent heat of sublimation of water vapor in J/kg
         * Source: Polynomial curve fit to Table 2.1 of Rogers and Yau (1989)*/

        /* INPUT PARAMETERS:
         * - double T :: temperature in K
         *
         * OUTPUT PARAMETERS:
         * - double :: Latent heat of sublimation in J/kg*/

        return ( 2.8341E+06 - 2.90E+02 * ( T - physConst::TEMP_REF ) - 4.00E+00 * ( T - physConst::TEMP_REF ) * ( T - physConst::TEMP_REF ) );

    } /* End of LHeatSubl_H2O */

    double Kelvin( const double r )
    {
        
        /* DESCRIPTION:
         * Returns the saturation pressure correction factor to account for the Kelvin effect */

        /* INPUT PARAMETERS:
         * - double r :: particle radius expressed in m
         *
         * OUTPUT PARAMETERS:
         * - double :: Dimensionless Kelvin factor */

        /* The parameter a_k is taken from:
         * (J. Picot et al., Large-eddy simulation of contrail evolution in the vortex phase 
         * and its interaction with atmospheric turbulence, Atmospheric Chemistry and Physics, 2015)*/

        static const double a_k = 5.00E-10; /* [m] */

        return exp( a_k / r );

    } /* End of Kelvin */

    Vector_2D RHi_Field(const Vector_2D& H2O, const Vector_2D& T, const Vector_1D& P) {
        Vector_2D RHi(H2O.size(), Vector_1D(H2O[0].size(), 0));
        double pH2O;
        #pragma omp parallel for \
        private(pH2O)
        for(std::size_t jNy = 0; jNy < H2O.size(); jNy++){
            for (std::size_t iNx = 0; iNx < H2O[0].size(); iNx++){
                pH2O = physConst::R * T[jNy][iNx] * H2O[jNy][iNx] /(physConst::Na*1e-6);
                RHi[jNy][iNx] = pH2O / pSat_H2Os(T[jNy][iNx]);
            }
        }
        return RHi;
    } //End of RH_Field


}

/* End of PhysFunction.cpp */

