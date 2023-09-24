/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* PhysFunction Header File                                         */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : PhysFunction.hpp                          */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef PHYSFUNCTION_H_INCLUDED
#define PHYSFUNCTION_H_INCLUDED

#include <cmath>

#include "ForwardDecl.hpp"
#include "PhysConstant.hpp"

namespace physFunc
{

    /* Defined physical functions */

    /* Saturation pressures [Pa] */
    double pSat_H2Ol( const double T );
    double pSat_H2Os( const double T );
    double pSat_H2SO4( const double T );
    double pSat_HNO3( const double T , const double PPH2O );

    /* Derivative of the saturation pressures wrt to the temperature [Pa/K] */
    double dpSat_H2Os( const double T );

    /* Density of air [kg/m^3] */
    double rhoAir( const double T, const double P );

    /* Dynamic viscosity [kg/(m s)] */
    double dynVisc( const double T );

    /* Kinematic viscosity [m^2/s] */
    double kinVisc( const double T, const double P );

    /* Thermal speed of an air molecule/particle [m/s] */
    double thermalSpeed( const double T, \
                             const double m = physConst::M_Air );

    /* Mean free path in air [m] */
    double lambda( const double T, const double P );

    /* Mass of spherical particle [kg] */
    double mass_sphere( const double r, const double rho );
    
    /* Terminal fall speed of a spherical particle [m/s] */
    double vFall( const double r, const double rho, \
                      const double T, const double P );

    /* Knudsen number [-] */
    double Kn( const double r, const double T, \
                   const double P );

    /* Particle diffusion coefficient [m^2/s] */
    double partDiffCoef( const double r, const double T, \
                             const double P );

    /* Cunningham slip-flow correction factor [-] */
    double slip_flowCorrection( const double Kn );
    
    /* Mean free path of particles in air [m] */
    double lambda_p( const double r, const double m, \
                         const double T, const double P );
    
    /* Mean distance [m] from the center of a sphere reach by particles leaving */
    double delta_p( const double r, const double m, \
                        const double T, const double P );
    
    /* Particle Reynolds number [-] */
    double Reynolds_p( const double r, const double rho, \
                           const double T, const double P );
    
    /* Particle Schmidt number [-] */
    double Schmidt_p( const double r, const double T, \
                          const double P );

    /* Particle Stokes number [-] */
    double Stokes_p( const double r_1, const double rho_1, \
                         const double r_2, const double rho_2, \
                         const double T, const double P );
    
    /* Aggregation efficiency for liquid particles [-] */
    double E_agg( const double r_1, const double rho_1, \
                      const double r_2, const double rho_2, \
                      const double T, const double P );
    
    /* H2O gas phase diffusion coefficient in [m^2/s] */
    double DiffCoef_H2O( const double T, const double P );
    
    /* H2SO4 gas phase diffusion coefficient in [m^2/s] */
    double DiffCoef_H2SO4( const double T, const double P );
    
    /* HNO3 gas phase diffusion coefficient in [m^2/s] */
    double DiffCoef_HNO3( const double T, const double P );

    /* Corrected H2O gas phase diffusion coefficient in [m^2/s] */
    double CorrDiffCoef_H2O( const double r, const double T, \
                                 const double P );
    
    /* Corrected H2SO4 gas phase diffusion coefficient in [m^2/s] */
    double CorrDiffCoef_H2SO4( const double r, const double T, \
                                   const double P );

    /* Corrected HNO3 gas phase diffusion coefficient in [m^2/s] */
    double CorrDiffCoef_HNO3( const double r, const double T, \
                                  const double P );
    
    /* Thermal conductivity of dry air in [J/(msK)] */
    double ThermalCond( const double r, const double T, \
                            const double P );

    /* Latent heat of sublimation of water vapor in [J/kg] */
    double LHeatSubl_H2O( const double T );

    /* Kelvin factor [-] */
    double Kelvin( const double r );

    /* RH Field */
    Vector_2D RHi_Field(const Vector_2D& H2O, const Vector_2D& T, const Vector_1D& P);

    inline double RHwToRHi(double rhw, double temp) {
        return rhw * pSat_H2Ol( temp ) / pSat_H2Os( temp );
    }

    inline double RHiToRHw(double rhi, double temp) {
        return rhi * pSat_H2Os( temp ) / pSat_H2Ol( temp );
    }

    inline double H2OToRHw(double h2o, double temp) {
        //h2o is in molec / cm3
        return h2o / 1.0e-6 * physConst::kB * temp / pSat_H2Ol(temp) * 100;
    }

    inline double RHwToH2O(double rhw, double temp) {
        //returns h2o conc in molec / cm3
        return rhw / 100 * pSat_H2Ol(temp) / (physConst::kB * temp) * 1e-6;
    }

    inline double RHiToH2O(double rhi, double temp) {
        //returns h2o conc in molec / cm3
        return rhi / 100 * pSat_H2Os(temp) / (physConst::kB * temp) * 1e-6;
    }

}

#endif /* PHYSFUNCTION_H_INCLUDED */
