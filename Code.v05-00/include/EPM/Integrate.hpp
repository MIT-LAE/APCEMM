/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*                      Early Plume Microphysics                    */
/*                              (EPM)                               */
/*                                                                  */
/* Integrate Header File                                            */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 9/27/2018                                 */
/* File                 : Integrate.hpp                             */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef INTEGRATE_H_INCLUDED
#define INTEGRATE_H_INCLUDED

#include <iostream>
#include <cmath>
#include <vector>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include "omp.h"

#include "Util/ForwardDecl.hpp"
#include "Util/PhysConstant.hpp"
#include "Util/PhysFunction.hpp"
#include "Core/Interface.hpp"
#include "Core/Parameters.hpp"
#include "Core/Monitor.hpp"
#include "Core/Aircraft.hpp"
#include "Core/Emission.hpp"
#include "Core/Status.hpp"
#include "AIM/Coagulation.hpp"
#include "AIM/Nucleation.hpp"
#include "AIM/Aerosol.hpp"
#include "odeSolver.hpp"

namespace EPM
{
    struct EPMOutput {
        double finalTemp;
        double iceRadius;
        double iceDensity;
        double sootDensity;
        double H2O_mol;
        double SO4g_mol;
        double SO4l_mol;
        AIM::Aerosol SO4Aer;
        AIM::Aerosol IceAer;
        double area;
        double bypassArea;
        double coreExitTemp;
    };

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

    SimStatus Integrate( double &temperature_K, double pressure_Pa, double relHumidity_w, double varArray[], \
                   const Vector_2D& aerArray, const Aircraft &AC, const Emission &EI, \
                   double &Ice_rad, double &Ice_den, double &Soot_den, double &H2O_mol, \
                   double &SO4g_mol, double &SO4l_mol, AIM::Aerosol &SO4Aer, AIM::Aerosol &IceAer, \
                   double &Area, double &Ab0, double &Tc0, const bool CHEMISTRY, double ambientLapseRate, std::string micro_data_out );

    std::pair<EPMOutput, SimStatus> Integrate(double tempInit_K, double pressure_Pa, double rhw, double bypassArea, double coreExitTemp, double varArray[], 
                            const Vector_2D& aerArray, const Aircraft& AC,const Emission& EI, bool CHEMISTRY, double ambientLapseRate, std::string micro_data_out);

    SimStatus RunMicrophysics( double &temperature_K, double pressure_Pa, double relHumidity_w, \
                         double varArray[], const Vector_2D& aerArray, \
                         const Aircraft &AC, const Emission &EI, double delta_T_ad, double delta_T, \
                         double &Ice_rad, double &Ice_den, double &Soot_den, double &H2O_mol, \
                         double &SO4g_mol, double &SO4l_mol, AIM::Aerosol &SO4Aer, AIM::Aerosol &IceAer, \
                         double &Area, double &Ab0, double &Tc0, const bool CHEMISTRY, double ambientLapseRate, std::string micro_data_out );
    double dT_Vortex( const double time, const double delta_T, bool deriv = 0 );
    double entrainmentRate( const double time );
    double depositionRate( const double r, const double T, const double P, const double H2O, \
                               const double r_0,  const double theta );
    void odeRHS( const Vector_1D &x, Vector_1D &dxdt, const double t = 0.0 );
    bool isFreezable( const double r, const double T, const double H2O, const double r0 );
    double condensationRate( const double r, const double T, const double P, const double H2O, \
                                 const double theta );


}


#endif /* INTEGRATE_H_INCLUDED */
