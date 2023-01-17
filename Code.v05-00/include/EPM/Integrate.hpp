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
#include "AIM/Coagulation.hpp"
#include "AIM/Nucleation.hpp"
#include "AIM/Aerosol.hpp"
#include "odeSolver.hpp"

namespace EPM
{

    static const int EPM_SUCCESS = 1;
    static const int EPM_FAILURE = 0;
    // Stopped early (no ice)
    static const int EPM_EARLY   = 2;

    /* Vortex sinking timescales, taken from Unterstrasser et al., 2008 */
    const RealDouble t_Vortex_0 = 8.00E+00;
    const RealDouble t_Vortex_1 = 1.10E+02;

    /* Dilution timescales for a B747, taken from:
     * B. Kärcher, "A trajectory box model for aircraft exhaust plumes", Journal of Geophysical Research, 1995 */
    const RealDouble t_0 = 1.00E-04; /* [s],  */
    const RealDouble t_j = 1.00E-02; /* [s],  */
    const RealDouble t_1 = 8.00E+00; /* [s], Transition to vortex regime */
    const RealDouble t_2 = 6.60E+01; /* [s], Transition to dispersion regime */
   
    const RealDouble m = 2.0;
    const RealDouble n = 50.0;
    const RealDouble Cv = 3.0;
    
    /* Engine exit plane characteristics for a B747, taken from:
     * B. Kärcher, "A trajectory box model for aircraft exhaust plumes", Journal of Geophysical Research, 1995 */
    /* Engine exit core area im m^2 */
    const RealDouble Ac0 = 0.604;
    /* Engine exit core velocity in m/s */
    const RealDouble uc0 = 475.7;
    /* Engine exit core temperature in K */
    /* const RealDouble Tc0 = 547.3; */
    /* RealDouble Tc0 */
    /* Engine exit bypass area in m^2 */
    /* const RealDouble Ab0 = 1.804; */
    /* RealDouble Ab0 */

    int Integrate( RealDouble &temperature_K, RealDouble pressure_Pa, RealDouble relHumidity_w, RealDouble varArray[], \
                   RealDouble fixArray[], RealDouble aerArray[][2], const Aircraft &AC, const Emission &EI, \
                   RealDouble &Ice_rad, RealDouble &Ice_den, RealDouble &Soot_den, RealDouble &H2O_mol, \
                   RealDouble &SO4g_mol, RealDouble &SO4l_mol, AIM::Aerosol &SO4Aer, AIM::Aerosol &IceAer, \
                   RealDouble &Area, RealDouble &Ab0, RealDouble &Tc0, const bool CHEMISTRY, std::string micro_data_out );
    int RunMicrophysics( RealDouble &temperature_K, RealDouble pressure_Pa, RealDouble relHumidity_w, \
                         RealDouble varArray[], RealDouble fixArray[], RealDouble aerArray[][2], \
                         const Aircraft &AC, const Emission &EI, RealDouble delta_T_ad, RealDouble delta_T, \
                         RealDouble &Ice_rad, RealDouble &Ice_den, RealDouble &Soot_den, RealDouble &H2O_mol, \
                         RealDouble &SO4g_mol, RealDouble &SO4l_mol, AIM::Aerosol &SO4Aer, AIM::Aerosol &IceAer, \
                         RealDouble &Area, RealDouble &Ab0, RealDouble &Tc0, const bool CHEMISTRY, std::string micro_data_out );
    RealDouble dT_Vortex( const RealDouble time, const RealDouble delta_T, bool deriv = 0 );
    RealDouble entrainmentRate( const RealDouble time );
    RealDouble depositionRate( const RealDouble r, const RealDouble T, const RealDouble P, const RealDouble H2O, \
                               const RealDouble r_0,  const RealDouble theta );
    void odeRHS( const Vector_1D &x, Vector_1D &dxdt, const RealDouble t = 0.0 );
    bool isFreezable( const RealDouble r, const RealDouble T, const RealDouble H2O, const RealDouble r0 );
    RealDouble condensationRate( const RealDouble r, const RealDouble T, const RealDouble P, const RealDouble H2O, \
                                 const RealDouble theta );


}


#endif /* INTEGRATE_H_INCLUDED */
