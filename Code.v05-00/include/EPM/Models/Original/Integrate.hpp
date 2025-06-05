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

#ifndef EPM_MODELS_ORIGINAL_INTEGRATE_H_INCLUDED
#define EPM_MODELS_ORIGINAL_INTEGRATE_H_INCLUDED

#include <variant>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include "Util/ForwardDecl.hpp"
#include "Core/Aircraft.hpp"
#include "Core/Emission.hpp"
#include "Core/Status.hpp"
#include "AIM/Coagulation.hpp"
#include "AIM/Aerosol.hpp"
#include "EPM/EPM.hpp"

namespace EPM::Models::OriginalImpl
{
    std::variant<EPM::Output, SimStatus> Integrate(
        double tempInit_K, double pressure_Pa, double rhw,
        double bypassArea, double coreExitTemp,
        double varArray[], const Vector_2D& aerArray,
        const Aircraft& AC,const Emission& EI, bool CHEMISTRY,
        double ambientLapseRate, std::string micro_data_out);

    SimStatus RunMicrophysics(
        double &temperature_K, double pressure_Pa, double relHumidity_w,
        double varArray[], const Vector_2D& aerArray,
        const Aircraft &AC, const Emission &EI, double delta_T_ad, double delta_T,
        double &Ice_rad, double &Ice_den, double &Soot_den, double &H2O_mol,
        double &SO4g_mol, double &SO4l_mol, AIM::Aerosol &SO4Aer, AIM::Aerosol &IceAer,
        double &Area, double &Ab0, double &Tc0, const bool CHEMISTRY,
        double ambientLapseRate, std::string micro_data_out);
}

#endif

