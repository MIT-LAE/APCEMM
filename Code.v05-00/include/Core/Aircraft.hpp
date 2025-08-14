/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Aircraft Header File                                             */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : Aircraft.hpp                              */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef AIRCRAFT_H_INCLUDED
#define AIRCRAFT_H_INCLUDED

#include <string>
#include <cstring>

#include "Core/Input.hpp"
#include "Core/Engine.hpp"
#include "Core/Vortex.hpp"
#include "Core/Meteorology.hpp"

class Aircraft 
{
    public:

        /* Constructors */

        Aircraft( );
        Aircraft( const Input& input, std::string engineFilePath, std::string engineName = "GEnx-2B67B");

        /* Debug */
        void Debug( ) const;

        /* Compute vortex losses */

        double VortexLosses( const double N_postjet, const double WV_exhaust, const double N0_ref );

        /* Getters: */

        /* Flight conditions */
        inline double T_CA_K() const { return T_CA_K_; }
        inline double RHW_CA_PC() const { return RHW_CA_PC_; }
        inline double RHi_CA_PC() const { return RHi_CA_PC_; }
        inline double nBV_Hz() const { return nBV_Hz_; }
        inline double p_CA_Pa() const { return p_CA_Pa_; }
        /* Aircraft name */
        inline std::string Name() const { return Name_; }
        /* Flight velocity */
        inline double VFlight() const { return vFlight_ms_; }
        /* Mach number */
        inline double Mach() const { return machNumber_; }
        /* Fuel consumption per distance */
        inline double fuel_per_dist() const { return fuel_per_dist_; }
        /* Wingspan */
        inline double Wingspan() const { return wingspan_; }
        /* Current mass */
        inline double currMass() const { return currMass_; }
        /* Fuel flow */
        inline double FuelFlow() const { return engine_.getFuelFlow() * engNumber_; }
        /* Engine number */
        inline UInt EngNumber() const { return engNumber_; }
        /* Engine */
        inline const Engine& engine() const { return engine_; }
        /* Vortex */
        inline const Vortex& vortex() const { return vortex_; }

        /* Setters for engine properties */
        void setEI_NOx(const double NOx);
        void setEI_CO(const double CO);
        void setEI_HC(const double HC);
        void setEI_Soot(const double Soot);
        void setSootRad(const double sootRad);
        void setFuelFlow(const double ff);
        void setVFlight(const double Vf, double temperature_K);
        void setEngNumber(const double nEng);
        void setWingspan(const double span);
        void setMass(double mass);

        /* Engine */
        Engine engine_;

    protected:

        /* Flight conditions */
        double T_CA_K_;
        double RHW_CA_PC_;
        double RHi_CA_PC_;
        double nBV_Hz_;
        double p_CA_Pa_;

        /* Aircraft name */
        std::string Name_;

        /* Flight speed & mach Number */
        double vFlight_ms_;
        double machNumber_;

        /* Fuel consumption */
        double fuel_per_dist_; /* [kg/m] */

        /* Dimensions */
        double wingspan_; /* [m] */

        /* Weight */
        double currMass_; /* [kg] */

        /* Number of engines */
        UInt engNumber_;

        /* Vortex */
        Vortex vortex_;

    private:

};

#endif /* AIRCRAFT_H_INCLUDED */
