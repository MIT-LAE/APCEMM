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

class Aircraft 
{
    public:

        /* Constructors */

        Aircraft( );
        Aircraft( const char *aircraftName, std::string engineFilePath, double aircraftMass, \
                  double temperature_K, double pressure_Pa,  \
                  double relHumidity_w, double nBV );
        Aircraft( const Input& input, std::string engineFilePath, std::string engineName = "GEnx-2B67B");

        /* Debug */
        void Debug( ) const;

        /* Compute vortex losses */

        double VortexLosses( const double EI_Soot, const double EI_SootRad, \
                             const double wetDepth );

        /* Getters: */

        /* Aircraft name */
        inline std::string Name() const { return Name_; }
        /* Flight velocity */
        inline double VFlight() const { return vFlight_ms_; }
        /* Mach number */
        inline double Mach() const { return machNumber_; }
        /* Wingspan */
        inline double Wingspan() const { return wingspan_; }
        /* Current mass */
        inline double currMass() const { return currMass_; }
        /* Fuel flow */
        inline double FuelFlow() const { return engine_.getFuelFlow() * engNumber_; }
        /* Engine number */
        inline UInt EngNumber() const { return engNumber_; }
        /* Mean vertical displacement */
        inline double deltaz1() const { return vortex_.delta_z1(); }
        /* Maximum vertical displacement */
        inline double deltazw() const { return vortex_.delta_zw(); }
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

        /* Aircraft name */
        std::string Name_;

        /* Flight speed & mach Number */
        double vFlight_ms_;
        double machNumber_;

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
