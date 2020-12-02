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
#include <iostream>
#include <iomanip>
#include <cstring>

#include "Core/Engine.hpp"
#include "Core/Vortex.hpp"
#include "Util/PhysConstant.hpp"

class Aircraft 
{
    public:

        /* Constructors */

        Aircraft( );
        Aircraft( const char *aircraftName, RealDouble aircraftMass, \
                  RealDouble temperature_K, RealDouble pressure_Pa,  \
                  RealDouble relHumidity_w );

        /* Destructor */

        ~Aircraft( );

        /* Copy */

        Aircraft( const Aircraft &ac );

        /* Copy */

        Aircraft& operator=( const Aircraft &ac ); 

        /* Debug */
        void Debug( ) const;

        /* Compute vortex losses */

        RealDouble VortexLosses( const RealDouble EI_Soot, const RealDouble EI_SootRad, \
                             const RealDouble wetDepth );

        /* Getters: */

        /* Aircraft name */
        std::string Name() const { return Name_; }
        /* Flight velocity */
        RealDouble VFlight() const { return vFlight_ms_; }
        /* Mach number */
        RealDouble Mach() const { return machNumber_; }
        /* Wingspan */
        RealDouble Wingspan() const { return wingspan_; }
        /* Max take-off weight */
        RealDouble MTOW() const { return MTOW_; }
        /* Current mass */
        RealDouble currMass() const { return currMass_; }
        /* Fuel flow */
        RealDouble FuelFlow() const { return engine_.getFuelFlow() * engNumber_; }
        /* Engine number */
        UInt EngNumber() const { return engNumber_; }
        /* Mean vertical displacement */
        RealDouble deltaz1() const { return vortex_.delta_z1(); }
        /* Maximum vertical displacement */
        RealDouble deltazw() const { return vortex_.delta_zw(); }
        /* Engine */
        Engine engine() const { return engine_; }

        /* Setters for engine properties */
        void setEI_NOx(const RealDouble NOx);
        void setEI_CO(const RealDouble CO);
        void setEI_HC(const RealDouble HC);
        void setEI_Soot(const RealDouble Soot);
        void setSootRad(const RealDouble sootRad);
        void setFuelFlow(const RealDouble ff);
        void setVFlight(const RealDouble Vf, RealDouble temperature_K);
        void setEngNumber(const RealDouble nEng);
        void setWingspan(const RealDouble span);

        /* Engine */
        Engine engine_;

    protected:

        /* Aircraft name */
        std::string Name_;

        /* Flight speed & mach Number */
        RealDouble vFlight_ms_;
        RealDouble machNumber_;

        /* Dimensions */
        RealDouble wingspan_; /* [m] */

        /* Weight */
        RealDouble MTOW_; /* [kg] */
        RealDouble currMass_; /* [kg] */

        /* Number of engines */
        UInt engNumber_;

        /* Vortex */
        Vortex vortex_;

    private:

};

#endif /* AIRCRAFT_H_INCLUDED */
