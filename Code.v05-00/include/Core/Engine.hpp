/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Engine Header File                                               */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : Engine.hpp                                */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef ENGINE_H_INCLUDED
#define ENGINE_H_INCLUDED

#include <string>
#include <fstream>
#include <cstring>
#include <vector>
#include <cmath>


class Engine
{

    /* Forward declaration */
    class Aircraft;

    friend class Aircraft;

    public:
        
        Engine( );
        Engine( const char *engineName, std::string engineFileName, double tempe_K, double pres_Pa, double relHum_w, double machNumber );
        Engine( const Engine &e );
        Engine& operator=( const Engine &e );
        ~Engine( );
        void OpenFile( const char *fileName, std::ifstream &file );
        void CloseFile( std::ifstream &file );
        bool GetEDB( std::ifstream &file, const char *engineName, std::string &idle, std::string &approach, std::string &climbout, std::string &takeoff );
        std::string getName() const;
        double getEI_NOx() const;
        double getEI_NO() const;
        double getEI_NO2() const;
        double getEI_HNO2() const;
        double getEI_CO() const;
        double getEI_HC() const;
        double getEI_Soot() const;
        double getSootRad() const;
        double getFuelFlow() const;
        double getTheta() const;
        double getDelta() const;
        void setEI_NOx(const double NOx);
        void setEI_CO(const double CO);
        void setEI_HC(const double HC);
        void setEI_Soot(const double Soot);
        void setSootRad(const double sootRad);
        void setFuelFlow(const double ff);

    protected:

        /* Engine name */
        std::string Name;

        /* Emission indices */
        /** NOx emission indices */
        double EI_NOx; /* [g/kg fuel] */
        double EI_NO;
        double EI_NO2;
        double EI_HNO2;
        /** CO emission index */
        double EI_CO;  /* [g/kg fuel] */
        /** UHC emission index */
        double EI_HC;  /* [g/kg fuel] */

        /* Soot */
        double EI_Soot; /* [g/kg fuel] */
        double SootRad; /* [m] */

        /* Fuel flow */
        double fuelflow; /* [kg fuel/s] */

        /* Atmospheric characteristics */
        /** Ratio of ambient to sea-level conditions */
        double theta; /* [-] */
        double delta; /* [-] */

    private:

        static const char * const engineFileName;

        std::vector<double> ratedThrust = std::vector<double>(4);
        std::vector<double> LTO_fuelflow = std::vector<double>(4);
        std::vector<double> LTO_NOx = std::vector<double>(4);
        std::vector<double> LTO_CO = std::vector<double>(4);
        std::vector<double> LTO_HC = std::vector<double>(4);

        double NOxtoHNO2 = 0.015;
        double NOxtoNO2 = 0.20 * ( 1.0 - NOxtoHNO2 );
        double NOxtoNO = 1.0 - NOxtoNO2 - NOxtoHNO2;

};

#endif /* ENGINE_H_INCLUDED */
