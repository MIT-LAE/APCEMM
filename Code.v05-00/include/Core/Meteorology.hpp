/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Meteorology Header File                                          */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 11/2/2018                                 */
/* File                 : Meteorology.hpp                           */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef METEOROLOGY_H_INCLUDED
#define METEOROLOGY_H_INCLUDED

#include "APCEMM.h"
#ifdef OMP
    #include "omp.h"
#endif /* OMP */

#include "Util/ForwardDecl.hpp"
#include "Core/Input_Mod.hpp"
#include "Util/MetFunction.hpp"
#include "Util/PhysFunction.hpp"
#include <netcdf>

using namespace netCDF;
using namespace netCDF::exceptions;

enum class MetVarLoadType : unsigned char {
    NoMetInput,
    MetInputDefault,
    TimeSeries
};

class Meteorology
{

    public:
        Meteorology( ) = default;
        Meteorology( const OptInput &optInput,      \
                     const double pressure_Pa,   \
                     const Vector_1D& yCoords,
                     const Vector_1D& yEdges);


        void applyUpdraft(double dt);
        void recalculateSimulationGrid();
        void regenerate(const Vector_1D& yCoord_new, const Vector_1D& yEdges_new, int nx_new );
        void Update( const double dt, const double solarTime_h, \
                     const double simTime_h, const double dTrav_x = 0, const double dTrav_y = 0);
        
        void updateTempPerturb();
        inline double alt( int j ) const { return altitude_[j]; }
	    inline double press( int j ) const { return pressure_[j]; }
        inline double shear( int j ) const { return shear_[j]; }
        inline double shear() const { return i_Zp_ == -1 ? shear_[0] : shear_[i_Zp_]; }

        inline double temp( int j, int i ) const { return tempTotal_[j][i]; }
        inline double airMolecDens( int j, int i ) const { return airMolecDens_[j][i]; } // molecules/cm3
        inline double H2O( int j, int i ) const { return H2O_[j][i]; }

        //For getting the temp, rhw, and satdepth corresponding to initial pressure when using met input
        //TODO: Fix these functions and delete the _user variables, just calculate it from the reference altitude.

        inline double tempAtAlt(double alt) const { 
            return met::linInterpMetData(altitude_, tempBase_, alt);
        }
        inline double tempRef() const { return tempAtAlt(altitudeRef_); }
        inline double rhwAtAlt(double alt) const {
            Vector_1D rhwVec(yCoords_.size());
            Vector_1D h2o1D = H2O_1D();
            for(std::size_t j = 0; j < yCoords_.size(); j++) {
                rhwVec[j] = physFunc::H2OToRHw(h2o1D[j], tempBase_[j]);
            }
            return met::linInterpMetData(altitude_, rhwVec, alt);
        }
        inline double rhwRef() const { return rhwAtAlt(altitudeRef_); }
        inline double rhiAtAlt(double alt) const { return physFunc::RHwToRHi(rhwAtAlt(altitudeRef_), tempAtAlt(altitudeRef_)); }
        inline double rhiRef() const { return rhiAtAlt(altitudeRef_); }
        inline double shearAtAlt(double alt) const { 
            return met::linInterpMetData(altitude_, shear_, alt);
        }
        inline double shearRef() const { return shearAtAlt(altitudeRef_); }
        inline double lastOmega() const { return oldPressureVelocity_; }
        inline double satdepthUser() const { return satdepth_estimate_; }
        inline double referenceAlt() const { return altitudeRef_; } //Altitude at y = 0
        inline double referencePress() const { return pressureRef_; } //Pressure at y = 0

        inline const Vector_1D& tempBase() const { return tempBase_; }
        inline const Vector_1D H2O_1D() const {
            Vector_1D h2o1D(yCoords_.size());
            for(std::size_t j = 0; j < yCoords_.size(); j++) {
                h2o1D[j] = H2O_[j][0];
            }
            return h2o1D;
        }
        inline const Vector_2D& Temp() const { return tempTotal_; }
        inline const Vector_1D& Press() const { return pressure_; }
        inline const Vector_1D& Shear() const { return shear_; }
        inline const Vector_2D& H2O_field() const { return H2O_; }
        inline const Vector_2D& AirND_field() const { return airMolecDens_; }
        inline const Vector_1D& VertVeloc() const { return vertVeloc_; }
        inline const Vector_1D& AltEdges() const { return altitudeEdges_; }
        inline const Vector_1D& PressEdges() const { return pressureEdges_; }
        inline const Vector_1D& Altitude() const { return altitude_; }
        inline const Vector_1D& yCoords() const { return yCoords_; }
        inline const Vector_1D& yEdges() const { return yEdges_; }

        inline Vector_1D dy_vec() const {
            Vector_1D dy(ny_);
            for(int j = 0; j < ny_; j++) {
                dy[j] = yEdges_[j+1] - yEdges_[j];
            }  
            return dy;
        }
        
    private:
        inline void zeroVectors() { 
            tempTotal_ = Vector_2D(ny_, Vector_1D (nx_, 0));
            tempPerturbation_ = Vector_2D(ny_, Vector_1D (nx_, 0));
            airMolecDens_ = Vector_2D(ny_, Vector_1D (nx_, 0));
            H2O_ = Vector_2D(ny_, Vector_1D (nx_, 0));

            tempBase_ = Vector_1D(ny_, 0);
            shear_ = Vector_1D(ny_, 0);
            rhi_ = Vector_1D(ny_, 0);
            vertVeloc_ = Vector_1D(ny_, 0);
            altitude_ = Vector_1D(ny_, 0);
            pressure_ = Vector_1D(ny_, 0);
            altitudeEdges_ = Vector_1D(ny_ + 1, 0);
            pressureEdges_ = Vector_1D(ny_ + 1, 0);
        }
        inline void initMetLoadTypes(const OptInput& optInput) {
            tempLoadType_ = optInput.MET_TEMPTIMESERIES ? MetVarLoadType::TimeSeries : MetVarLoadType::MetInputDefault;
            rhLoadType_ = optInput.MET_RHTIMESERIES ? MetVarLoadType::TimeSeries : MetVarLoadType::MetInputDefault;
            shearLoadType_ = optInput.MET_SHEARTIMESERIES ? MetVarLoadType::TimeSeries : MetVarLoadType::MetInputDefault;
            vertVelocLoadType_ = optInput.MET_VERTVELOCTIMESERIES ? MetVarLoadType::TimeSeries : MetVarLoadType::MetInputDefault;
        }

        void readMetVar( const NcFile& dataFile, std::string varName, Vector_2D& vec_ts, bool timeseries);

        // APCEMM v1.2
        void estimateMetDataAltitudes();
        void readMetDataFromFile( const NcFile& dataFile );
        void interpolateMetToSimulationGrid();
        void updateMetData(double simTime_h);
        void estimateSimulationGridPressures();
        void updateSimulationGridProperties();
        double estimateSatDepth();

        Vector_1D interpMetTimeseriesData(double simTime_h, const Vector_2D& ts_data, bool timeseries) const;

        void updateAirMolecDens();
        //For passing params to EPM (initial params at specified altitude)
        int i_Zp_ = -1; //Index for initial variables in altitude vector/coordinates if using met input
        double temp_user_;
        double shear_user_;
        double rhw_user_;
        double satdepth_estimate_;

        //Grid data
        int nx_;
        int ny_;
        Vector_1D yCoords_;
        Vector_1D yEdges_;


        /* For processing met input */
        double met_dt_h_;
        int altitudeDim_;
        int timeDim_;
        Vector_1D pressureInit_; 
        Vector_1D altitudeInit_;
        Vector_1D tempInit_;
        Vector_1D shearInit_;
        Vector_1D rhiInit_;
        Vector_1D vertVelocInit_;
        Vector_1D lapseInit_; 
        Vector_1D pressureEdgesInit_; 
        Vector_1D altitudeEdgesInit_; 

	    Vector_2D tempTimeseriesData_;
        Vector_2D shearTimeseriesData_;
        Vector_2D rhiTimeseriesData_;
        Vector_2D vertVelocTimeseriesData_;

        /* Ambient input parameters */
        double altitudeRef_;
        double pressureRef_;

        /* Temperature variations */
        double diurnalAmplitude_; // [K]
        double diurnalPhase_; // [hours]
        double diurnalPert_; // [K]
        double turbTempPertAmplitude_; // [K]

        /* Assume that pressure only depends on the vertical coordinate */

        /* Temperature, air density and humidity fields can potentially be
         * 2D fields */
        Vector_2D tempTotal_;
        Vector_2D tempPerturbation_;
        Vector_1D tempBase_; //Temp without the perturbations
        Vector_2D airMolecDens_;
        Vector_2D H2O_;
        Vector_1D shear_;
        Vector_1D rhi_;
        Vector_1D vertVeloc_; // [m/s]

        Vector_1D altitude_;
	    Vector_1D pressure_;
        Vector_1D altitudeEdges_;
        Vector_1D pressureEdges_;


        //Information on loading
        bool useMetFileInput_;
        MetVarLoadType tempLoadType_;
        MetVarLoadType rhLoadType_;
        MetVarLoadType shearLoadType_;
        MetVarLoadType vertVelocLoadType_;
        bool interpTemp_;
        bool interpRH_;
        bool interpShear_;
        bool interpVertVeloc_;

        // For handling vertical advection due to updrafts
        double oldPressureVelocity_;

};

double hydrostaticDeltaAltitude(double lapseRate, double refTemperature, double refPressure, double targPressure);

#endif /* METEOROLOGY_H_INCLUDED */
