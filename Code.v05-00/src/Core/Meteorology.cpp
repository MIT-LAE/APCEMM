/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Meteorology Program File                                         */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 11/2/2018                                 */
/* File                 : Meteorology.cpp                           */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <algorithm>
#include "Util/PhysFunction.hpp"
#include "Util/PhysConstant.hpp"
#include "Core/Parameters.hpp"
#include "Core/Meteorology.hpp"
#include "Util/MC_Rand.hpp"

Meteorology::Meteorology( const OptInput &optInput,
                          const AmbientMetParams& ambParams,
                          const Vector_1D& yCoords,
                          const Vector_1D& yEdges):
    rhi_far_(optInput.MET_SUBSAT_RHI),
    nx_(optInput.ADV_GRID_NX),
    ny_(optInput.ADV_GRID_NY),
    yCoords_(yCoords),
    yEdges_(yEdges),
    met_dt_h_(optInput.MET_DT),
    ambParams_(ambParams),
    pressureRef_(ambParams.press_Pa),
    turbTempPertAmplitude_(optInput.MET_TEMP_PERTURB_AMPLITUDE),
    useMetFileInput_(optInput.MET_LOADMET),
    interpTemp_(optInput.MET_INTERPTEMPDATA),
    interpRH_(optInput.MET_INTERPRHDATA),
    interpShear_(optInput.MET_INTERPSHEARDATA),
    interpVertVeloc_(optInput.MET_INTERPVERTVELOC)
{

    zeroVectors();

    initMetLoadTypes(optInput);

    //Set reference altitude from reference pressure
    met::ISA_pAlt(altitudeRef_, pressureRef_);
    //Steps:
    // 1) Calculate lapse rate based on params
    // 2) Initalize the initial and time-dependent temperature, h2o, and shear

    //Set lapse rate (either calculated using saturation depth, or specified.)
    lapseRate_ = optInput.MET_FIXDEPTH ? met::ComputeLapseRate( ambParams_.temp_K, ambParams_.rhi, optInput.MET_DEPTH )
                                       : optInput.MET_LAPSERATE;
                                    
    //Set met depth for cases where we don't load rh from met input
    met_depth_ =  optInput.MET_FIXDEPTH ? optInput.MET_DEPTH : 200.0;

    /* The temperature can vary depending on the time of day. (Diurnal variations)
     * The data on the diurnal temperature amplitude and phase is gathered
     * from: 
     * Seidel, D. J., M. Free, and J. Wang (2005), Diurnal cycle of 
     * upper-air temperature estimated from radiosondes,J. Geophys.Res.,
     * 110, D09102, doi:10.1029/2004JD005526.*/


    //TODO: Fix this diurnal phase hardcode.
    if ( optInput.MET_DIURNAL ) {
        /* Based on values around 220 hPa */
        diurnalAmplitude_  = 0.2 ; /* [K] */
        diurnalPhase_ = 12.0; /* [hrs] */
    } else {
        diurnalAmplitude_  = 0.0 ; /* [K] */
        diurnalPhase_ = 12.0; /* [hrs] */
    }

    diurnalPert_ = diurnalAmplitude_ * cos( 2.0E+00 * physConst::PI * ( ambParams_.solarTime_h - diurnalPhase_ ) / 24.0E+00 );

    //Declare met input file here, keeps ownership of file to constructor only.
    //Dont want to deal with using a ptr and worry about what functions deference it.
    NcFile dataFile;
    if( optInput.MET_LOADMET ) {
        dataFile.open( optInput.MET_FILENAME.c_str(), NcFile::read );
    }

    try {
        initAltitudeAndPress( dataFile );
    }
    catch (NcException& e) {
        throw std::runtime_error("Could not parse altitude and pressure data from specified met input file");
    }

    try {
        initTemperature( dataFile );
    }
    catch (NcException& e) {
        throw std::runtime_error("Could not parse temperature data from specified met input file");
    }

    try {
        initH2O( dataFile, optInput );
    }
    catch (NcException& e) {
        throw std::runtime_error("Could not parse relative humidity data from specified met input file");
    }

    try {
        initShear( dataFile );
    }
    catch (NcException& e) {
        throw std::runtime_error("Could not parse shear data from specified met input file");
    }

    try {
        initVertVeloc( dataFile );
    }
    catch (NcException& e) {
        throw std::runtime_error("Could not parse vertical velocity data from specified met input file");
    }


    double invkB = 1.00E-06 / physConst::kB;

    #pragma omp parallel for if ( !PARALLEL_CASES ) 
    for (int jNy = 0; jNy < ny_; jNy++ ) {
        airMolecDens_[jNy].assign(nx_, pressure_[jNy] / tempBase_[jNy] * invkB);
    }

} /* End of Meteorology::Meteorology */

void Meteorology::regenerate( const Vector_1D& yCoord_new, const Vector_1D& yEdges_new, int nx_new ) {
    double dy_new = yEdges_new[1] - yEdges_new[0];
    int ny_new = yCoord_new.size();
    double alt_y0 = altitudeEdges_[0] + (yEdges_new[0] - yEdges_[0]);
    Vector_1D alt_new(ny_new);
    Vector_1D altEdges_new(ny_new + 1);
    Vector_1D press_new(ny_new);
    Vector_1D pressEdges_new(ny_new + 1);

    std::generate(altEdges_new.begin(), altEdges_new.end(), [&, j = 0] () mutable { return alt_y0 + dy_new * j++; });
    std::generate(alt_new.begin(), alt_new.end(), [&, j = 0] () mutable { return alt_y0 + dy_new * (0.5 + j++); });

    met::ISA(alt_new, press_new);
    met::ISA(altEdges_new, pressEdges_new);

    tempBase_.resize(ny_new);
    tempTotal_ = Vector_2D(ny_new, Vector_1D(nx_new));
    tempPerturbation_ = Vector_2D(ny_new, Vector_1D(nx_new));
    H2O_ = Vector_2D(ny_new, Vector_1D(nx_new));
    airMolecDens_ = Vector_2D(ny_new, Vector_1D(nx_new));
    shear_.resize(ny_new);
    yCoords_ = yCoord_new;
    yEdges_ = yEdges_new;
    ny_ = ny_new;
    nx_ = nx_new;
    altitude_ = alt_new;
    pressure_ = press_new;
    altitudeEdges_ = altEdges_new;
    pressureEdges_ = pressEdges_new;
    if(vertVelocLoadType_ != MetVarLoadType::NoMetInput) {
        vertVeloc_.resize(ny_new);
    }

    //Regenerate temp field
    if(tempLoadType_ == MetVarLoadType::NoMetInput) {
        initTempNoMet(yCoord_new);
    }
    else {
        for(int j = 0; j < ny_new; j++) {
            int i_Z = met::nearestNeighbor( altitudeInit_, alt_new[j]);
            double tempInterp = met::linInterpMetData(altitudeInit_, tempInit_, alt_new[j]);
            tempBase_[j] = interpTemp_ ? tempInterp : tempInit_[i_Z];
            tempTotal_[j].assign(nx_new, tempBase_[j]);
        }
    }
    
    //Regenerate H2O
    if(rhLoadType_ == MetVarLoadType::NoMetInput) {
        initH2ONoMet(yCoord_new);
    }
    else {
        for(int j = 0; j < ny_new; j++) {
            int i_Z = met::nearestNeighbor( altitudeInit_, alt_new[j]);
            double rhiInterp = met::linInterpMetData(altitudeInit_, rhiInit_, alt_new[j]);
            double rhiToUse = interpRH_ ? rhiInterp : rhiInit_[i_Z];
            H2O_[j].assign(nx_new, physFunc::RHiToH2O(rhiToUse, tempBase_[j]));
        }
    }

    // Regenerate shear
    for ( int j = 0;  j < ny_new; j++ ) {
        int i_Z = met::nearestNeighbor( altitudeInit_, alt_new[j]);
        if(shearLoadType_ == MetVarLoadType::NoMetInput) {
            shear_[j] = ambParams_.shear;
        }
        else {
            double shear_local = met::linInterpMetData(altitudeInit_, shearInit_, alt_new[j]);
            shear_[j] = interpShear_ ? shear_local : shearInit_[i_Z];
        }
    }

    // Regenerate vert veloc
    if(vertVelocLoadType_ != MetVarLoadType::NoMetInput) {
        for ( int j = 0;  j < ny_new; j++ ) {
            int i_Z = met::nearestNeighbor( altitudeInit_, alt_new[j]);
            double w_local = met::linInterpMetData(altitudeInit_, vertVelocInit_, alt_new[j]);
            vertVeloc_[j] = interpVertVeloc_ ? w_local : vertVelocInit_[i_Z];
        }
    }
    updateAirMolecDens();
}

void Meteorology::Update( const double dt, const double solarTime_h, \
                          const double simTime_h, const double dTrav_x, const double dTrav_y)
{
    diurnalPert_ = diurnalAmplitude_ * cos( 2.0E+00 * physConst::PI * ( solarTime_h - diurnalPhase_ ) / 24.0E+00 );

    //Update altitude if dTrav_y is nonzero 
    //dTrav_x, dTrav_y are for velocities not accounted for in the meteorological vertical velocity.
    if( dTrav_y != 0 ) { 
        for (int jNy = 0; jNy < ny_; jNy++ )
            altitude_[jNy] +=  dTrav_y;
        met::ISA( altitude_, pressure_ );
    }

    //First, we take the vertical velocity at the simtime specifed outside, typically halfway into the timestep.
    //Then advect to the new altitude based on the pressure velocity at the reference altitude at time simTime.
    updateVertVeloc(simTime_h);
    vertAdvectAltPress(dt); 

    /* User defined fields can be set here ! */
    updateTemperature(solarTime_h, simTime_h);
    updateShear(simTime_h);
    updateH2O(simTime_h);
    updateAirMolecDens();
} /* End of Meteorology::UpdateMet */

void Meteorology::initAltitudeAndPress( const NcFile& dataFile ) {
    //Must call this before the other initialize functions!
    if( !useMetFileInput_ ) {
            
        for (int i = 0; i < ny_; i++ ) {
            altitude_[i] = altitudeRef_ + yCoords_[i];
        }

        for (int i = 0; i < ny_ + 1; i++ ) {
            altitudeEdges_[i] = altitudeRef_ + yEdges_[i];
        }
        
        met::ISA( altitude_, pressure_ );
        met::ISA( altitudeEdges_, pressureEdges_ );
        return;
    }

    /*
        Met data format:
        2 dimensions: altitude and time
        RHw, Temp, Shear given as time series.
    */

    /* Identify the length of variables in input file */
    altitudeDim_ = dataFile.getDim("altitude").getSize();
    timeDim_ = dataFile.getDim("time").getSize();

    /* Extract pressure and altitude from input file. */
    altitudeInit_.resize(altitudeDim_);
    pressureInit_.resize(altitudeDim_);

    //The netcdf API is garbage and makes you pass in a C style array despite being branded as a "C++" library
    NcVar altitude_ncVar = dataFile.getVar("altitude");
    NcVar pressure_ncVar = dataFile.getVar("pressure");
    altitude_ncVar.getVar(altitudeInit_.data());
    pressure_ncVar.getVar(pressureInit_.data());

    //Cache initial pressure and altitude values for generating later met vars.
    for (int i = 0; i < altitudeDim_; i++ ) {
        pressureInit_[i] *= 100.0; //convert from hPa to Pa
        altitudeInit_[i] *= 1000.0; //convert from km to m
    }

    for ( int j = 0; j < ny_; j++ ) {
        altitude_[j] = altitudeRef_ + yCoords_[j];
    }

    for ( int j = 0; j < ny_ + 1; j++ ) {
        altitudeEdges_[j] = altitudeRef_ + yEdges_[j];
    }
    met::ISA(altitude_, pressure_);
    met::ISA(altitudeEdges_, pressureEdges_);
    i_Zp_ = met::nearestNeighbor( pressure_, pressureRef_); 
}

void Meteorology::readMetVar( const NcFile& dataFile, std::string varName, Vector_2D& vec_ts, bool timeseries) {
    NcVar ncvar = dataFile.getVar(varName.c_str());
    bool supportsTimeseries = ncvar.getDimCount() == 2;
    if( !supportsTimeseries && timeseries ) {
        throw std::runtime_error("Variable\"" + varName + "\" in met input file does not support time series input! Please set the corresponding time series input option to false.");
    }
    int dataLength = supportsTimeseries ? altitudeDim_ * timeDim_ : altitudeDim_;
    double data[dataLength]; //flattened array to hold values for all altitude and time
    ncvar.getVar(data);

    vec_ts = Vector_2D(altitudeDim_, Vector_1D(timeDim_, 0));
    

    for ( int i = 0; i < altitudeDim_; i++ ) {
        for ( int itime = 0; itime < timeDim_; itime++ ) {
            int data_idx;

            if(!supportsTimeseries) data_idx = i;
            else data_idx = timeseries ? i*timeDim_ + itime : i*timeDim_; 

            vec_ts[i][itime] = data[data_idx];
        }
    }
}

void Meteorology::initTempNoMet (const Vector_1D& yCoords) {
    #pragma omp parallel for default(shared)
    for ( std::size_t j = 0; j < yCoords.size(); j++ ) {
        //Convention: lower altitude than reference= negative y, higher = positive y
        double temp_local = ambParams_.temp_K + yCoords[j] * lapseRate_ + diurnalPert_;
        tempBase_[j] = temp_local;
        tempTotal_[j].assign(nx_, tempBase_[j]);
    }
}
void Meteorology::initTemperature( const NcFile& dataFile ) {

    if ( tempLoadType_ == MetVarLoadType::NoMetInput ) {
        initTempNoMet(yCoords_);
        return;
    }

    readMetVar(dataFile, "temperature", tempTimeseriesData_, tempLoadType_ == MetVarLoadType::TimeSeries); 
    
    tempInit_.resize(altitudeDim_);
    for (int i = 0; i < altitudeDim_; i++) {
        tempInit_[i] = tempTimeseriesData_[i][0];
    }

    /* Identify closest temperature to given pressure */
    /* Loop round each vertical layer to estimate temperature */

    #pragma omp parallel for if(!PARALLEL_CASES)
    for ( int j = 0; j < ny_; j++ ) {

        /* Find the closest values above and below the central pressure */
        int i_Z = met::nearestNeighbor( altitudeInit_, altitude_[j]);
        double tempInterp = met::linInterpMetData(altitudeInit_, tempInit_, altitude_[j]);
        tempBase_[j] = interpTemp_ ? tempInterp : tempInit_[i_Z];
        tempTotal_[j].assign(nx_, tempBase_[j]);
    }
}

void Meteorology::initH2ONoMet (const Vector_1D& yCoords) {

    //Unfortunately, need to set some arbitrary values for thickness of supersaturated layer above, as well as the farfield RH. 
    //We set the moist layer depth to be the same on the top and bottom, essentially assuming the contrail spawns in the middle of the layer.
    //Statistically, ^ might be inaccurate and could have room for calibration.
    /* RHi layer centered on y = 0 */
    double moist_layer_top_y = met_depth_;
    double moist_layer_bot_y = -met_depth_;
    satdepth_user_ = met_depth_;

    double RH_star = ambParams_.rhi;
    double RH_far = rhi_far_;
    #pragma omp parallel for default(shared)
    for ( std::size_t j = 0; j < yCoords.size(); j++ ) {
        double H2O_local = yCoords[j] > moist_layer_bot_y && yCoords[j] < moist_layer_top_y
                            ? physFunc::RHiToH2O(RH_star, tempBase_[j])
                            : physFunc::RHiToH2O(RH_far, tempBase_[j]);
        H2O_[j].assign(nx_, H2O_local);
    }
}

void Meteorology::initH2O( const NcFile& dataFile, const OptInput& optInput ) { 
    //Cannot call this before initTemperature!

    if( rhLoadType_ == MetVarLoadType::NoMetInput ) {
        initH2ONoMet(yCoords_);
        return;
    }

    readMetVar(dataFile, "relative_humidity_ice", rhiTimeseriesData_, rhLoadType_ == MetVarLoadType::TimeSeries); 
    
    rhiInit_.resize(altitudeDim_);
    
    for (int i = 0; i < altitudeDim_; i++) {
        //Scale RHi if specified
        if (optInput.MET_HUMIDSCAL_MODIFICATION_SCHEME == "scaling") {
            rhiInit_[i] = met::rhiCorrection(rhiTimeseriesData_[i][0], optInput.MET_HUMIDSCAL_SCALING_A, optInput.MET_HUMIDSCAL_SCALING_B);
        }
        else if (optInput.MET_HUMIDSCAL_MODIFICATION_SCHEME == "constant") {
            rhiInit_[i] = optInput.MET_HUMIDSCAL_CONST_RHI;
        }
        else {
            rhiInit_[i] = rhiTimeseriesData_[i][0];
        }
    }
    Vector_1D localRHi(ny_);
    /* Identify closest RH to given pressure at level i_Z */
    /* Loop round each vertical layer to estimate RH */
    #pragma omp parallel for\
    if      ( !PARALLEL_CASES ) \
    default ( shared          )
    for ( int jNy = 0; jNy < ny_; jNy++ ) {

        /* Find the closest values above and below the central pressure */
        int i_Z = met::nearestNeighbor( altitudeInit_, altitude_[jNy]);

        double rhiInterp = met::linInterpMetData(altitudeInit_, rhiInit_, altitude_[jNy]);
        double rhiToUse = interpRH_ ? rhiInterp : rhiInit_[i_Z];
        localRHi[jNy] = rhiToUse;
        H2O_[jNy].assign(nx_, physFunc::RHiToH2O(rhiToUse, tempBase_[jNy]));

    }

    int i_Z = met::nearestNeighbor( pressure_, ambParams_.press_Pa);
    double dy = yCoords_[1] - yCoords_[0];

    //Throws exception if domain not big enough.
    //TODO: Decide on what to do about variable saturation depths (i.e. advecting contrail through space).
    Vector_1D localRHw(ny_);
    for (int j = 0; j < ny_; j++) {
        localRHw[j] = physFunc::RHiToRHw(localRHi[j], tempBase_[j]);
    }
    //Yes, satdepth_calc converts this RHw right back into RHi to calculate it. Whatever it doesn't affect performance. -MX
    try {
        satdepth_user_ = met::satdepth_calc(localRHw, tempBase_, altitude_, i_Z, std::abs(yCoords_[0]) + dy/2);
    }
    catch (std::out_of_range &e){
        if(optInput.MET_HUMIDSCAL_MODIFICATION_SCHEME != "constant")
            std::cout << "WARNING: end of initial domain does not cover saturation depth." << std::endl;
        satdepth_user_ = 1000.0;
    }
}

void Meteorology::initShear ( const NcFile& dataFile ) {

    if ( shearLoadType_ == MetVarLoadType::NoMetInput ) {
        shear_.assign(ny_, ambParams_.shear);
        return;
    }

    readMetVar(dataFile, "shear", shearTimeseriesData_, shearLoadType_ == MetVarLoadType::TimeSeries); 
    shearInit_.resize(altitudeDim_);
    for (int i = 0; i < altitudeDim_; i++) {
        shearInit_[i] = shearTimeseriesData_[i][0];
    }

    for ( int jNy = 0;  jNy < ny_; jNy++ ) {

        int i_Z = met::nearestNeighbor( altitudeInit_, altitude_[jNy]);
        double shear_local = met::linInterpMetData(altitudeInit_, shearInit_, altitude_[jNy]);
        shear_[jNy] = interpShear_ ? shear_local : shearInit_[i_Z];

    }
}

void Meteorology::initVertVeloc ( const NcFile& dataFile ) {
    if ( vertVelocLoadType_ == MetVarLoadType::NoMetInput ) {
        vertVeloc_.assign(ny_, 0);
        return;
    }

    //Vert veloc is assumed default as timeseries input.
    readMetVar(dataFile, "w", vertVelocTimeseriesData_, vertVelocLoadType_ == MetVarLoadType::TimeSeries);
    
    vertVelocInit_.resize(altitudeDim_);
    for (int i = 0; i < altitudeDim_; i++) {
        vertVelocInit_[i] = vertVelocTimeseriesData_[i][0];
    }

    for ( int jNy = 0;  jNy < ny_; jNy++ ) {

        int i_Z = met::nearestNeighbor( altitudeInit_, altitude_[jNy]);
        double w_local = met::linInterpMetData(altitudeInit_, vertVelocInit_, altitude_[jNy]);
        vertVeloc_[jNy] = interpVertVeloc_ ? w_local : vertVelocInit_[i_Z];
    }
}

Vector_1D Meteorology::interpMetTimeseriesData(double simTime_h, const Vector_2D& ts_data, bool timeseries) const {

    int itime = timeseries 
              ? std::min(static_cast<int>(simTime_h / met_dt_h_), timeDim_ - 1)
              : 0;

    double before;
    double after;
    Vector_1D interp(altitudeDim_);

    /* Extract temperature data before and after current time, and interpolate */
    for ( int i = 0; i < altitudeDim_; i++ ) {
        before = ts_data[i][itime];
        if ( itime >= timeDim_ - 1 ) {
            std::cout <<  "WARNING: Simulation time exceeded extent of timeseries data. Using last entry provided in timeseries data." << std::endl;
            after = ts_data[i][itime];
        }
        else {
            after = ts_data[i][itime+1];
        }
        interp[i] = before + ( after - before ) * ( simTime_h - itime * met_dt_h_ );
    }
    return interp;
}

void Meteorology::updateTemperature(double solarTime_h, double simTime_h) {

    double diurnalPertPrev = diurnalPert_;
    diurnalPert_ = diurnalAmplitude_ * cos( 2.0 * physConst::PI * (solarTime_h - diurnalPhase_) / 24.0);
    double deltaDiurnalPert = diurnalPert_ - diurnalPertPrev;
    
    if( tempLoadType_ == MetVarLoadType::NoMetInput ) {
        #pragma omp parallel for if(!PARALLEL_CASES)
        for (int j = 0; j < ny_; j++ ) {
            tempBase_[j] += deltaDiurnalPert;
            /* Unit check: Y [m] * LapseRate [K/m] */
            for (int i = 0; i < nx_; i++ ) {
                tempTotal_[j][i] = tempBase_[j] + tempPerturbation_[j][i];
            }
        }
        return;
    }
    bool timeseries = (tempLoadType_ == MetVarLoadType::TimeSeries);
    tempInit_ = interpMetTimeseriesData(simTime_h, tempTimeseriesData_, timeseries);

    #pragma omp parallel for if (!PARALLEL_CASES)
    for ( int j = 0; j < ny_; j++ ) {
        int i_Z = met::nearestNeighbor( altitudeInit_, altitude_[j] );
        double temp_local = met::linInterpMetData(altitudeInit_, tempInit_, altitude_[j]);
        tempBase_[j] = interpTemp_ ? temp_local : tempInit_[i_Z];
        for ( int i = 0; i < nx_; i++ ) {
            tempTotal_[j][i] =  tempBase_[j] + tempPerturbation_[j][i];
        }

    }

}

void Meteorology::updateH2O(double simTime_h) { 
    /*
        Updating the RH from the data at all inherently breaks the idea of performing ice growth
        on a given H2O field and then reusing the H2O field on the next timestep. Therefore,
        if RH timeseries is not specified, the RH field will not be changed by the update function.
     */
    if (rhLoadType_ == MetVarLoadType::NoMetInput) return;
    rhiInit_ = interpMetTimeseriesData(simTime_h, rhiTimeseriesData_, true);

    #pragma omp parallel for if (!PARALLEL_CASES)
    for ( int j = 0; j < ny_; j++ ) {
        double rh_local = met::linInterpMetData(altitudeInit_, rhiInit_, altitude_[j]);
        double h2o_local = physFunc::RHiToH2O(rh_local, tempBase_[j]);
        H2O_[j].assign(nx_, h2o_local);
    }
}


void Meteorology::updateShear(double simTime_h) {
    if( shearLoadType_ == MetVarLoadType::NoMetInput )  return;

    bool timeseries = (shearLoadType_ == MetVarLoadType::TimeSeries);
    shearInit_ = interpMetTimeseriesData(simTime_h, shearTimeseriesData_, timeseries);

    for ( int jNy = 0; jNy < ny_; jNy++ ) {
        int i_Z = met::nearestNeighbor( altitudeInit_, altitude_[jNy] );
        double shear_local = met::linInterpMetData(altitudeInit_, shearInit_, altitude_[jNy]);
        shear_[jNy] = interpShear_ ? shear_local : shearInit_[i_Z];

    }
}

void Meteorology::updateVertVeloc(double simTime_h) {
    if( vertVelocLoadType_ == MetVarLoadType::NoMetInput )  return;

    bool timeseries = (vertVelocLoadType_ == MetVarLoadType::TimeSeries);
    vertVelocInit_ = interpMetTimeseriesData(simTime_h, vertVelocTimeseriesData_, timeseries);

    for ( int jNy = 0; jNy < ny_; jNy++ ) {
        int i_Z = met::nearestNeighbor( altitudeInit_, altitude_[jNy] );
        double w_local = met::linInterpMetData(altitudeInit_, vertVelocInit_, altitude_[jNy]);
        vertVeloc_[jNy] = interpVertVeloc_ ? w_local : vertVelocInit_[i_Z];

    }
}

void Meteorology::updateTempPerturb() {
    #pragma omp parallel for\
    if(!PARALLEL_CASES) \
    default(shared)
    for (int j = 0; j < ny_; j++){
        for(int i = 0; i < nx_; i++){
            double epsilon1 = fRand(-1.0, 1.0);
            double epsilon2 = fRand(-1.0, 1.0);
            tempPerturbation_[j][i] = epsilon1 * epsilon2 * turbTempPertAmplitude_;
            tempTotal_[j][i] = tempBase_[j] + tempPerturbation_[j][i]; //Bad practice of having 1 function update both the temp perturb and the total temp but whatever
        }
    }
}

void Meteorology::updateAirMolecDens() {
    double invkB = 1.00E-06 / physConst::kB;
    #pragma omp parallel for        \
    if      ( !PARALLEL_CASES ) \
    default ( shared          )
    for ( int jNy = 0; jNy < ny_; jNy++ ) {
        for ( int iNx = 0; iNx < nx_; iNx++ )
            airMolecDens_[jNy][iNx] = pressure_[jNy] / tempTotal_[jNy][iNx] * invkB;
    }
}

void Meteorology::vertAdvectAltPress(double dt) {
    /* For now, it only makes sense to vertically advect by a constant velocity in pressure.
       Otherwise, we would just violate conservation of mass, but because
       the continuity equation in pressure coordinates is not a function of density, this is fine.
       -MX */
    
    if(vertVelocLoadType_ == MetVarLoadType::NoMetInput) {
        return;
    }
    //Save old altitude edges and centers to update y coordinates
    Vector_1D altEdges_old = altitudeEdges_;
    Vector_1D altCenters_old = altitude_;

    //Convert w at current alt to pressure velocity (omega)
    //Use 4th order central finite difference to estimate derivative
    double p_jm2, p_jm1, p_jp1, p_jp2;
    double dp = 0.1;
    met::ISA(altitudeRef_ - 2*dp, p_jm2);
    met::ISA(altitudeRef_ - dp, p_jm1);
    met::ISA(altitudeRef_ + 2*dp, p_jp2);
    met::ISA(altitudeRef_ + dp, p_jp1);
    double dp_dz = (p_jm2 - 8*p_jm1 + 8*p_jp1 - p_jp2) / (12 * dp);

    double w_local = met::linInterpMetData(altitudeInit_, vertVelocInit_, altitudeRef_);    //dp/dt = dp/dz * dz/dt [Pa/s]

    double omega = dp_dz * w_local;

    // SDE 2024-08-12: Disabled temporarily while pressure velocity issues are fixed
    if (std::abs(omega) > 1.0e-10 ){
        std::cout << " -- WARNING: PRESSURE VELOCITIES CURRENTLY DISABLED -- " << std::endl;
    }
    return;

    /* Update edges instead of centers (altitude_, pressure_) because altitude-pressure relationship
       is nonlinear. Therefore, the relative position of the center will move hwen advecting at some
       constant pressure velocity. If we then assume a constant dx without moving the center, we will
       fail to conserve mass at the edges, and this effect will snowball in time. -MX */
    for(double& pE: pressureEdges_) {
        pE += omega * dt;
    }
    
    // double ytop_old = altitudeEdges_[altitudeEdges_.size() - 1];
    //Update altitude edges to go with updated pressure edges
    met::ISA_pAlt(altitudeEdges_, pressureEdges_);
    //Re-calculate geometric altitude centers
    for(int i = 0; i < ny_; i++) { 
        altitude_[i] = (altitudeEdges_[i] + altitudeEdges_[i+1]) / 2;
    }

    //Update pressure centers (based on altitude centers)
    met::ISA(altitude_, pressure_);

    //Update reference pressures/altitudes for next timestep
    pressureRef_ += omega * dt;
    met::ISA_pAlt(altitudeRef_, pressureRef_);

    //Update y coordinates based on these changes
    for (int i = 0; i < ny_; i++) {
        yCoords_[i] += (altitude_[i] - altCenters_old[i]);
    }
    for (int i = 0; i < ny_ + 1; i++) {
        yEdges_[i] += (altitudeEdges_[i] - altEdges_old[i]);
    }
}

/* End of Meteorology.cpp */

