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
    oldPressureVelocity_ = 0.0;
    zeroVectors();
    initMetLoadTypes(optInput);

    //Set reference altitude from reference pressure using ISA
    met::ISA_pAlt(altitudeRef_, pressureRef_);
    std::cout << "Starting pressure of " << pressureRef_ << " Pa. Pressure altitude is " << (altitudeRef_/1000.0) << " km" << std::endl;

    // These should eventually be deleted
    lapseRate_ = -6.5;
    met_depth_ =  200.0;
    if ( optInput.MET_DIURNAL ) {
        throw std::runtime_error("Diurnal temperature variations now disabled");
    }
    diurnalAmplitude_ = 0.0;
    diurnalPert_ = 0.0;
    diurnalPhase_ = 0.0;

    //Declare met input file here, keeps ownership of file to constructor only.
    //Dont want to deal with using a ptr and worry about what functions deference it.
    NcFile dataFile;
    dataFile.open( optInput.MET_FILENAME.c_str(), NcFile::read );

    // Read in the meteorological data
    readMetDataFromFile( dataFile );

    // Set temperature, shear etc. to be the first entry in the file
    updateMetData(0.0);

    // Estimate met data altitude edges and mid-points
    // y should be relative to altitudeRef_, and pressure should be relative to pressureRef_
    // altitudeRef_ in m, pressureRef_ in Pa
    estimateMetDataAltitudes();

    updateSimulationGridProperties();

    // Get the saturation depth
    satdepth_user_ = estimateSatDepth();
} /* End of Meteorology::Meteorology */

void Meteorology::updateSimulationGridProperties(){
    // Develop the simulation grid
    for ( int j = 0; j < ny_; j++ ) {
        altitude_[j] = altitudeRef_ + yCoords_[j];
        altitudeEdges_[j] = altitudeRef_ + yEdges_[j];
    }
    altitudeEdges_[ny_] = altitudeRef_ + yEdges_[ny_];

    // Estimate the pressure at each of these altitudes
    estimateSimulationGridPressures();

    // With the met data altitudes known, estimate the cell mean values for the simulation grid
    interpolateMetToSimulationGrid();

    // Update any derived quantities for the simulation grid
    updateAirMolecDens();
}

void Meteorology::estimateSimulationGridPressures(){
    // Iterate through the meteorolgical grid, calculating the
    // pressure appropriate to the given altitude
    double pressureBase, pressureCeil, pressureTop, zBase, zCeil;
    // This will be copied to pressureBase immediately - so we need it to be the bottom pressure
    pressureCeil = pressureEdgesInit_[0];
    pressureTop = pressureEdgesInit_[altitudeDim_];
    zCeil = altitudeEdgesInit_[0];

    // Index of the simulation grid (NOT the met data grid)
    int iCell = 0;
    bool complete = false;
    bool isEdge = true;
    double zNext = altitudeEdges_[iCell];
    if (zNext <= altitudeInit_[0]) {
        throw std::runtime_error("Lowest altitude is outside of the range of the met data");
    }

    // We are iterating over the spaces between met data points. This is because, between
    // points, we are assuming a constant lapse rate (K/m)
    const double constFactor = physConst::R_Air / physConst::g;
    for (int i = 0; i < altitudeDim_; i++) {
        zBase = zCeil;
        zCeil = altitudeInit_[i];
        while (zNext > zCeil){
            continue;
        }
        pressureBase = pressureCeil;
        pressureCeil = (i == (altitudeDim_)) ? pressureTop : pressureInit_[i];
        double lapseRate = lapseInit_[i];
        // If i < 1 this will cause problems, but that should not be possible
        // as that condition is caught by the check of zNext before the loop
        double tBase = tempInit_[i-1];
        double constFactorLocalInv = 1.0 / (constFactor * lapseRate);
        while (zNext <= zCeil){
            double pressureTarg = pressureBase * pow(tBase/(tBase + lapseRate * (zNext-zBase)),constFactorLocalInv);
            // We always do the lower edge, then the mid point. Therefore
            // only increment the counter when completing the mid point.
            // Once the final edge is complete, don't move on to the next
            // mid point as there isn't one!
            if (isEdge) {
                pressureEdges_[iCell] = pressureTarg;
                if (iCell == ny_) {
                    complete = true;
                    break;
                }
                zNext = altitude_[iCell];
            } else {
                pressure_[iCell] = pressureTarg;
                iCell += 1;
                zNext = altitudeEdges_[iCell];
            }
            isEdge = (!isEdge);
        }
        if (complete) { break; }
    }
}

void Meteorology::applyUpdraft(double dt){
    // Vertical velocity is read from "w" in the netCDF; now assume this to be pressure velocity in Pa s-1
    // For now, we are still using the "reference altitude" rather than a measure of the contrail
    // centroid
    double omega = met::linInterpMetData(altitudeInit_, vertVelocInit_, altitudeRef_);

    // Represent the effect of a pressure velocity by moving all grid edges the same amount (!)
    for(double& pE: pressureEdges_) {
        pE += omega * dt;
    }
    for(double& pC: pressure_) {
        pC += omega * dt;
    }
}

void Meteorology::regenerate( const Vector_1D& yCoord_new, const Vector_1D& yEdges_new, int nx_new ) {
    // Set up new y coordinate, starting from the (generally) increased depth
    double dy_new = yEdges_new[1] - yEdges_new[0];
    int ny_new = yCoord_new.size();
    double alt_y0 = altitudeEdges_[0] + (yEdges_new[0] - yEdges_[0]);
    Vector_1D alt_new(ny_new);
    Vector_1D altEdges_new(ny_new + 1);
    Vector_1D press_new(ny_new);
    Vector_1D pressEdges_new(ny_new + 1);

    std::generate(altEdges_new.begin(), altEdges_new.end(), [&, j = 0] () mutable { return alt_y0 + dy_new * j++; });
    std::generate(alt_new.begin(), alt_new.end(), [&, j = 0] () mutable { return alt_y0 + dy_new * (0.5 + j++); });

    // Grid size has changed - resize all the relevant vectors    
    tempBase_.resize(ny_new);
    tempTotal_ = Vector_2D(ny_new, Vector_1D(nx_new));
    tempPerturbation_ = Vector_2D(ny_new, Vector_1D(nx_new));
    H2O_ = Vector_2D(ny_new, Vector_1D(nx_new));
    airMolecDens_ = Vector_2D(ny_new, Vector_1D(nx_new));
    shear_.resize(ny_new);
    rhi_.resize(ny_new);
    yCoords_ = yCoord_new;
    yEdges_ = yEdges_new;
    ny_ = ny_new;
    nx_ = nx_new;
    altitude_ = alt_new;
    pressure_ = press_new;
    altitudeEdges_ = altEdges_new;
    pressureEdges_ = pressEdges_new;
    vertVeloc_.resize(ny_new);

    // Reinterpolate met data to the new simulation grid
    updateSimulationGridProperties();
}

void Meteorology::Update( const double dt, const double solarTime_h, \
                          const double simTime_h, const double dTrav_x, const double dTrav_y)
{
    //Update altitude if dTrav_y is nonzero 
    //dTrav_x, dTrav_y are for velocities not accounted for in the meteorological vertical velocity.
    if( dTrav_y != 0 ) { 
        throw std::runtime_error("Non-zero dTrav_y not permitted.");
        //for (int jNy = 0; jNy < ny_; jNy++ )
        //    altitude_[jNy] +=  dTrav_y;
        //met::ISA( altitude_, pressure_ );
    }

    // Update the values held in xInit_ - interpolate in time only
    updateMetData(simTime_h);
    // Calculate lapse rates and map altitudes to pressures
    estimateMetDataAltitudes();
} /* End of Meteorology::UpdateMet */

void Meteorology::readMetDataFromFile( const NcFile& dataFile ){
    // Reads in the met data and stores it in memory. This should
    // only ever be called once; the relevant data is stored in
    // the xTimeseriesData_ variables on the original (native)
    // grid, which has a vertical extent of altitudeDim_.
    // This routine does NOT populate the xInit_ variables,
    // which are interpolated in time but not space; and
    // it does NOT populate the x_ variables, which are
    // interpolated in both.
    /*
        Met data format:
        2 dimensions: altitude and time
        RHw, Temp, Shear given as time series.
    */

    /* Identify the length of variables in input file */
    // The dimension is called altitude but could as easily
    // be called level, layer, or z
    altitudeDim_ = dataFile.getDim("altitude").getSize();
    timeDim_ = dataFile.getDim("time").getSize();

    /* Extract pressure and altitude from input file. */
    altitudeInit_.resize(altitudeDim_);
    pressureInit_.resize(altitudeDim_);

    // Read in the pressure variable. This is what we derive our vertical grid from
    // If the altitude variable is present, it is not used
    NcVar pressure_ncVar = dataFile.getVar("pressure");
    pressure_ncVar.getVar(pressureInit_.data());

    //Cache initial pressure and altitude values for generating later met vars.
    for (int i = 0; i < altitudeDim_; i++ ) {
        pressureInit_[i] *= 100.0; //convert from hPa to Pa
    }

    // Resize all variables to match the vertical extent of the met data
    tempInit_.resize(altitudeDim_);
    rhiInit_.resize(altitudeDim_);
    shearInit_.resize(altitudeDim_);
    vertVelocInit_.resize(altitudeDim_);

    // Populate time series data from file
    // Temperature (K), RHi (% - H2O calculated), shear (1/s), and vertical velocity (Pa/s)
    readMetVar(dataFile, "temperature", tempTimeseriesData_, tempLoadType_ == MetVarLoadType::TimeSeries); 
    readMetVar(dataFile, "relative_humidity_ice", rhiTimeseriesData_, rhLoadType_ == MetVarLoadType::TimeSeries); 
    readMetVar(dataFile, "shear", shearTimeseriesData_, shearLoadType_ == MetVarLoadType::TimeSeries); 
    readMetVar(dataFile, "w", vertVelocTimeseriesData_, vertVelocLoadType_ == MetVarLoadType::TimeSeries);
}

void Meteorology::updateMetData(double simTime_h) {
    // Update the vertical profiles of meteorological variables vs pressure
    // Store the old pressure velocity
    if (simTime_h > 0.0) {
        oldPressureVelocity_ = met::linInterpMetData(altitudeInit_, vertVelocInit_, altitudeRef_); // dp/dt
    }
    // The below interpolates in time only
    // We now either interpolate everything or nothing (in time)
    bool timeseries = (tempLoadType_ == MetVarLoadType::TimeSeries);
    tempInit_      = interpMetTimeseriesData(simTime_h, tempTimeseriesData_, timeseries);
    shearInit_     = interpMetTimeseriesData(simTime_h, shearTimeseriesData_, timeseries);
    vertVelocInit_ = interpMetTimeseriesData(simTime_h, vertVelocTimeseriesData_, timeseries);
    rhiInit_       = interpMetTimeseriesData(simTime_h, rhiTimeseriesData_, timeseries);
}

void Meteorology::interpolateMetToSimulationGrid() {
    // Interpolate met data (spatially) from the native vertical grid to the simulation grid
    // This ONLY populates the fields of the met data object - it does not overwrite anything
    // for the main simulation object (i.e. H2O_ here is the meteorological H2O, not the
    // H2O for the main simulation)
    #pragma omp parallel for if (!PARALLEL_CASES)
    for ( int j = 0; j < ny_; j++ ) {
        // Closest grid cell (if using nearest neighbor rather than interpolating)
        int i_Z = met::nearestNeighbor( altitudeInit_, altitude_[j] );
        tempBase_[j] = interpTemp_ ? met::linInterpMetData(altitudeInit_, tempInit_, altitude_[j]) : tempInit_[i_Z];
        shear_[j] = interpShear_ ? met::linInterpMetData(altitudeInit_, shearInit_, altitude_[j]) : shearInit_[i_Z];
        vertVeloc_[j] = interpVertVeloc_ ? met::linInterpMetData(altitudeInit_, vertVelocInit_, altitude_[j]) : vertVelocInit_[i_Z];
        // Water is a bit different - calculate molec/cm3 based on met-data temperature, interpolated (!)
        double rhi_local = interpRH_ ? met::linInterpMetData(altitudeInit_, rhiInit_, altitude_[j]) : rhiInit_[i_Z];
        rhi_[j] = rhi_local;
        double h2o_local = physFunc::RHiToH2O(rhi_local, tempBase_[j]);
        H2O_[j].assign(nx_, h2o_local);
        // Add turbulent fluctuations to temperature if selected
        for ( int i = 0; i < nx_; i++ ) {
            tempTotal_[j][i] =  tempBase_[j] + tempPerturbation_[j][i];
        }
    }
}

void Meteorology::estimateMetDataAltitudes() {
    // Calculates the vertical spacing of pressure levels in the met data.
    // This assumes that temperature varies linearly between levels,
    // and that the pressure varies hydrostatically. We also assume
    // that the air is dry.
    const double constFactor = physConst::R_Air / physConst::g;
    // Calculate lapse rate (K/m) between each point
    // Indexing is a bit tricky: 0 is below the first point, last is above the last point
    lapseInit_.resize(altitudeDim_+1);
    altitudeEdgesInit_.resize(altitudeDim_+1);
    pressureEdgesInit_.resize(altitudeDim_+1);
    for (int i = 0; i < altitudeDim_ - 1; i++) {
        double t0 = tempInit_[i];
        double t1 = tempInit_[i+1];
        double dT = t1 - t0;
        double tempFactor;
        if (std::abs(dT) > 1.0e-10){
            tempFactor = dT/log(t1/t0);
        } else {
            tempFactor = t0;
        }
        double dz = log(pressureInit_[i]/pressureInit_[i+1]) * constFactor * tempFactor;
        lapseInit_[i+1] = dT / dz;
    }
    // Extrapolate below and above the grid
    lapseInit_[0] = lapseInit_[1];
    lapseInit_[altitudeDim_+1] = lapseInit_[altitudeDim_];
    // Calculate the edge spacing, as well as the altitude delta between the lowermost
    // point in the met data and the reference pressure
    double lowerEdge = 0.0;
    double zOffset = 0.0;
    bool zFound = false;
    double pressureBase, pressureCeil, pressureTop;
    // This will be copied to pressureBase immediately - so we need it to be the bottom pressure
    // Do NOT limit to 101325 as the local surface pressure may be greater than that
    pressureCeil = pressureInit_[0] - 0.5 * (pressureInit_[1] - pressureInit_[0]);
    // Pressure at the top of the met data grid, extrapolated
    pressureTop = std::max(0.0,pressureInit_[altitudeDim_-1] + 0.5 * (pressureInit_[altitudeDim_-1] - pressureInit_[altitudeDim_-2]));
    pressureEdgesInit_[0] = pressureCeil;
    altitudeEdgesInit_[0] = lowerEdge;

    for (int i = 0; i < altitudeDim_; i++) {
        pressureBase = pressureCeil;
        pressureCeil = (i == (altitudeDim_-1)) ? pressureTop : 0.5*(pressureInit_[i+1] + pressureInit_[i]);
        pressureEdgesInit_[i+1] = pressureCeil;
        double lowerLapse = lapseInit_[i];
        double upperLapse = lapseInit_[i+1];
        double pMid = pressureInit_[i];
        double tMid = tempInit_[i];
        // Two calculations: one for the lower half of the cell, one for the upper half
        // Lapse rate changes these are not equal
        double dzLower = hydrostaticDeltaAltitude(lowerLapse,tMid,pMid,pressureBase);
        double dzUpper = -1.0*hydrostaticDeltaAltitude(upperLapse,tMid,pMid,pressureCeil);
        altitudeInit_[i] = lowerEdge + dzLower;
        // Does this "cell" contain the reference pressure? If so, figure out how far 
        // into it that pressure occurs
        if ((!zFound) && (pressureBase >= pressureRef_) && (pressureCeil < pressureRef_)){
            if (pressureRef_ > pMid){
                // Lower half of the cell
                zOffset = lowerEdge + hydrostaticDeltaAltitude(lowerLapse,tMid,pMid,pressureRef_);
            } else {
                // Upper half of the cell
                zOffset = lowerEdge + dzLower - hydrostaticDeltaAltitude(upperLapse,tMid,pMid,pressureRef_);
            }
            zFound = true;
        }
        // Update to the top of this cell (and bottom of the next)
        lowerEdge += dzLower + dzUpper;
        altitudeEdgesInit_[i+1] = lowerEdge;
    }
    if (!zFound) {
        throw std::runtime_error("Reference pressure outside meteorological data bounds.");
    }
    // Now update all the altitude centers so that they place the reference pressure and 
    // altitude in the same location on the grid
    for (int i = 0; i < altitudeDim_; i++) {
        altitudeInit_[i] += (altitudeRef_ - zOffset);
        altitudeEdgesInit_[i] += (altitudeRef_ - zOffset);
    } 
}

// Not a class method
double hydrostaticDeltaAltitude(double lapseRate, double refTemperature, double refPressure, double targPressure) {
    // Given a lapse rate, calculates the distance from the reference point to the given pressure
    // Lapse rate is in K/m. If the target pressure is greater than the reference pressure the
    // result will be positive. NB: The "reference" pressure here means the pressure given at the
    // mid-point of the cell, and is NOT the simulation reference pressure.
    const double constFactor = physConst::R_Air / physConst::g;
    double dz;
    // Singularity at lapse = 0; use the limiting form for very small/zero lapse rates
    if (std::abs(lapseRate) < 1.0e-10){
        dz = -1.0 * refTemperature * constFactor * log(refPressure/targPressure);
    } else {
        dz = (refTemperature/lapseRate) * (1.0 - pow(refPressure/targPressure,lapseRate*constFactor));
    }
    return dz;
}

void Meteorology::recalculateSimulationGrid(){
    // Given a set of pressure edges, calculate their altitudes.
    const double constFactor = physConst::R_Air / physConst::g;
    double pressureBase, pressureCeil, pressureTop;
    // Start from the bottom of the met data grid
    double lowerEdge = altitudeEdgesInit_[0];
    pressureCeil = pressureEdgesInit_[0];

    // Pressures for which we are trying to calculate the altitude
    int iSim = 0;
    int ny = pressure_.size();
    double pressureTarg = pressureEdges_[0];
    altitudeEdges_.resize(ny+1);
    yEdges_.resize(ny+1);

    // First entry is for the space between the lowermost
    // met data pressure edge and the first met data pressure.
    // After that, we loop over the spaces between the met data
    // points, as the lapse rate is constant over this space
    for (int i = 0; i < altitudeDim_; i++) {
        pressureBase = pressureCeil;
        pressureCeil = pressureInit_[i];
        double lapseRate = lapseInit_[i];
        double tMid = tempInit_[i];
        // Loop over all entries which are between this pressure and the prior one
        // The same lapse rate (remembering lapse rates are indexed like pressure edges,
        // not like pressure mid-points) will apply throughout
        while ((pressureTarg >= pressureCeil) && (pressureTarg < pressureBase)){
            double dz = hydrostaticDeltaAltitude(lapseRate,tMid,pressureTarg,pressureBase);
            altitudeEdges_[iSim] = dz + lowerEdge;
            yEdges_[iSim] = dz + lowerEdge - altitudeRef_;
            iSim++;
            if (iSim > ny) { break; }
            pressureTarg = pressureEdges_[iSim];
        }
        if (iSim > ny) { break; }
        lowerEdge = altitudeInit_[i];
    }
    // Don't bother with updating altitude_ - it will not be used again
    // before the remapping is performed
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

double Meteorology::estimateSatDepth() {
    // Calculates the saturation depth
    int i_Z = met::nearestNeighbor( pressure_, ambParams_.press_Pa);
    double dy = yCoords_[1] - yCoords_[0];
    try {
        return met::satdepth_calc(rhi_, altitude_, i_Z, std::abs(yCoords_[0]) + dy/2);
    }
    catch (std::out_of_range &e){
        std::cout << "WARNING: end of initial domain does not cover saturation depth." << std::endl;
        return 1000.0;
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
    // Use ideal gas law to estimate molec/cm3 of air
    // We are here including the temperature perturbation,
    // which may or may not be correct depending on the
    // situation - caveat emptor
    const double invkB = 1.00E-06 / physConst::kB;
    const double constFactor = physConst::Na / (physConst::g * MW_Air);
    #pragma omp parallel for
    for ( int jNy = 0; jNy < ny_; jNy++ ) {
        double integralND = constFactor * (pressureEdges_[jNy] - pressureEdges_[jNy+1]) / (altitudeEdges_[jNy+1] - altitudeEdges_[jNy]);
        for ( int iNx = 0; iNx < nx_; iNx++ ) {
            airMolecDens_[jNy][iNx] = integralND;
        }
    }
}

/* End of Meteorology.cpp */

