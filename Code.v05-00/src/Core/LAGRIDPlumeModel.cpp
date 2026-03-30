#include <chrono>
#include <filesystem>
#include <memory>
#include "AIM/Settling.hpp"
#include "Util/PlumeModelUtils.hpp"
#include "Core/Status.hpp"
#include "KPP/KPP_Parameters.h"
#include "Core/LAGRIDPlumeModel.hpp"
#include "EPM/Models/Original.hpp"

using physConst::Na, physConst::PI, physConst::R_Air;

/*
These global variables are required because of the KPP auto-generated code.
In the LAGRID Plume Model (no plume chemistry), the KPP is used only when
spinning up the background ambient conditions during the initialization
step of the EPM. For now declare these variables at the top level.
*/
double RCONST[NREACT];       /* Rate constants (global) */
double NOON_JRATES[NPHOTOL]; /* Noon-time photolysis rates (global) */
double PHOTOL[NPHOTOL];      /* Photolysis rates (global) */
double HET[NSPEC][3];        /* Heterogeneous chemistry rates (global) */
double TIME;                 /* Current integration time (global) */
double SZA_CST[3];           /* Require this for adjoint integration */

LAGRIDPlumeModel::LAGRIDPlumeModel(const OptInput &optInput, Input &input) :
    optInput_(optInput), input_(input),
    numThreads_(optInput.SIMULATION_OMP_NUM_THREADS),
    jetA_(Fuel("C12H24")),
    yCoords_(optInput_.ADV_GRID_NY),
    yEdges_(optInput_.ADV_GRID_NY + 1)
{
    /* Multiply by 500 since it gets multiplied by 1/500 within the Emission object ... */
    jetA_.setFSC(input.EI_SO2() * 500.0);

    // Set up initial grid coordinates: required for meteorology setup, which
    // is needed for the EPM.
    double dy = (optInput_.ADV_GRID_YLIM_UP + optInput_.ADV_GRID_YLIM_DOWN) /
                optInput_.ADV_GRID_NY;
    double y0 = -optInput_.ADV_GRID_YLIM_DOWN;
    std::generate(yEdges_.begin(), yEdges_.end(),
                  [dy, y0, j = 0]() mutable { return y0 + dy * j++; });
    std::generate(yCoords_.begin(), yCoords_.end(),
                  [dy, y0, j = 0]() mutable { return y0 + dy * (0.5 + j++); });

    met_ = Meteorology(optInput_, input.pressure_Pa(), yCoords_, yEdges_);
    std::cout << "Temperature      = " << met_.tempRef() << " K" << std::endl;
    std::cout << "RHw              = " << met_.rhwRef() << " %" << std::endl;
    std::cout << "RHi              = " << met_.rhiRef() << " %" << std::endl;
    std::cout << "Saturation depth = " << met_.satdepthEstimate() << " m" << std::endl;
    
    // Set the met variables in the input object
    input.set_temperature_K(met_.tempRef());
    input.set_relHumidity_w(met_.rhwRef());
    input.set_relHumidity_i(met_.rhiRef());
    input_ = input;

    simVars_ = MPMSimVarsWrapper(input_, optInput_, met_.satdepthEstimate());
    timestepVars_ = TimestepVarsWrapper(input_, optInput_);
    aircraft_ = Aircraft(input_, optInput_.SIMULATION_INPUT_ENG_EI);
    EI_ = Emission(aircraft_.engine(), jetA_, input_.EI_SO2TOSO4());
    WV_exhaust_ = EI_.getH2O() * aircraft_.fuel_per_dist();
    std::cout << "H2O EI     = " << EI_.getH2O() << std::endl;

    timestepVars_.setTimeArray(PlumeModelUtils::BuildTime (
            timestepVars_.tInitial_s, timestepVars_.tFinal_s, timestepVars_.dt));
    createOutputDirectories();
}

SimStatus LAGRIDPlumeModel::runFullModel() {
    auto start = std::chrono::high_resolution_clock::now();
    omp_set_num_threads(numThreads_);

    std::variant<EPM::Output, SimStatus> epm_result = runEPM();
    if (std::holds_alternative<SimStatus>(epm_result)) {
        return std::get<SimStatus>(epm_result);
    }
    EPM::Output epmOut = std::get<EPM::Output>(epm_result);

    //Initialize aerosol into grid and init H2O
    initializeGrid(epmOut);
    initH2O();
    saveTSAerosol();

    //Setup settling velocities
    if ( simVars_.GRAVSETTLING ) {
        vFall_ = AIM::SettlingVelocity( iceAerosol_.getBinCenters(), \
                                       met_.tempRef(), simVars_.pressure_Pa );
    }

    bool EARLY_STOP = false;
    SimStatus status = SimStatus::Incomplete;
    //Start time loop
    while ( timestepVars_.curr_Time_s < timestepVars_.tFinal_s ) {
        #ifdef ENABLE_TIMING
        auto timestep_start = std::chrono::high_resolution_clock::now();
        #endif

        /* Print message */
        std::cout << "\n";
        std::cout << "\n - Time step: " << timestepVars_.nTime + 1 << " out of " << timestepVars_.timeArray.size();
        std::cout << "\n -> Solar time: " << std::fmod( timestepVars_.curr_Time_s/3600.0, 24.0 ) << " [hr]" << std::endl;

        // Run Transport
        std::cout << "Running Transport..." << std::endl;
        bool timeForTransport = (simVars_.TRANSPORT && (timestepVars_.nTime == 0 || timestepVars_.checkTimeForTransport()));
        if (timeForTransport) {

            #ifdef ENABLE_TIMING
            auto transport_start = std::chrono::high_resolution_clock::now();
            #endif

            runTransport(timestepVars_.TRANSPORT_DT);

            #ifdef ENABLE_TIMING
            auto transport_end = std::chrono::high_resolution_clock::now();
            auto transport_duration = std::chrono::duration_cast<std::chrono::milliseconds>(transport_end - transport_start);
            std::cout << "  Ran transport in " << transport_duration.count() << " ms" << std::endl;
            #endif
        }

        /*  With LAGRID remapping every transport timestep, it fundamentally only makes physical sense to update
            the temperature perturbations at the same interval as the transport timestep. Turbulence timestep is one
            tool used to tune the intensity of the simulated turbulence, but we can also just vary the amplitude.
        */
        if (simVars_.TEMP_PERTURB){
            met_.updateTempPerturb();
        }

        solarTime_h_ = ( timestepVars_.curr_Time_s + timestepVars_.TRANSPORT_DT / 2 ) / 3600.0;
        simTime_h_ = ( timestepVars_.curr_Time_s + timestepVars_.TRANSPORT_DT / 2 - timestepVars_.timeArray[0] ) / 3600;

        // Run Ice Growth
        if (simVars_.ICE_GROWTH && timestepVars_.checkTimeForIceGrowth()) {
            std::cout << "Running ice growth..." << std::endl;

            #ifdef ENABLE_TIMING
            auto icegrowth_start = std::chrono::high_resolution_clock::now();
            #endif

            timestepVars_.lastTimeIceGrowth = timestepVars_.curr_Time_s + timestepVars_.dt;
            iceAerosol_.Grow( timestepVars_.ICE_GROWTH_DT, H2O_, met_.Temp(), met_.Press());

            #ifdef ENABLE_TIMING
            auto icegrowth_end = std::chrono::high_resolution_clock::now();
            auto icegrowth_duration = std::chrono::duration_cast<std::chrono::milliseconds>(icegrowth_end - icegrowth_start);
            std::cout << "  Ran ice growth in " << icegrowth_duration.count() << " ms" << std::endl;
            #endif
        }

        #ifdef ENABLE_TIMING
        auto tracer_start = std::chrono::high_resolution_clock::now();
        #endif

        // Update the tracer of contrail influence to include all locations where we have ice
        // Set it to 1 when there's at least 1 particle per m3 
        auto number = iceAerosol_.TotalNumber();
        for (std::size_t j=0; j<yCoords_.size(); j++){
            for (std::size_t i=0; i<xCoords_.size(); i++){
                Contrail_[j][i] = std::max(0.0,std::min(1.0,Contrail_[j][i] + number[j][i]*1.0e6));
            }
        }

        #ifdef ENABLE_TIMING
        auto tracer_end = std::chrono::high_resolution_clock::now();
        auto tracer_duration = std::chrono::duration_cast<std::chrono::milliseconds>(tracer_end - tracer_start);
        std::cout << "  Ran tracer update in " << tracer_duration.count() << " ms" << std::endl;

        auto met_update_start = std::chrono::high_resolution_clock::now();
        #endif

        // Create the mask of cells to retain, and terminate if none left
        auto dataMask = ContrailMask(0.90);
        auto& mask = dataMask.first;
        auto& maskInfo = dataMask.second;
        if (maskInfo.count == 0){
            std::cout << "No remaining grid cells marked as contrail-containing." << std::endl;
            EARLY_STOP = true;
            status = SimStatus::Complete;
            break;
        }

        // Convert all quantities from #/cm3 to #/#
        // delta-P/delta-Z is proportional to number density, and accounts
        // for temperature variation (because the calculation of yEdge
        // included it). This uses the hydrostatic assumption:
        //    dp/dz = -rho*g = -(n/V)Mg
        Vector_3D& pdfRef = iceAerosol_.getPDF();
        auto pressureEdges = met_.PressEdges();
        double localND;
        #pragma omp parallel for
        for (std::size_t j=0; j<yCoords_.size(); j++){
            localND = (pressureEdges[j] - pressureEdges[j+1])/(yEdges_[j+1] - yEdges_[j]);
            for (std::size_t i=0; i<xCoords_.size(); i++){
                Contrail_[j][i] = Contrail_[j][i] / localND; // parts per trillion
                H2O_[j][i] = H2O_[j][i] / localND; // parts per trillion
                for(UInt n = 0; n < iceAerosol_.getNBin(); n++) {
                    pdfRef[n][j][i] = pdfRef[n][j][i] / localND;
                }
            }
        }

        // Run vertical advection to get the new pressure edges
        if (std::abs(met_.lastOmega()) > 1.0e-10){
            met_.applyUpdraft(timestepVars_.TRANSPORT_DT);
        }

        // Read in the meteorology for the next time step, interpolating in time 
        // This only affects the xInit_ values in the meteorology object, but does
        // update altitudeInit_ and altitudeEdgesInit_ (pressureInit_ unchanged)
        std::cout << "Updating Met..." << std::endl;
        met_.Update( timestepVars_.TRANSPORT_DT, solarTime_h_, simTime_h_);


        #ifdef ENABLE_TIMING
        auto met_update_end = std::chrono::high_resolution_clock::now();
        auto met_update_duration = std::chrono::duration_cast<std::chrono::milliseconds>(met_update_end - met_update_start);
        std::cout << "  Ran met update in " << met_update_duration.count() << " ms" << std::endl;

        auto regrid_start = std::chrono::high_resolution_clock::now();
        #endif
        
        // Determine the altitude of each pressure edge after vertical advection and the 
        // application of new temperatures
        met_.recalculateSimulationGrid();
        
        // Resynchronize the y coordinates
        yEdges_ = met_.yEdges();
        //yCoords_ = met_.yCoords(); // We didn't update met_.yCoords

        //Remap the grid to account for changes in shape due to vertical advection and the growth of the contrail
        std::cout << "Remapping... " << std::endl;
        remapAllVars(timestepVars_.TRANSPORT_DT, mask, maskInfo);

        #ifdef ENABLE_TIMING
        auto regrid_end = std::chrono::high_resolution_clock::now();
        auto regrid_duration = std::chrono::duration_cast<std::chrono::milliseconds>(regrid_end - regrid_start);
        std::cout << "  Ran remapping in " << regrid_duration.count() << " ms" << std::endl;
        #endif

        Vector_2D areas = VectorUtils::cellAreas(xEdges_, yEdges_);
        double numparts = iceAerosol_.TotalNumber_sum(areas);
        std::cout << "Num Particles: " << numparts << std::endl;
        std::cout << "Ice Mass: " << iceAerosol_.TotalIceMass_sum(areas) << std::endl;
        if(numparts / initNumParts_ < 1e-5) {
            std::cout << "Less than 0.001% of the particles remain, stopping sim" << std::endl;
            EARLY_STOP = true;
        }

        // Advance time
        timestepVars_.curr_Time_s += timestepVars_.dt;
        timestepVars_.nTime++;

        #ifdef ENABLE_TIMING
        auto save_start = std::chrono::high_resolution_clock::now();
        #endif

        saveTSAerosol();

        #ifdef ENABLE_TIMING
        auto save_end = std::chrono::high_resolution_clock::now();
        auto save_duration = std::chrono::duration_cast<std::chrono::milliseconds>(save_end - save_start);
        std::cout << "  Ran save in " << save_duration.count() << " ms" << std::endl;

        auto timestep_end = std::chrono::high_resolution_clock::now();
        auto timestep_duration = std::chrono::duration_cast<std::chrono::milliseconds>(timestep_end - timestep_start);
        std::cout << "Ran timestep in " << timestep_duration.count() << " ms" << std::endl;
        #endif

        if(EARLY_STOP) {
            status = SimStatus::Complete;
            break;
        }
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
    std::cout << "APCEMM LAGRID Plume Model Run Finished! Run time: " << duration.count() << "ms" << std::endl;
    return status;
}

std::variant<EPM::Output, SimStatus> LAGRIDPlumeModel::runEPM() {
    // Calculate the number of emitted soot particles
    double EI_getSootRad = EI_.getSootRad();
    const double volParticle  = 4.0 / 3.0 * physConst::PI * EI_getSootRad * EI_getSootRad * EI_getSootRad; //EI_SootRad in m -> volume in m3
    const double massParticle = volParticle * physConst::RHO_SOOT * 1.0E+03; //Gives mass of a particle in grams
    const double EI_icenum = EI_.getSoot() / massParticle; /* [#/kg_fuel] */
    const double N0 = EI_icenum * aircraft_.fuel_per_dist();
    std::cout << "EI_icenum prejet: " << EI_icenum << " [#/kg]" << std::endl;
    std::cout << "Emitted soot: " << N0 << " [#/m]" << std::endl;

    std::unique_ptr<EPM::Models::Base> epm = EPM::make_epm(
        optInput_, input_, aircraft_, EI_, met_, simVars_);
    std::variant<EPM::Output, SimStatus> epmResult = epm->run();
    if (std::holds_alternative<SimStatus>(epmResult)) {
        return std::get<SimStatus>(epmResult);
    }
    EPM::Output &epmOutput = std::get<EPM::Output>(epmResult);
    epmOutput.write(optInput_.SIMULATION_OUTPUT_FOLDER + "/epm-output.nc");

    /* The output area from the EPM is principally used just to scale from
    number densities (eg #/cm3) to totals (eg #/m). Increasing the area
    is the easiest way to represent the fact that multiple engines are in
    use - recalling that the EPM is only actually run for one engine. */
    epmOutput.area *= aircraft_.EngNumber();


    // Calculate the number of ice particles past the jet regime
    double N_postjet = epmOutput.IceAer.Moment(0) * epmOutput.area * 1e6;

    if (optInput_.ADV_EP_N_POSTJET_OVERRIDE) {
        epmOutput.IceAer.scalePdf(optInput_.ADV_EP_N_POSTJET / N_postjet);
        N_postjet = epmOutput.IceAer.Moment(0) * epmOutput.area * 1e6; // recalculate after scaling
        std::cout << "Overriding post-jet ice crystal count to: " << N_postjet;
        std::cout << " [#/m] (Requested value:" << optInput_.ADV_EP_N_POSTJET << " [#/m])" << std::endl;
    } else {
        std::cout << "Post-jet ice crystal count: " << N_postjet << " [#/m]" << std::endl;
    }
    std::cout << "EI_icenum postjet: " << N_postjet / aircraft_.fuel_per_dist() << " [#/kg]" << std::endl;

    // Run vortex phase parameterization
    const double N0_ref = optInput_.ADV_EP_N_REF;
    const double wingspan_ref = optInput_.ADV_EP_WINGSPAN_REF;
    const double icenum_survfrac = aircraft_.VortexLosses(N_postjet, WV_exhaust_, N0_ref, wingspan_ref);
    std::cout << "Parameterized vortex sinking survival fraction: " << icenum_survfrac
                << std::endl;
    if (icenum_survfrac <= 0.00E+00) {
        std::cout << "EndSim: vortex sinking" << std::endl;
        return SimStatus::NoSurvivalVortex;
    }
    epmOutput.IceAer.scalePdf(icenum_survfrac);

    std::cout << "Post-vortex ice particle count: " << epmOutput.IceAer.Moment(0) * epmOutput.area * 1e6 << " [#/m]" << std::endl;

    return epmOutput;
}

void LAGRIDPlumeModel::initializeGrid(const EPM::Output &epmOut) {
    auto& epmIceAer = epmOut.IceAer;

    /* TODO: Fix the initial boundaries */
    int nx_init = optInput_.ADV_GRID_NX;
    double dx_init = (optInput_.ADV_GRID_XLIM_RIGHT + optInput_.ADV_GRID_XLIM_LEFT)/nx_init;
    double x0 = -optInput_.ADV_GRID_XLIM_LEFT;

    int ny_init = optInput_.ADV_GRID_NY;
    double dy = (optInput_.ADV_GRID_YLIM_DOWN + optInput_.ADV_GRID_YLIM_UP)/ny_init;
    double y0 = -optInput_.ADV_GRID_YLIM_DOWN;

    yCoords_ = Vector_1D(ny_init);
    yEdges_ =  Vector_1D(ny_init + 1);
    std::generate(yEdges_.begin(), yEdges_.end(), [dy, y0, j = 0] () mutable { return y0 + dy * j++; } );
    std::generate(yCoords_.begin(), yCoords_.end(), [dy, y0, j = 0] () mutable { return y0 + dy * (0.5 + j++); } );

    xCoords_ = Vector_1D(nx_init);
    xEdges_ = Vector_1D(nx_init + 1);
    std::generate(xEdges_.begin(), xEdges_.end(), [dx_init, x0, i = 0] () mutable { return x0 + dx_init * i++; } );
    std::generate(xCoords_.begin(), xCoords_.end(), [dx_init, x0, i = 0] () mutable { return x0 + dx_init * (0.5 + i++); } );

    iceAerosol_ = AIM::Grid_Aerosol(nx_init, yCoords_.size(), epmIceAer.getBinCenters(), epmIceAer.getBinEdges(), 0, 1, 1.6);

    Vector_3D pdf_init;
    
    //Initialize the mature plume as a rectangle
    const double initWidth = aircraft_.vortex().width_rect_mature();
    const double initDepth = aircraft_.vortex().depth_mature();
    const double xcenter = 0;
    const double ycenter = -aircraft_.vortex().z_center();

    for (UInt n = 0; n < iceAerosol_.getNBin(); n++) {
        double logBinRatio = log(iceAerosol_.getBinEdges()[n+1] / iceAerosol_.getBinEdges()[n]);
        double EPM_nPart_bin = epmIceAer.binMoment(n) * epmOut.area * logBinRatio;
        pdf_init.push_back( LAGRID::initVarToGridRectangular(EPM_nPart_bin, xEdges_, yEdges_, xcenter, ycenter, initWidth, initDepth, logBinRatio) );

    }
    iceAerosol_.updatePdf(std::move(pdf_init));
    Vector_2D areas = VectorUtils::cellAreas(xEdges_, yEdges_);
    initNumParts_ = iceAerosol_.TotalNumber_sum(areas);

    std::cout << "Mature Plume Initial Num Particles: " << initNumParts_ << std::endl;
    std::cout << "Mature Plume Initial Ice Mass: " << iceAerosol_.TotalIceMass_sum(areas) << std::endl;
}

void LAGRIDPlumeModel::initH2O() {
    Contrail_ = Vector_2D(yCoords_.size(), Vector_1D(xCoords_.size()));
    H2O_ = met_.H2O_field();

    //Add emitted plume H2O. This function is called after releasing the initial crystals into the grid,
    //so we can use that as a "mask" for where to emit the H2O.

    auto maskInfo = VectorUtils::mask(iceAerosol_.TotalNumber(), [](double val) { return val > 1e-4; } );
    auto& mask = maskInfo.first;
    int nonMaskCount = maskInfo. second;

    auto areas = VectorUtils::cellAreas(xEdges_, yEdges_);
    const double icemass = iceAerosol_.TotalIceMass_sum(areas);

    double mass_WV = WV_exhaust_ - icemass;
    double E_H2O = mass_WV / (MW_H2O * 1e3) * Na;

    // Spread the emitted water evenly over the cells that contain ice crystals
    auto localPlumeEmission = [&](std::size_t j, std::size_t i) -> double {
        if(mask[j][i] == 0) return 0;
        return E_H2O * 1.0E-06 / ( nonMaskCount * areas[j][i] );
    };
    //Set the H2O field to the plume H2O, and mark locations with 1 or 0 to indicate
    //the presence of contrail ice
    for(std::size_t j = 0; j < mask.size(); j++) {
        for(std::size_t i = 0; i < mask[0].size(); i++) {
            if(mask[j][i] == 1) {
                double localPlumeH2O = localPlumeEmission(j, i);
                H2O_[j][i] += localPlumeH2O;
                Contrail_[j][i] = 1.0;
            } else {
                Contrail_[j][i] = 0.0;
            }
        }
    }
}

void LAGRIDPlumeModel::updateDiffVecs() {
    double dh_enhanced, dv_enhanced;
    // Update Diffusion
    PlumeModelUtils::DiffParam( timestepVars_.curr_Time_s - timestepVars_.tInitial_s + timestepVars_.TRANSPORT_DT / 2.0,
                                dh_enhanced, dv_enhanced, input_.horizDiff(), input_.vertiDiff() );
    auto number = iceAerosol_.TotalNumber();
    auto num_max = VectorUtils::max(number);
    
    diffCoeffX_ = Vector_2D(yCoords_.size(), Vector_1D(xCoords_.size()));
    diffCoeffY_ = Vector_2D(yCoords_.size(), Vector_1D(xCoords_.size()));
    
    #pragma omp parallel for
    for(std::size_t j = 0; j < yCoords_.size(); j++) {
        for(std::size_t i = 0; i < xCoords_.size(); i++) {
            diffCoeffX_[j][i] = number[j][i] > num_max * 1e-4 ? dh_enhanced : input_.horizDiff();
            diffCoeffY_[j][i] = number[j][i] > num_max * 1e-4 ? dv_enhanced : input_.vertiDiff();
        }
    }
}
void LAGRIDPlumeModel::runTransport(double timestep) {

    #ifdef ENABLE_TIMING
    auto start = std::chrono::high_resolution_clock::now();
    #endif
    //Update the zero bc to reflect grid size changes
    auto ZERO_BC = FVM_ANDS::bcFrom2DVector(iceAerosol_.getPDF()[0], true);

    //TODO: Implement height dependent shear. For now, just taking shear of y coordinate with highest xOD to avoid bugs.
    auto xOD = iceAerosol_.xOD(Vector_1D(xCoords_.size(), xCoords_[1] - xCoords_[0]));
    int maxIdx = 0;
    for (std::size_t i = 0; i < xOD.size(); i++) {
        if(xOD[i] > xOD[maxIdx]) maxIdx = i;
    }
    shear_rep_ = met_.shear(maxIdx);

    const FVM_ANDS::AdvDiffParams fvmSolverInitParams(0, 0, shear_rep_, input_.horizDiff(), input_.vertiDiff(), timestepVars_.TRANSPORT_DT);
    const FVM_ANDS::BoundaryConditions ZERO_BC_INIT = FVM_ANDS::bcFrom2DVector(iceAerosol_.getPDF()[0], true);

    #ifdef ENABLE_TIMING
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "      Initialized transport in " << duration.count() << " ms" << std::endl;

    start = std::chrono::high_resolution_clock::now();
    #endif

    updateDiffVecs();

    #ifdef ENABLE_TIMING
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "      Updated DiffVecs in " << duration.count() << " ms" << std::endl;

    start = std::chrono::high_resolution_clock::now();
    #endif

    // Convert this once for all ~38 aerosol size bins as it reused for all of them
    // The real solution is to use Eigen::VectorXd from the start...
    Eigen::VectorXd H2O_eigen = FVM_ANDS::std2dVec_to_eigenVec(H2O_);
    Eigen::VectorXd diffCoeffX_eigen = FVM_ANDS::std2dVec_to_eigenVec(diffCoeffX_);
    Eigen::VectorXd diffCoeffY_eigen = FVM_ANDS::std2dVec_to_eigenVec(diffCoeffY_);

    #ifdef ENABLE_TIMING
    end = std::chrono::high_resolution_clock::now();
    auto duration_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "      Copies vec2Eigen in " << duration_us.count() << " us" << std::endl;

    start = std::chrono::high_resolution_clock::now();
    #endif


    /* PRECOMPUTE THE COEFFICIENT MATRIX FOR ALL AEROSOL BINS */

    // Cache the coefficients matrix because it is expensive to build and is identical
    // across all aerosol bins.
    // This is only true for operator splitting because the matrix is only used
    // for diffusion. The boolean here is to make this explicit, it gets compiled away.
    bool useOperatorSplit = true;
    // Need to instantiate a solver instance to compute the matrix, this solver is not used otherwise
    FVM_ANDS::FVM_Solver template_solver(fvmSolverInitParams, xCoords_, yCoords_, ZERO_BC_INIT, H2O_eigen);
    std::shared_ptr<const Eigen::SparseMatrix<double, Eigen::RowMajor>> shared_totalCoefMatrixPtr;
    if (useOperatorSplit) {
        template_solver.updateTimestep(timestep);

        #ifdef ENABLE_TIMING
        auto start_diff = std::chrono::high_resolution_clock::now();
        #endif

        template_solver.updateDiffusion(diffCoeffX_eigen, diffCoeffY_eigen);

        #ifdef ENABLE_TIMING
        auto end_diff = std::chrono::high_resolution_clock::now();
        duration_us = std::chrono::duration_cast<std::chrono::microseconds>(end_diff - start_diff);
        std::cout << "      updateDiffusion " << duration_us.count() << " us" << std::endl;
        #endif

        #ifdef ENABLE_TIMING
        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "      Template solver instantiation " << duration.count() << " ms" << std::endl;

        start = std::chrono::high_resolution_clock::now();
        #endif

        template_solver.buildCoeffMatrix(useOperatorSplit);

        #ifdef ENABLE_TIMING
        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "      Template build matrix creation " << duration.count() << " ms" << std::endl;

        start = std::chrono::high_resolution_clock::now();
        #endif

        shared_totalCoefMatrixPtr = template_solver.coefMatrixPtr();
    }

    //Transport the Ice Aerosol PDF

    /* PRE-ALLOCATE A SOLVER POOL */

    // Pre-allocate one solver per thread to creating a new solver each iteration
    // Each thread will reuse its solver for all its assigned bins.
    // Use vector<unique_ptr<solver>> and not directly a vector<solver> because of unique_pointers
    // in the solver object which cannot be moved...
    // This forces the instantiation of the solvers to be serial which is a slowdown
    // relative to the naive each iteration creates a new solver...
    // This slowdown is mostly compensated for by the gains of caching the matrix:
    // For 8 cpus there is a slight slowdown between naive and preallocating solvers:
    // -> preallocating solvers + caching the matrix (this)
    // -> Naive each iteration of the aerosol bins loop creates its own solver (previous)
    // For 1 cpu however we get a 15-35% speedup on transport (which gets better as the grid is larger
    // e.g. the older the contrail the better the speedup for each transport step)
    std::vector<std::unique_ptr<FVM_ANDS::FVM_Solver>> solver_pool;
    solver_pool.reserve(numThreads_);
    
    for (int i = 0; i < numThreads_; ++i) {
        solver_pool.emplace_back(std::make_unique<FVM_ANDS::FVM_Solver>(
            fvmSolverInitParams, xCoords_, yCoords_, ZERO_BC_INIT, H2O_eigen));
        solver_pool.back()->updateTimestep(timestep);
        solver_pool.back()->updateDiffusion(diffCoeffX_eigen, diffCoeffY_eigen);
        if (useOperatorSplit) {
            solver_pool.back()->setPrebuiltMatrix(shared_totalCoefMatrixPtr);
        }
    }

    #ifdef ENABLE_TIMING
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "      Init solver pool " << duration.count() << " ms" << std::endl;

    start = std::chrono::high_resolution_clock::now();
    #endif
    
    #pragma omp parallel for default(shared) schedule(dynamic)
    for ( UInt n = 0; n < iceAerosol_.getNBin(); n++ ) {
        /* Transport particle number and volume for each bin and
            * recompute centers of each bin for each grid cell
            * accordingly */
        
        // Fetch the solver for this thread
        int tid = omp_get_thread_num();
        auto& solver = *solver_pool[tid];
        
        //Update solver params for this bin
        solver.updateAdvection(0, -vFall_[n], shear_rep_);

        //passing in "false" to the "parallelAdvection" param to not spawn more threads
        solver.operatorSplitSolve2DVec(iceAerosol_.getPDF()[n], ZERO_BC, false);
    }

    #ifdef ENABLE_TIMING
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "      Ran ice aerosol transport in " << duration.count() << " ms" << std::endl;

    start = std::chrono::high_resolution_clock::now();
    #endif

    //Transport H2O
    {   
        //Dont use enhanced diffusion on the H2O (and zero settling velocity)
        FVM_ANDS::FVM_Solver solver(fvmSolverInitParams, xCoords_, yCoords_, ZERO_BC_INIT, H2O_eigen);
        solver.updateTimestep(timestep);
        solver.updateDiffusion(input_.horizDiff(), input_.vertiDiff());
        solver.updateAdvection(0, 0, shear_rep_);

        // Calculate diffusion relative to a vertically-varying background H2O field
        // This prevents APCEMM from smoothing out pre-existing meteorological gradients
        // which will remain in the background/boundary conditions.
        Vector_2D H2O_Delta;
        H2O_Delta = Vector_2D(yCoords_.size(), Vector_1D(xCoords_.size()));
        auto H2O_Background = met_.H2O_field();
        for (std::size_t j=0; j<yCoords_.size(); j++){
            for (std::size_t i=0; i<xCoords_.size(); i++){
                H2O_Delta[j][i] = H2O_[j][i] - H2O_Background[j][i];
            }
        }
        // BC is zero, since we're calculating the difference relative to background.
        solver.operatorSplitSolve2DVec(H2O_Delta, ZERO_BC);
        for (std::size_t j=0; j<yCoords_.size(); j++){
            for (std::size_t i=0; i<xCoords_.size(); i++){
                H2O_[j][i] = H2O_Delta[j][i] + H2O_Background[j][i];
            }
        }
    }

    #ifdef ENABLE_TIMING
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "      Ran H2O transport in " << duration.count() << " ms" << std::endl;

    start = std::chrono::high_resolution_clock::now();
    #endif

    //Transport the contrail tracer
    {   
        //Identical settings to H2O
        FVM_ANDS::FVM_Solver solver(fvmSolverInitParams, xCoords_, yCoords_, ZERO_BC_INIT, FVM_ANDS::std2dVec_to_eigenVec(Contrail_));
        solver.updateTimestep(timestep);
        solver.updateDiffusion(input_.horizDiff(), input_.vertiDiff());
        solver.updateAdvection(0, 0, shear_rep_);

        solver.operatorSplitSolve2DVec(Contrail_, ZERO_BC);
    }

    #ifdef ENABLE_TIMING
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "      Ran contrail tracer transport in " << duration.count() << " ms" << std::endl;
    #endif
}

std::pair<LAGRID::twoDGridVariable,LAGRID::twoDGridVariable> LAGRIDPlumeModel::remapVariable(const VectorUtils::MaskInfo& maskInfo, const BufferInfo& buffers, const Vector_2D& phi, const std::vector<std::vector<int>>& mask) {
    double dx_grid_old = xCoords_[1] - xCoords_[0];
    auto boxGrid = LAGRID::rectToBoxGrid(met_.dy_vec(), dx_grid_old, xEdges_[0], yEdges_[0], phi, mask);

    //Enforce at least x many points in the contrail while limiting minimum/maximum dx and dy
    double dx_grid_new =  std::max(20.0, std::min((maskInfo.maxX - maskInfo.minX) / 50.0, 50.0));
    double dy_grid_new = std::max(5.0, std::min((maskInfo.maxY - maskInfo.minY) / 50.0, 7.0));
    //Need 2 extra points account for the buffer
    int nx_new = floor((maskInfo.maxX - maskInfo.minX) / dx_grid_new) + 2;
    int ny_new = floor((maskInfo.maxY - maskInfo.minY) / dy_grid_new) + 2;
    LAGRID::Remapping remapping(maskInfo.minX - dx_grid_new, maskInfo.minY - dy_grid_new, dx_grid_new, dy_grid_new, nx_new, ny_new);

    auto remappedGrid = LAGRID::mapToStructuredGrid(boxGrid, remapping);
    auto unusedFraction = LAGRID::getUnusedFraction(boxGrid, remapping);

    remappedGrid.addBuffer(buffers.leftBuffer, buffers.rightBuffer, buffers.topBuffer, buffers.botBuffer,0.0);
    unusedFraction.addBuffer(buffers.leftBuffer, buffers.rightBuffer, buffers.topBuffer, buffers.botBuffer,1.0);

    // std::cout << "remap successful" << std::endl;
    return std::make_pair(remappedGrid,unusedFraction);
}

void LAGRIDPlumeModel::remapAllVars(double remapTimestep, const std::vector<std::vector<int>>& mask, const VectorUtils::MaskInfo& maskInfo) {
    //Generate free box grid from met post-advection

    // Step 1: Figure out what the new output grid will be relative to this one
    // Step 2: Generate the weights relating this grid to the next one, using
    //         a simple first-order, conservative remapping
    // Step 3: Apply the weights to each tracer

    // Step 1: Determine the new grid size
    // The below lines look forward and predict, based on the maximum possible
    // fall speed, how far a crystal could fall in one outer time step. Horizontal
    // limits are determined based on maximum shear
    double vertDiffLengthScale = sqrt(VectorUtils::max(diffCoeffY_) * remapTimestep);
    double horizDiffLengthScale = sqrt(VectorUtils::max(diffCoeffX_) * remapTimestep);

    double settlingLengthScale = vFall_[vFall_.size() - 1] * remapTimestep;

    double shear_abs = abs(shear_rep_);

    double shearTop = shear_abs * maskInfo.maxY * remapTimestep; 
    // Bot: Takes into account how much the largest crystals will fall
    double shearBot = shear_abs * -(maskInfo.minY - settlingLengthScale) * remapTimestep;  
    

    // Top shears left if shear > 0, right if shear < 0. Bottom, the opposite.
    double shearLengthScaleLeft = shear_rep_ > 0 ? shearTop : shearBot;
    double shearLengthScaleRight = shear_rep_ > 0 ? shearBot : shearTop; 

    // shearLengthScaleLeft and shearLengthScaleRight can become negative if the contrail sinks a lot.
    // This prevents addBuffer() from throwing an error.
    shearLengthScaleLeft = std::max(shearLengthScaleLeft, 0.0);
    shearLengthScaleRight = std::max(shearLengthScaleRight, 0.0);

    BufferInfo buffers;
    buffers.leftBuffer = (shearLengthScaleLeft + horizDiffLengthScale) * LEFT_BUFFER_SCALING;
    buffers.rightBuffer = (shearLengthScaleRight + horizDiffLengthScale) * RIGHT_BUFFER_SCALING;
    buffers.topBuffer = std::max(vertDiffLengthScale * TOP_BUFFER_SCALING, 100.0);
    buffers.botBuffer = std::min((vertDiffLengthScale + settlingLengthScale) * BOT_BUFFER_SCALING, 300.0);
    //std::cout << buffers.botBuffer << std::endl;
    
    // Step 2: Generate remapping weights
    Vector_1D xEdgesNew, yEdgesNew, xCoordsNew, yCoordsNew;
    auto remapWeights = createRegriddingWeightsSparse(maskInfo, buffers, mask, xEdgesNew, yEdgesNew, xCoordsNew, yCoordsNew);

    // Step 3: Apply them
    int nx_old = xCoords_.size();
    int ny_old = yCoords_.size();
    int nx_new = xCoordsNew.size();
    int ny_new = yCoordsNew.size();

    Vector_3D& pdfRef = iceAerosol_.getPDF();
    Vector_3D volume = iceAerosol_.Volume();

    #pragma omp parallel for default(shared)
    for(UInt n = 0; n < iceAerosol_.getNBin(); n++) {
        //Update pdf and volume
        pdfRef[n] = applyWeights(remapWeights, nx_old,  ny_old,  nx_new,  ny_new, pdfRef[n]);
        volume[n] = applyWeights(remapWeights, nx_old,  ny_old,  nx_new,  ny_new, volume[n]);
    }

    //Only update nx and ny of iceAerosol after the loop, otherwise functions will get messed up if we later add other calls in the loop above
    iceAerosol_.updateNx(pdfRef[0][0].size());
    iceAerosol_.updateNy(pdfRef[0].size());

    //Recalculate VCenters
    iceAerosol_.UpdateCenters(volume, pdfRef);

    //Remap the tracer of contrail presence
    Contrail_ = applyWeights(remapWeights, nx_old,  ny_old,  nx_new,  ny_new, Contrail_);

    //Remap H2O
    H2O_ = applyWeights(remapWeights, nx_old,  ny_old,  nx_new,  ny_new, H2O_);

    //Need to update bottom-of-domain altitude before updating coordinates
    xCoords_ = std::move(xCoordsNew);
    yCoords_ = std::move(yCoordsNew);
    xEdges_ = std::move(xEdgesNew);
    yEdges_ = std::move(yEdgesNew);
    double dx = xEdges_[1] - xEdges_[0];
    double dy = yEdges_[1] - yEdges_[0];
    std::cout << "dx: " << dx << ", dy: " << dy << std::endl;

    // Synchronize the grid desscribed by the met data object
    // and re-interpolate variables such as temperature
    met_.regenerate(yCoords_, yEdges_, xCoords_.size());

    // Convert back from v/v to number density.
    // delta-P/delta-Z is proportional to number density, and accounts
    // for temperature variation (because the calculation of yEdge
    // included it). This uses the hydrostatic assumption:
    //    dp/dz = -rho*g = -(n/V)Mg
    auto pressureEdges = met_.PressEdges(); 
    double localND;
    #pragma omp parallel for
    for (std::size_t j=0; j<yCoords_.size(); j++){
        localND = (pressureEdges[j] - pressureEdges[j+1])/(yEdges_[j+1] - yEdges_[j]);
        for (std::size_t i=0; i<xCoords_.size(); i++){
            Contrail_[j][i] = Contrail_[j][i] * localND;
            H2O_[j][i] = H2O_[j][i] * localND;
            for(UInt n = 0; n < iceAerosol_.getNBin(); n++) {
                pdfRef[n][j][i] = pdfRef[n][j][i] * localND;
            }
        }
    }

    // Set boundary locations to use meteorological H2O
    auto met_H2O = met_.H2O_field();
    Vector_2D uniform = Vector_2D(ny_old,Vector_1D(nx_old));
    for (int j=0; j<ny_old; j++){
        for (int i=0; i<nx_old; i++){
            uniform[j][i] = 1.0;
        }
    }
    Vector_2D usedArea = applyWeights(remapWeights, nx_old,  ny_old,  nx_new,  ny_new, uniform);
    for (int j=0; j<ny_new; j++){
        for (int i=0; i<nx_new; i++){
            if (usedArea[j][i] < 1.0){
                H2O_[j][i] = H2O_[j][i] + std::max(std::min(1.0,(1.0 - usedArea[j][i])),0.0) * met_H2O[j][i];
            }
        }
    }
}

double LAGRIDPlumeModel::totalAirMass() {
    auto numberMask = iceNumberMask();
    auto& mask = numberMask.first;
    double totalAirMass = 0;
    double cellArea = (xEdges_[1] - xEdges_[0]) * (yEdges_[1] - yEdges_[0]);
    double conversion_factor = 1.0e6 * MW_Air / Na; // molec/cm3 * cm3/m3 * kg/mol * mol/molec = kg/m3
    for(std::size_t j = 0; j < yCoords_.size(); j++) {
        for(std::size_t i = 0; i < xCoords_.size(); i++) {
            totalAirMass += mask[j][i] * met_.airMolecDens(j, i) * cellArea * conversion_factor;
        }
    }
    return totalAirMass; // kg/m 
}

// Create an array which maps between 2D grids
Eigen::SparseMatrix<double> LAGRIDPlumeModel::createRegriddingWeightsSparse(const VectorUtils::MaskInfo& maskInfo, const BufferInfo& buffers, const std::vector<std::vector<int>>& mask, Vector_1D& xEdgesNew, Vector_1D& yEdgesNew, Vector_1D& xCoordsNew, Vector_1D& yCoordsNew) {
    Vector_2D phi = Vector_2D(yCoords_.size(), Vector_1D(xCoords_.size()));
    double dx_grid_old = xCoords_[1] - xCoords_[0];
    auto boxGrid = LAGRID::rectToBoxGrid(met_.dy_vec(), dx_grid_old, xEdges_[0], yEdges_[0], phi, mask);

    //Enforce at least x many points in the contrail while limiting minimum/maximum dx and dy
    double dx_grid_new =  std::max(20.0, std::min((maskInfo.maxX - maskInfo.minX) / 50.0, 50.0));
    double dy_grid_new = std::max(5.0, std::min((maskInfo.maxY - maskInfo.minY) / 50.0, 7.0));
    //Need 2 extra points account for the buffer
    int nx_new = floor((maskInfo.maxX - maskInfo.minX) / dx_grid_new) + 2;
    int ny_new = floor((maskInfo.maxY - maskInfo.minY) / dy_grid_new) + 2;
    LAGRID::Remapping remapping(maskInfo.minX - dx_grid_new, maskInfo.minY - dy_grid_new, dx_grid_new, dy_grid_new, nx_new, ny_new);

    auto remappedGrid = LAGRID::mapToStructuredGrid(boxGrid, remapping);
    remappedGrid.addBuffer(buffers.leftBuffer, buffers.rightBuffer, buffers.topBuffer, buffers.botBuffer, 0.0);
    // Contains: phi, xCoords, yCoords

    // Extract the new x and y coordinates
    double dy = remappedGrid.dy;
    double dx = remappedGrid.dx;

    //Update coordinates, accounting for buffers
    yCoordsNew = std::move(remappedGrid.yCoords);
    xCoordsNew = std::move(remappedGrid.xCoords);
    ny_new = yCoordsNew.size();
    nx_new = xCoordsNew.size();
    yEdgesNew = Vector_1D(ny_new + 1);
    xEdgesNew = Vector_1D(nx_new + 1);
    double y0New = yCoordsNew[0];
    double x0New = xCoordsNew[0];
    std::generate(yEdgesNew.begin(), yEdgesNew.end(), [dy, y0New, j = 0.0]() mutable { return y0New + dy*(j++ - 0.5); });
    std::generate(xEdgesNew.begin(), xEdgesNew.end(), [dx, x0New, i = 0.0]() mutable { return x0New + dx*(i++ - 0.5); });

    
    int nx_old = phi[0].size();
    int ny_old = phi.size();

    //Eigen::MatrixXd weights = Eigen::MatrixXd::Zero(nx_old*ny_old,nx_new*ny_new);
    std::vector<Eigen::Triplet<double>> sparseIJV;

    // Create a matrix of areas
    Eigen::MatrixXd areas1 = Eigen::MatrixXd::Zero(ny_new,nx_new);
    for (int iy1=0; iy1<ny_new; iy1++){
        double y_bottom_to = yEdgesNew[iy1];
        double y_top_to = yEdgesNew[iy1+1];
        for (int ix1=0; ix1<nx_new; ix1++){
            double x_left_to = xEdgesNew[ix1];
            double x_right_to = xEdgesNew[ix1+1];
            areas1(iy1,ix1) = (y_top_to - y_bottom_to) * (x_right_to - x_left_to);
        }
    }

    int i0_max=0;
    int i1_max=0;
    // Loop over the source grid
    for (int iy0=0; iy0<ny_old; iy0++){
        for (int ix0=0; ix0<nx_old; ix0++){
            int i0 = ix0 + (iy0 * nx_old);
            if (mask[iy0][ix0] == 0.0){
                continue;
            }
            double x_left   = xEdges_[ix0];
            double x_right  = xEdges_[ix0+1];
            double y_bottom = yEdges_[iy0];
            double y_top    = yEdges_[iy0+1];
            // Always start from the left-most output cell
            int ix1_start = 0;
            int iy1_start = 0; // Revise this later
            for (int iy1=iy1_start; iy1<ny_new; iy1++){
                double y_top_to    = yEdgesNew[iy1+1];
                if (y_top_to <= y_bottom){
                    continue;
                }
                double y_bottom_to = yEdgesNew[iy1];
                if (y_bottom_to >= y_top){
                    break;
                }
                for (int ix1=ix1_start; ix1<nx_new; ix1++){
                    double x_right_to = xEdgesNew[ix1+1];
                    if (x_right_to <= x_left){
                        continue;
                    }
                    double x_left_to = xEdgesNew[ix1];
                    if (x_left_to >= x_right){
                        break;
                    }
                    // 1D index of output cell
                    int i1 = ix1 + (iy1 * nx_new);
                    // How much overlap is there between
                    // the old cell and the new?
                    double x_left_overlap = std::max(x_left,x_left_to);
                    double x_right_overlap = std::min(x_right,x_right_to);
                    double y_bottom_overlap = std::max(y_bottom,y_bottom_to);
                    double y_top_overlap = std::min(y_top,y_top_to);
                    double area_overlap = (y_top_overlap - y_bottom_overlap) * (x_right_overlap - x_left_overlap);
                    //weights(i0,i1) = area_overlap / areas1(iy1,ix1);
                    sparseIJV.push_back(Eigen::Triplet<double>(i0,i1,area_overlap/areas1(iy1,ix1)));
                    i0_max = std::max(i0,i0_max);
                    i1_max = std::max(i1,i1_max);
                }
            }
        }
    }
    Eigen::SparseMatrix<double> weights(nx_old*ny_old,nx_new*ny_new);
    weights.setFromTriplets(sparseIJV.begin(),sparseIJV.end());
    return weights;
}

Vector_2D LAGRIDPlumeModel::applyWeights(const Eigen::SparseMatrix<double>& weights, int nx_old, int ny_old, int nx_new, int ny_new, const Vector_2D& data2D) {
    // Applies weights derived from createRegriddingWeightsSparse and remaps between
    // grids, using a first-order area-conservative approach.
    // Expect oldConc1D and newConc1D to have been previously set up (working arrays only):
    Eigen::VectorXd oldConc1D(ny_old*nx_old);
    Eigen::VectorXd newConc1D(ny_new*nx_new);
    for (int j=0; j<ny_old; j++){
        for (int i=0; i<nx_old; i++){
            oldConc1D(i + j*nx_old) = data2D[j][i];
        }
    }
    newConc1D = oldConc1D.transpose() * weights;
    Vector_2D newConc2D = Vector_2D(ny_new,Vector_1D(nx_new));
    for (int j=0; j<ny_new; j++){
        for (int i=0; i<nx_new; i++){
            newConc2D[j][i] = newConc1D(i + j*nx_new);
        }
    }
    return newConc2D;
}

void LAGRIDPlumeModel::printVector2D(const std::string fieldName, const Vector_2D& dataIn){
    std::cout << " -- BEGIN " << fieldName << " -- " << std::endl;
    int ny = dataIn.size();
    int nx = dataIn[0].size();
    std::cout << ny << " " << nx << std::endl;
    for (int j=0; j<ny; j++) {
        for (int i=0; i<nx; i++) {
            printf("%10.4e ",dataIn[j][i]);
        }
        std::cout << std::endl;
    }
    std::cout << " -- END " << fieldName << " -- " << std::endl;
    return;
}

void LAGRIDPlumeModel::createOutputDirectories() {
    if (!simVars_.TS_AERO ) return;

    std::filesystem::path tsAeroFolderPath(simVars_.TS_FOLDER);
    std::cout << simVars_.TS_FOLDER << std::endl;
    std::cout << "Saving TS_AERO files to: " << simVars_.TS_AERO_FILEPATH << "\n";
    if (std::filesystem::exists(tsAeroFolderPath)) return;

    std::cout << "Creating directory " <<  simVars_.TS_FOLDER << "\n";
    std::filesystem::create_directory(tsAeroFolderPath);
}

void LAGRIDPlumeModel::saveTSAerosol() {
    const double MOD_EPS = 1e-3;
    if ( simVars_.TS_AERO && \
        (( simVars_.TS_AERO_FREQ == 0 ) || \
        ( std::fmod((timestepVars_.curr_Time_s - timestepVars_.timeArray[0])/60.0, simVars_.TS_AERO_FREQ) < MOD_EPS )) ) 
    {
        std::cout << "Saving aerosol data.." << std::endl;    
        int hh = (int) (timestepVars_.curr_Time_s - timestepVars_.timeArray[0])/3600;
        int mm = (int) (timestepVars_.curr_Time_s - timestepVars_.timeArray[0])/60   - 60 * hh;
        int ss = (int) (timestepVars_.curr_Time_s - timestepVars_.timeArray[0])      - 60 * ( mm + 60 * hh );

        Diag::Diag_TS_Phys( simVars_.TS_AERO_FILEPATH.c_str(), hh, mm, ss, \
                        iceAerosol_, H2O_, xCoords_, yCoords_, xEdges_, yEdges_, met_);
        std::cout << "Save Complete" << std::endl;    
    }

}
