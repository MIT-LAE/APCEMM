/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*                      Early Plume Microphysics                    */
/*                              (EPM)                               */
/*                                                                  */
/* Integrate Program File                                           */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 9/27/2018                                 */
/* File                 : Integrate.cpp                             */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <boost/range/algorithm.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp>
#include <boost/numeric/odeint/stepper/controlled_runge_kutta.hpp>
#include <boost/numeric/odeint/iterator/adaptive_iterator.hpp>

#include "AIM/Nucleation.hpp"
#include "KPP/KPP_Parameters.h"
#include "Util/MolarWeights.hpp"
#include "Util/PhysConstant.hpp"
#include "Util/PhysFunction.hpp"
#include "Core/Monitor.hpp"
#include "Core/OutputFilenames.hpp"
#include "Core/Parameters.hpp"
#include "Core/Status.hpp"
#include "EPM/EPM.hpp"
#include "EPM/Models/Original.hpp"
#include "EPM/Models/Original/Indexes.hpp"
#include "EPM/Models/Original/Integrate.hpp"
#include "EPM/Models/Original/RHS.hpp"
#include "EPM/Models/Original/StateObserver.hpp"


using physConst::Na, physConst::R, physConst::PI,
    physConst::RHO_SOOT, physConst::RHO_SULF;

using physFunc::pSat_H2Ol, physFunc::pSat_H2Os;

using namespace EPM::Models::OriginalImpl;

typedef boost::numeric::odeint::runge_kutta_fehlberg78<Vector_1D> error_stepper_type;
typedef boost::numeric::odeint::controlled_runge_kutta<error_stepper_type> controlled_stepper_type;


namespace EPM::Models
{
    std::variant<Output, SimStatus> Original::Integrate(double sootDensity) {
        Output out;
        out.finalTemp = met_.tempRef();
        out.bypassArea = input_.bypassArea();
        out.coreExitTemp = input_.coreExitTemp();

        /* The EPM integrates the plume at flight level, and
         * the vortex parametrisation afterwards accounts for the
         * sinking and its ice crystal loss. */
        SimStatus returnCode = RunMicrophysics(out, sootDensity);

        if (returnCode != SimStatus::EPMSuccess)
            return returnCode;
        return out;
    }

    SimStatus Original::RunMicrophysics(
        Output &state, double sootDensity) {
        double &temperature_K = state.finalTemp;
        double &Ice_rad = state.iceRadius;
        double &Ice_den = state.iceDensity;
        double &Soot_den = state.sootDensity;
        double &H2O_mol = state.H2O_mol;
        double &SO4g_mol = state.SO4g_mol;
        double &SO4l_mol = state.SO4l_mol;
        AIM::Aerosol &SO4Aer = state.SO4Aer;
        AIM::Aerosol &IceAer = state.IceAer;
        double &Area = state.area;
        double Ab0 = state.bypassArea;
        double Tc0 = state.coreExitTemp;

        /* Persistence check: ambient RHi at flight level. rhiRef is in %. */
        const double relHumidity_i_Final = met_.rhiRef() / 100.0;

      /* Homogeneous NAT, using supercooling requirement *
       * Heterogeneous NAT, no supercooling requirement */

      /* Number of SO4 size distribution bins based on specified min and max
       * radii * and volume ratio between two adjacent bins */
      const UInt SO4_NBIN =
          static_cast<UInt>(
              std::floor(1 + 3.0*log(LA_R_HIG / LA_R_LOW) / log(LA_VRAT)));

      if (LA_VRAT <= 1.0) {
        std::cout << "\nVolume ratio of consecutive bins for SO4 has to be "
                     "greater than 1.0 ( LA_VRAT = "
                  << LA_VRAT << " )";
      }

      /* SO4 bin radius and volume centers */
      Vector_1D SO4_rJ( SO4_NBIN    , 0.0 ); // bin centers, radius
      Vector_1D SO4_rE( SO4_NBIN + 1, 0.0 ); // bin edges, radius
      Vector_1D SO4_vJ( SO4_NBIN    , 0.0 ); // bin centers, volume

      /* Adjacent bin radius ratio */
      const double LA_RRAT = cbrt(LA_VRAT);

        /* Initialize bin center and edge radii, as well as volume */
        SO4_rE[0] = LA_R_LOW;
        for ( UInt iBin = 1; iBin < SO4_NBIN + 1; iBin++ )
            SO4_rE[iBin] = SO4_rE[iBin-1] * LA_RRAT;                                                           /* [m] */

        for ( UInt iBin = 0; iBin < SO4_NBIN; iBin++ ) {
            SO4_rJ[iBin] = 0.5 * ( SO4_rE[iBin] + SO4_rE[iBin+1] );                                            /* [m] */
            SO4_vJ[iBin] = 4.0 / double(3.0) * PI * SO4_rJ[iBin] * SO4_rJ[iBin] * SO4_rJ[iBin]; /* [m^3] */
        }

        /* Number of ice size distribution bins based on specified min and max radii *
        * and volume ratio between two adjacent bins */
        const double ice_vrat = optInput_.ADV_AERO_ICE_BIN_VRAT;
        const double ice_r_low = optInput_.ADV_AERO_ICE_BIN_R_LOW;
        const double ice_r_hig = optInput_.ADV_AERO_ICE_BIN_R_HIG;
        const UInt Ice_NBIN = static_cast<UInt>(
            std::floor(1 + 3.0*log(ice_r_hig / ice_r_low) / log(ice_vrat)));

        /* Adjacent bin radius ratio */
        const double PA_RRAT = cbrt( ice_vrat );

        /* Ice bin center and edge radii */
        Vector_1D Ice_rJ( Ice_NBIN    , 0.0 );
        Vector_1D Ice_rE( Ice_NBIN + 1, 0.0 );
        Ice_rE[0] = ice_r_low;
        for ( UInt iBin = 1; iBin < Ice_NBIN + 1; iBin++ )
            Ice_rE[iBin] = Ice_rE[iBin-1] * PA_RRAT;                                                           /* [m] */

        for ( UInt iBin = 0; iBin < Ice_NBIN; iBin++ )
            Ice_rJ[iBin] = 0.5 * ( Ice_rE[iBin] + Ice_rE[iBin+1] );                                            /* [m] */


        /* Create coagulation kernels */
        const AIM::Coagulation Kernel( AIM::CoagulationPair::LiquidLiquid, SO4_rJ, SO4_vJ, RHO_SULF, temperature_K, simVars_.pressure_Pa );
        const AIM::Coagulation Kernel_SO4_Soot( AIM::CoagulationPair::LiquidSoot, SO4_rJ, RHO_SULF, EI_.getSootRad(), RHO_SOOT, temperature_K, simVars_.pressure_Pa );

        const Vector_1D KernelSO4Soot = Kernel_SO4_Soot.getKernel_1D();

        /* Create SO4 aerosol number distribution.
         * We allocate the PDF with a very small number of existing particles */

        double nSO4 = 1.00E-10; /* [#/cm^3] */
        double rSO4 = 5.00E-10; /* [m] */
        double sSO4 = 1.4; /* For lognormal distributions, sSO4 = sigma */
//        double sSO4 = rSO4/10; /* For normal distributions, sSO4 */

        /* Type of sulfate distribution */
        /* lognormal distribution */

        AIM::Aerosol nPDF_SO4( SO4_rJ, SO4_rE, nSO4, rSO4, sSO4, "lognormal" );

        double n_air_amb = Na * simVars_.pressure_Pa / (R * temperature_K * 1.0e6) ;
        double n_air_eng = Na * simVars_.pressure_Pa / (R * Tc0 * 1.0e6) ;


        /* Store ambient mixing ratios */
        double H2O_amb  = VAR_[ind_H2O]/n_air_amb;     /* [mol/mol] */
        double SO4_amb  = VAR_[ind_SO4]/n_air_amb;     /* [mol/mol] */
        double SO4l_amb, SO4g_amb;                     /* [mol/mol] */
        double HNO3_amb = VAR_[ind_HNO3]/n_air_amb;    /* [mol/mol] */
        double Soot_amb = sootDensity/n_air_amb;       /* [#/molec_air] */

        /* Add emissions of one engine using ambient mixing ratio plus
         * emitted exhaust divided by engine-exit air number density n_air_eng.
         * fuel_conc has units [kg_fuel / cm^3] at the engine exit plane. */
        const double fuel_conc = aircraft_.FuelFlow() / double(aircraft_.EngNumber()) / aircraft_.VFlight() / Ab0 * 1.00E-06;

        const double H2O_emis = EI_.getH2O() / ( MW_H2O * 1.0E+03 ) * Na * fuel_conc; /* [molec/cm^3] */
        VAR_[ind_H2O] = H2O_amb + H2O_emis / n_air_eng; /* [mol/mol] */

        /* Variable SO2 */
        const double SO4_emis = EI_.getSO2toSO4() * EI_.getSO2() / ( MW_H2SO4 * 1.0E+03 ) * Na * fuel_conc; /* [molec/cm^3] */
        VAR_[ind_SO4] = SO4_amb + SO4_emis / n_air_eng; /* [mol/mol] */

        const double soot_emis = EI_.getSoot() / ( 4.0 / double(3.0) * PI * RHO_SOOT * 1.00E+03 * EI_.getSootRad() * EI_.getSootRad() * EI_.getSootRad() ) * fuel_conc; /* [#/cm^3] */
        const double varSoot = Soot_amb + soot_emis / n_air_eng; /* [#/molec_air] */

        /* Unit check:
         *   fuel_conc = [kg_f/s] / [m/s] / [m^2] * [m^3/cm^3] = [kg_f/cm^3]
         *   soot_emis = [g/kg_f] / ([kg/m^3] * [g/kg] * [m^3]) * [kg_f/cm^3] = [#/cm^3]
         *   varSoot   = [#/molec_air] + [#/cm^3] / [molec_air/cm^3] = [#/molec_air] */

        /* Compute ambient share of liquid, gaseous sulfates
         * assuming that SO4 is in phase equilibrium */

        SO4g_amb = SO4_amb;
        SO4l_amb = 0.0E+00;

        SO4g_amb = ( SO4g_amb > 0.0 ) ? SO4g_amb : 0.0;
        SO4l_amb = ( SO4l_amb > 0.0 ) ? SO4l_amb : 0.0;

        /* Outer time grid: 301 log-spaced points from 1e-4s to 180s
        which is the step at which we start the mature plume model */
        const double t_3mins = 3.0 * 60.0;
        const Vector_1D timeArray = makeEPMTimeGrid(301, 1.0E-04, t_3mins, t_3mins);
        const UInt nTime = timeArray.size();
        UInt iTime = 0;

        // This is the last index of timeArray when tEnd = tInsert = 180s
        // Kept here a search in case we decide to run the EPM for longer (tEnd > tInsert in makeEPMTimeGrid)
        // and keep tInsert = 3 mins as the mature plume model starting point.
        const UInt iTime_3mins = std::lower_bound(timeArray.begin(), timeArray.end(), t_3mins) - timeArray.begin();

        /* State vector x:
         * x[EPM_ind_Trac] : Plume dilution tracer,          Unit [-]
         * x[EPM_ind_T]    : Temperature,                    Unit K
         * x[EPM_ind_P]    : Pressure,                       Unit Pa
         * x[EPM_ind_H2O]  : Water vapor mixing ratio,       Unit mol/mol
         * x[EPM_ind_SO4]  : Total SO4 mixing ratio,         Unit mol/mol
         * x[EPM_ind_SO4l] : Liquid SO4 mixing ratio,        Unit mol/mol
         * x[EPM_ind_SO4g] : Gaseous SO4 mixing ratio,       Unit mol/mol
         * x[EPM_ind_SO4s] : On-soot SO4 mixing ratio,       Unit mol/mol
         * x[EPM_ind_HNO3] : Gaseous HNO3 mixing ratio,      Unit mol/mol
         * x[EPM_ind_Part] : Particle mixing ratio,          Unit #/molec_air
         * x[EPM_ind_ParR] : Particle radius,                Unit m
         * x[EPM_ind_the1] : Soot coverage (adsorption),     Unit [-]
         * x[EPM_ind_the2] : Soot coverage (scavenging),     Unit [-]
         */

        const std::vector<UInt> EPM_ind = {
            EPM_ind_Trac, EPM_ind_T, EPM_ind_P, EPM_ind_H2O, EPM_ind_SO4,
            EPM_ind_SO4l, EPM_ind_SO4g, EPM_ind_SO4s, EPM_ind_HNO3,
            EPM_ind_Part, EPM_ind_ParR, EPM_ind_the1, EPM_ind_the2};

        Vector_1D x(13, 0.0);

        /* Initial conditions */
        x[EPM_ind_Trac] = 1.0;
        x[EPM_ind_T]    = Tc0;
        x[EPM_ind_P]    = simVars_.pressure_Pa;
        x[EPM_ind_H2O]  = VAR_[ind_H2O];
        x[EPM_ind_SO4]  = VAR_[ind_SO4];

        /* Assume that emitted SO4 is purely gaseous as temperature is high */
        x[EPM_ind_SO4l] = 0.0;
        x[EPM_ind_SO4g] = VAR_[ind_SO4];
        x[EPM_ind_SO4s] = 0.0;
        x[EPM_ind_HNO3] = HNO3_amb;
        x[EPM_ind_Part] = varSoot;
        x[EPM_ind_ParR] = EI_.getSootRad();
        x[11] = 0.0;
        x[12] = 0.0;

        StateHistory history(EPM_ind, simVars_.TS_FOLDER + OutputFiles::MICROPHYSICS);

        /* Set at any adaptive step where the plume reaches water saturation */
        bool reachedWaterSat = false;
        const WaterSaturationObserver waterSatObserver(reachedWaterSat);

        /* Adaptive stepper that can be shared by all integration steps */
        const auto stepper = boost::numeric::odeint::make_controlled< error_stepper_type >( EPM_ATOLS, EPM_RTOLS );

        /* Creating ode's right hand side */
        gas_aerosol_rhs rhs(
            temperature_K, simVars_.pressure_Pa,
            H2O_amb, SO4_amb, SO4l_amb, SO4g_amb, HNO3_amb, Soot_amb,
            EI_.getSootRad(), KernelSO4Soot, nPDF_SO4);

        /* Current time step in s */
        double currTimeStep;

        /* Freshly nucleated aerosol characteristics */
        double SO4l_prev, n_air_prev, n_air;
        double nPDF_new, x_star, nTot, radSO4;

        /* Get conditions @ 3 mins */
        double PartRad_3mins = 0;
        double PartDens_3mins = 0;
        double H2OMol_3mins = 0;
        double SO4l_3mins = 0;
        double SO4g_3mins = 0;
        double Tracer_3mins = 0;
        double n_air_3mins = 0;
        AIM::Aerosol pSO4pdf_3mins( nPDF_SO4 );

        // Record the initial state 
        history.record( x, timeArray[0] );

        while ( iTime + 1 < nTime ) {

            /* Update current time step */
            currTimeStep = timeArray[iTime+1] - timeArray[iTime];

            /* Coagulation of sulfate aerosols */
            nPDF_SO4.Coagulate( currTimeStep, Kernel );

            /* Liquid sulfate budget at the start of the interval, timeArray[iTime] */
            SO4l_prev  = x[EPM_ind_SO4l];
            n_air_prev = Na * x[EPM_ind_P] / ( R * x[EPM_ind_T] * 1.0e6 );

            /* Diffusion + Water uptake. On return, x is the state at the end
             * of the interval, timeArray[iTime+1]. */
            boost::numeric::odeint::integrate_adaptive(
                stepper, rhs, x, timeArray[iTime], timeArray[iTime+1], currTimeStep/100.0,
                waterSatObserver );
            history.record( x, timeArray[iTime+1] );

            /* Liquid sulfate budget at the end of the interval */
            n_air = Na * x[EPM_ind_P] / ( R * x[EPM_ind_T] * 1.0e6 );

            nPDF_new = x[EPM_ind_SO4l] * n_air - SO4l_prev * n_air_prev;

            if ( nPDF_new >= 1.0E-20 ) {
                double RH_liquid = n_air * x[EPM_ind_H2O] * 1.0e6 * physConst::kB * x[EPM_ind_T] / pSat_H2Ol( x[EPM_ind_T] ) ; 
                x_star   = AIM::x_star( x[EPM_ind_T], RH_liquid, std::max(x[EPM_ind_SO4g] * n_air, 0.0) );
                nTot     = AIM::nTot( x[EPM_ind_T], x_star, RH_liquid, std::max(x[EPM_ind_SO4g]*n_air, 0.0) );
                nTot     = ( nTot <= 1.0E-20 ) ? 1.0E-20 : nTot;
                radSO4   = AIM::radCluster( x_star, nTot );

                if ( radSO4 >= 1.0E-10 ) {
                    AIM::Aerosol nPDF_SO4_new( SO4_rJ, SO4_rE, nPDF_new, radSO4, sSO4, "lognormal" );

                    nPDF_SO4.addAerosolToPDF(nPDF_SO4_new);
                }
            }

            iTime++;

            /* Snapshot @ 3mins: x, n_air and the sulfate PDF are all at timeArray[iTime] */
            if ( iTime == iTime_3mins ) {
                PartRad_3mins  = x[EPM_ind_ParR];
                PartDens_3mins = x[EPM_ind_Part] * n_air;
                n_air_3mins    = n_air;
                H2OMol_3mins   = x[EPM_ind_H2O];
                Tracer_3mins   = x[EPM_ind_Trac];
                SO4l_3mins     = x[EPM_ind_SO4l];
                SO4g_3mins     = x[EPM_ind_SO4g];
                pSO4pdf_3mins.updatePdf( nPDF_SO4.getPDF() );
            }

        }

#pragma omp critical
        {
            history.print2File();
        }
        /* Output variables */
        /* Check if contrail is water supersaturated at some point during formation */
        if ( !simVars_.CHEMISTRY && !reachedWaterSat ) {
            std::cout << "EndSim: Never reaches water saturation... ending simulation" << std::endl;
            return SimStatus::NoWaterSaturation;
        }

        /* Persistent contrail */
        if ( relHumidity_i_Final >= 1.0 ) {

             Ice_rad  = PartRad_3mins;
             Ice_den  = PartDens_3mins;
             Soot_den = 0.0;

        }
        /* No persistent contrail */
        else {
             Ice_rad  = EI_.getSootRad();
             Ice_den  = 0.0;
             Soot_den = PartDens_3mins;
             std::cout << "No persistent contrail..." << std::endl;
             if (!simVars_.CHEMISTRY) return SimStatus::NoPersistence;
        }
	    std::cout << "Ice_den=" << Ice_den << std::endl;

        H2O_mol  = H2OMol_3mins;
        SO4l_mol = SO4l_3mins;
        SO4g_mol = SO4g_3mins;
        SO4Aer = pSO4pdf_3mins;

        const double expsIce = 1.15;
        AIM::Aerosol solidAer( Ice_rJ, Ice_rE, Ice_den, std::max(Ice_rad * exp ( -2.5 * log(expsIce) * log(expsIce) ), 1.5 * optInput_.ADV_AERO_ICE_BIN_R_LOW ), expsIce, "lognormal" );
        IceAer = solidAer;

        /* Compute plume area */
        Area = Ab0 * n_air_eng / (n_air_3mins * Tracer_3mins);

//        std::cout << " Soot end: " << ( varSoot - Soot_amb  ) * 1E6 * Ab0 << " [#/m]\n";
//        std::cout << " Soot end: " << ( Soot_den - Soot_amb ) * 1E6 * Ab0 / Tracer_3mins << "[#/m]\n";

        return SimStatus::EPMSuccess;
    } /* End of RunMicrophysics */

}
