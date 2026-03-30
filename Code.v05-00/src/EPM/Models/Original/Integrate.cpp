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
#include "Core/Parameters.hpp"
#include "Core/Status.hpp"
#include "EPM/EPM.hpp"
#include "EPM/Models/Original.hpp"
#include "EPM/Models/Original/Indexes.hpp"
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

        /* Get mean vortex displacement in [m] */
        double delta_z = VORTEX_SINKING ? aircraft_.vortex().z_center() : 0.0;

        /* Compute adiabatic and real temperature changes */
        double delta_T    = -optInput_.ADV_AMBIENT_LAPSERATE * delta_z * 1.0E-03;
        /*          [ K ] =          [ K/km ] *  [ m ]  * [ km/m ]
         * The minus sign is because delta_z is the distance pointing down */

        SimStatus returnCode = RunMicrophysics(out, sootDensity, delta_T);

        if (returnCode != SimStatus::EPMSuccess)
            return returnCode;
        return out;
    }

    SimStatus Original::RunMicrophysics(
        Output &state, double sootDensity, double delta_T) {
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

        double relHumidity_i_Final;

      /* RH is in % */
      relHumidity_i_Final = met_.rhwRef() * pSat_H2Ol(temperature_K + delta_T) /
                            pSat_H2Os(temperature_K + delta_T) / 100.0;

      double final_Temp = temperature_K + delta_T;

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
        const UInt Ice_NBIN = static_cast<UInt>(
            std::floor(1 + 3.0*log(PA_R_HIG / PA_R_LOW) / log(PA_VRAT)));

        /* Adjacent bin radius ratio */
        const double PA_RRAT = cbrt( PA_VRAT );

        /* Ice bin center and edge radii */
        Vector_1D Ice_rJ( Ice_NBIN    , 0.0 );
        Vector_1D Ice_rE( Ice_NBIN + 1, 0.0 );
        Ice_rE[0] = PA_R_LOW;
        for ( UInt iBin = 1; iBin < Ice_NBIN + 1; iBin++ )
            Ice_rE[iBin] = Ice_rE[iBin-1] * PA_RRAT;                                                           /* [m] */

        for ( UInt iBin = 0; iBin < Ice_NBIN; iBin++ )
            Ice_rJ[iBin] = 0.5 * ( Ice_rE[iBin] + Ice_rE[iBin+1] );                                            /* [m] */


        /* Create coagulation kernels */
        const AIM::Coagulation Kernel( "liquid", SO4_rJ, SO4_vJ, RHO_SULF, temperature_K, simVars_.pressure_Pa );
        const AIM::Coagulation Kernel_SO4_Soot( "liquid", SO4_rJ, RHO_SULF, EI_.getSootRad(), RHO_SOOT, temperature_K, simVars_.pressure_Pa );

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


        /* Store ambient concentrations */
        double H2O_amb  = VAR_[ind_H2O]/n_air_amb;     /* [molec/cm^3] */
        double SO4_amb  = VAR_[ind_SO4]/n_air_amb;     /* [molec/cm^3] */
        double SO4l_amb, SO4g_amb;               /* [molec/cm^3] */
        double HNO3_amb = VAR_[ind_HNO3]/n_air_amb;    /* [molec/cm^3] */
        double Soot_amb = sootDensity/n_air_amb; ///n_air_amb; /* [#/cm^3] */

        /* Add emissions of one engine to concentration arrays */
        VAR_[ind_H2O] += EI_.getH2O() / ( MW_H2O  * 1.0E+03 ) * aircraft_.FuelFlow() / double(aircraft_.EngNumber()) / aircraft_.VFlight() * Na / Ab0     * 1.00E-06 ;
        /* [ molec/cm^3 ] += [ g/kgf ]   / [ kg/mol ] * [ g/kg ] * [ kgf/s ]                                  / [ m/s ]      * [ molec/mol ] / [ m^2 ] * [ m^3/cm^3 ]
         *                += [ molec/cm^3 ] */

        VAR_[ind_H2O] = VAR_[ind_H2O] / n_air_eng ;

        /* Variable SO2 */
        VAR_[ind_SO4] += EI_.getSO2toSO4() * EI_.getSO2() / ( MW_H2SO4  * 1.0E+03 ) * aircraft_.FuelFlow() / double(aircraft_.EngNumber()) / aircraft_.VFlight() * Na / Ab0 * 1.00E-06;
        VAR_[ind_SO4] = VAR_[ind_SO4] / n_air_eng;

        double varSoot = Soot_amb + EI_.getSoot() / ( 4.0 / double(3.0) * PI * RHO_SOOT * 1.00E+03 * EI_.getSootRad() * EI_.getSootRad() * EI_.getSootRad() ) * aircraft_.FuelFlow() / double(aircraft_.EngNumber()) / aircraft_.VFlight() / Ab0 * 1.00E-06 ; /* [ #/cm^3 ] */
        varSoot = varSoot / n_air_eng;

        /* Unit check:
         *                 = [#/cm^3] + [g/kg_f]     / (                                         [kg/m^3]            * [g/kg]   * [m^3]                                               ) * [kg_f/s]                         / [m/s]            / [m^2] * [m^3/cm^3]
         *                 = [#/cm^3] */

        /* Compute ambient share of liquid, gaseous sulfates
         * assuming that SO4 is in phase equilibrium */

        SO4g_amb = SO4_amb;
        SO4l_amb = 0.0E+00;

        SO4g_amb = ( SO4g_amb > 0.0 ) ? SO4g_amb : 0.0;
        SO4l_amb = ( SO4l_amb > 0.0 ) ? SO4l_amb : 0.0;

        /* Adaptive time stepping? */
        const bool adaptiveStep = 1;

        UInt nTime = 301;
        UInt iTime;
        Vector_1D timeArray(nTime);
        double timeInitial = 1.0E-04;
        double log10_timeInitial = log10(timeInitial);
        double timeFinal   = 2.0E+03;
        double log10_timeFinal = log10(timeFinal);
        double denominator = 1/double( nTime - 1.0 );
        for ( iTime = 0; iTime < nTime; iTime++ ) {
            timeArray[iTime] = pow( 10.0, log10_timeInitial + iTime * ( log10_timeFinal - log10_timeInitial ) * denominator );
        }

        UInt iTime_3mins;
        auto const it = std::lower_bound(timeArray.begin(), timeArray.end(), 3.0 * 60.0);
        iTime_3mins = it - timeArray.begin();

        /* vars[0] : Temperature,                   Unit K
         * vars[1] : Water molecular concentration, Unit molecules/cm^3
         * vars[2] : Gaseous SO4 concentration,     Unit molecules/cm^3
         * vars[3] : Liquid SO4 concentration,      Unit molecules/cm^3
         * vars[4] : Gaseous HNO3 concentration,    Unit molecules/cm^3
         * vars[5] : Liquid HNO3 concentration,     Unit molecules/cm^3
         * vars[6] : Soot ambient concentration,    Unit #/cm^3
         * vars[7] : Ice particle radius,           Unit m
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
        x[EPM_ind_HNO3] = VAR_[ind_HNO3] / n_air_eng;
        x[EPM_ind_Part] = varSoot;
        x[EPM_ind_ParR] = EI_.getSootRad();
        x[11] = 0.0;
        x[12] = 0.0;

        Vector_2D obs_Var;
        Vector_1D obs_Time;

        StateObserver observer(obs_Var, obs_Time, EPM_ind, input_.fileName_micro(), 2);

        /* Creating ode's right hand side */
        gas_aerosol_rhs rhs(
            temperature_K, simVars_.pressure_Pa, delta_T,
            H2O_amb, SO4_amb, SO4l_amb, SO4g_amb, HNO3_amb, Soot_amb,
            EI_.getSootRad(), KernelSO4Soot, nPDF_SO4);

        /* Total number of steps */
        UInt totSteps = 0;

        /* Current time step in s */
        double currTimeStep;

        /* Dilution factor */
        // double dilFactor, dilFactor_b;

        /* Freshly nucleated aerosol characteristics */
        double SO4l, SO4l_b, T_b, P_b, n_air_prev, n_air;
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

        iTime = 0;
        while ( iTime + 1 < nTime ) {

            /* Update current time step */
            currTimeStep = timeArray[iTime+1] - timeArray[iTime];

            /* Coagulation of sulfate aerosols */
            nPDF_SO4.Coagulate( currTimeStep, Kernel );

            /* Get dilution factor */
            if ( iTime == 0 ) {
                // dilFactor_b = 1.0;
                SO4l_b = 0.0;
                T_b = Tc0; // Set initial temperature in plume to core exit temp
                P_b = simVars_.pressure_Pa; // Set initial pressure in plume to ambient pressure
            }
            else {
                // dilFactor_b = observer.m_states[observer.m_states.size()-1][EPM_ind_Trac];
                SO4l_b = observer.m_states[observer.m_states.size()-1][EPM_ind_SO4l];
                T_b = observer.m_states[observer.m_states.size()-1][EPM_ind_T];
                P_b = observer.m_states[observer.m_states.size()-1][EPM_ind_P];
            }

            /* Diffusion + Water uptake */
            if ( adaptiveStep == 1 ) {
                totSteps += boost::numeric::odeint::integrate_adaptive( boost::numeric::odeint::make_controlled< error_stepper_type >( EPM_ATOLS, EPM_RTOLS ), rhs, x, timeArray[iTime], timeArray[iTime+1], currTimeStep/100.0, observer );
            }
            else {
                totSteps += boost::numeric::odeint::integrate( rhs, x, timeArray[iTime], timeArray[iTime+1], currTimeStep/100.0, observer );
            }

            // dilFactor = observer.m_states[observer.m_states.size()-1][EPM_ind_Trac] / dilFactor_b;
            SO4l = observer.m_states[observer.m_states.size()-1][EPM_ind_SO4l] ;

            n_air = Na * x[EPM_ind_P]/(R * x[EPM_ind_T] * 1.0e6);
            n_air_prev = Na *  P_b/(R * T_b * 1.0e6);  // P_b and T_b seems to be used uninitialized during the first timestep with iTime = 0

            nPDF_new = ( SO4l*n_air - SO4l_b*n_air_prev);

            if ( nPDF_new >= 1.0E-20 ) {
                x_star   = AIM::x_star( observer.m_states[observer.m_states.size()-1][EPM_ind_T], observer.m_states[observer.m_states.size()-1][EPM_ind_H2O] * n_air, std::max(x[EPM_ind_SO4g] * n_air, 0.0) );
                nTot     = AIM::nTot( observer.m_states[observer.m_states.size()-1][EPM_ind_T], x_star, observer.m_states[observer.m_states.size()-1][EPM_ind_H2O] * n_air, std::max(x[EPM_ind_SO4g]*n_air, 0.0) );
                nTot     = ( nTot <= 1.0E-20 ) ? 1.0E-20 : nTot;
                radSO4   = AIM::radCluster( x_star, nTot );
//                rho_Sulf = AIM::rho( x_star, observer.m_states[observer.m_states.size()-1][EPM_ind_T]);

                if ( radSO4 >= 1.0E-10 ) {
                    AIM::Aerosol nPDF_SO4_new( SO4_rJ, SO4_rE, nPDF_new, radSO4, sSO4, "lognormal" );

                    nPDF_SO4.addAerosolToPDF(nPDF_SO4_new);
                }
            }

            /* For debugging purposes */
            if ( 0 ) {
                std::cout << "Time: " << timeArray[iTime+1] << ", Number: " << nPDF_SO4.Moment() << ", Radius: " << nPDF_SO4.Radius() << "\n";
            }

            totSteps += 1;
            iTime++;

            /* Aerosol PDF @ 3mins */
            if ( iTime == iTime_3mins ) {
                PartRad_3mins  = observer.m_states[observer.m_states.size()-1][EPM_ind_ParR];
                PartDens_3mins = observer.m_states[observer.m_states.size()-1][EPM_ind_Part] * n_air;
                n_air_3mins = n_air;
                H2OMol_3mins   = observer.m_states[observer.m_states.size()-1][EPM_ind_H2O]; //* observer.m_states[observer.m_states.size()-1][EPM_ind_P] / (kB * observer.m_states[observer.m_states.size()-1][EPM_ind_T]) * 1.0E-06;
                Tracer_3mins   = observer.m_states[observer.m_states.size()-1][EPM_ind_Trac];
//                SO4pdf_3mins   = nPDF_SO4;
                SO4l_3mins     = observer.m_states[observer.m_states.size()-1][EPM_ind_SO4l]; // * observer.m_states[observer.m_states.size()-1][EPM_ind_P] / ( kB * observer.m_states[observer.m_states.size()-1][EPM_ind_T] * 1.0E+06 ) ;
                SO4g_3mins     = observer.m_states[observer.m_states.size()-1][EPM_ind_SO4g]; // * observer.m_states[observer.m_states.size()-1][EPM_ind_P] / ( kB * observer.m_states[observer.m_states.size()-1][EPM_ind_T] * 1.0E+06 ) ;
//                pSO4pdf_3mins  = new AIM::Aerosol( nPDF_SO4 );
                pSO4pdf_3mins.updatePdf( nPDF_SO4.getPDF() );

            }

        }

#pragma omp critical
        {
            observer.print2File();
        }
        /* Output variables */
        /* Check if contrail is water supersaturated at some point during formation */
        if ( !simVars_.CHEMISTRY && !observer.checkwatersat() ) {
            std::cout << "EndSim: Never reaches water saturation... ending simulation" << std::endl;
            //exit(0);
            return SimStatus::NoWaterSaturation;
        }

        /* Persistent contrail */
        if ( relHumidity_i_Final >= 1.0 ) {

             Ice_rad  = PartRad_3mins;
             Ice_den  = PartDens_3mins;
             Soot_den = 0.0;
             H2O_mol  = pSat_H2Os( final_Temp ); // / ( kB * final_Temp * 1.0E+06 )  ;

        }
        /* No persistent contrail */
        else {
             Ice_rad  = EI_.getSootRad();
             Ice_den  = 0.0;
             Soot_den = PartDens_3mins;
             H2O_mol  = H2OMol_3mins;
             std::cout << "No persistent contrail..." << std::endl;
             if (!simVars_.CHEMISTRY) return SimStatus::NoPersistence;
        }
	    std::cout << "Ice_den=" << Ice_den << std::endl;

        SO4l_mol = SO4l_3mins;
        SO4g_mol = SO4g_3mins;
        SO4Aer = pSO4pdf_3mins;

        const double expsIce = 1.15;
        AIM::Aerosol solidAer( Ice_rJ, Ice_rE, Ice_den, std::max(Ice_rad * exp ( -2.5 * log(expsIce) * log(expsIce) ), 1.5 * PA_R_LOW ), expsIce, "lognormal" );
        IceAer = solidAer;

        /* Compute plume area */
        Area = Ab0 * n_air_eng / (n_air_3mins * Tracer_3mins);

//        std::cout << " Soot end: " << ( varSoot - Soot_amb  ) * 1E6 * Ab0 << " [#/m]\n";
//        std::cout << " Soot end: " << ( Soot_den - Soot_amb ) * 1E6 * Ab0 / Tracer_3mins << "[#/m]\n";

        /* Update temperature after vortex sinking */
        temperature_K = final_Temp;

        return SimStatus::EPMSuccess;
    } /* End of RunMicrophysics */

}
