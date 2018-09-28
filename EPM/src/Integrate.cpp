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

#include "Integrate.hpp"

namespace EPM
{
    int Integrate( double temperature_K, double pressure_Pa, double relHumidity_w, double varArray[], double fixArray[], const Aircraft &ac, const Emission &EI )
    {

        /* Get mean vortex displacement in [m] */
        double delta_z = ac.getVortexdeltaz1(); 

        /* Compute adiabatic and real temperature changes */
        double delta_T_ad = - physConst::GAMMA_AD * delta_z * 1.0E-03;
        double delta_T    = -            GAMMA    * delta_z * 1.0E-03;
        /*        [ K ]   = [ K/km ] *  [ m ]  * [ km/m ] 
         * The minus sign is because delta_z is the distance pointing down */
        std::cout << "Temperature lapse rate: " << GAMMA << " [K/km] \n";
        std::cout << "Vortex sinking: " << delta_z << " [m] \n";
        std::cout << "Temperature increase: " << delta_T << " [K] \n";

        RunMicrophysics( temperature_K, pressure_Pa, relHumidity_w, varArray, fixArray, ac, EI, delta_T_ad, delta_T );

        return EPM_SUCCESS;

    } /* End of Integrate */

    int RunMicrophysics( double temperature_K, double pressure_Pa, double relHumidity_w, double varArray[], double fixArray[], const Aircraft &ac, const Emission &EI, double delta_T_ad, double delta_T )
    {
    
        double relHumidity_i_Amb, relHumidity_i_postVortex, relHumidity_i_Final;
        /* relHumidity_w is in % */
        relHumidity_i_Amb        = relHumidity_w * physFunc::pSat_H2Ol( temperature_K ) \
                                                 / physFunc::pSat_H2Os( temperature_K ) / 100.0;
        relHumidity_i_postVortex = relHumidity_w * physFunc::pSat_H2Ol( temperature_K + delta_T_ad ) \
                                                 / physFunc::pSat_H2Os( temperature_K + delta_T_ad ) / 100.0;
        relHumidity_i_Final      = relHumidity_w * physFunc::pSat_H2Ol( temperature_K + delta_T ) \
                                                 / physFunc::pSat_H2Os( temperature_K + delta_T ) / 100.0;

        std::cout << "\nInitial ambient relative humidity w.r.t ice: " << 100.0 * relHumidity_i_Amb << " %" ;
        std::cout << "\nPost vortex relative humidity w.r.t ice    : " << 100.0 * relHumidity_i_postVortex << " %";
        std::cout << "\nFinal relative humidity w.r.t ice          : " << 100.0 * relHumidity_i_Final << " %";

        double partPHNO3_Hom, partPHNO3_Het, satHNO3_Hom, satHNO3_Het;
        double final_Temp = temperature_K + delta_T;
        double offset_Temp = final_Temp + T_NAT_SUPERCOOL;
        /* airDens = pressure_Pa / ( kB * temperature_K ) * 1.00E-06
         * partPHNO3 = [HNO3] / airDens * pressure_Pa = [HNO3] * kB * temperature_K * 1.00E+06 */

        /* Homogeneous NAT, using supercooling requirement *
         * Heterogeneous NAT, no supercooling requirement */
        partPHNO3_Hom = varArray[ind_HNO3] * physConst::kB * offset_Temp * 1.00E+06;
        partPHNO3_Het = varArray[ind_HNO3] * physConst::kB * final_Temp  * 1.00E+06;
        if ( relHumidity_i_Amb > 1.0 ) {
            /* partPH2O = pSat_H2Os( temperature_K ) */
            satHNO3_Hom = partPHNO3_Hom / physFunc::pSat_HNO3( offset_Temp, physFunc::pSat_H2Os( final_Temp ));
            satHNO3_Het = partPHNO3_Het / physFunc::pSat_HNO3( final_Temp , physFunc::pSat_H2Os( final_Temp ));
        } else {
            satHNO3_Hom = partPHNO3_Hom / physFunc::pSat_HNO3( offset_Temp, relHumidity_i_Final * physFunc::pSat_H2Os( final_Temp ) );
            satHNO3_Het = partPHNO3_Het / physFunc::pSat_HNO3( final_Temp , relHumidity_i_Final * physFunc::pSat_H2Os( final_Temp ) );
        }
        std::cout << "\nHomogeneous HNO3 saturation  : " << 100.0 * satHNO3_Hom << " % (Supercooling: " << T_NAT_SUPERCOOL << " [K])";
        std::cout << "\nHeterogeneous HNO3 saturation: " << 100.0 * satHNO3_Het << " %";

        const double SO4_NBIN = std::floor( 1 + log( pow( (SO4_R_HIG/SO4_R_LOW), 3.0 ) ) / log( SO4_VRAT ) );

        std::cout << "\nSulfate's distribution has " << SO4_NBIN << " bins";


        /* SO4 bin radius and volume centers */
        std::vector<double> SO4_rJ( SO4_NBIN, 0.0 );
        std::vector<double> SO4_vJ( SO4_NBIN, 0.0 );

        for ( unsigned int iBin = 0; iBin < SO4_NBIN; iBin++ ) {
            SO4_rJ[iBin] = SO4_R_LOW * pow( SO4_VRAT, iBin / double(3.0) ); /* [m] */
            SO4_vJ[iBin] = SO4_R_LOW * pow( SO4_VRAT, iBin );               /* [m^3] */
        }

        return EPM_SUCCESS;

    } /* End of RunMicrophysics */

    double temp_Vortex( double time, double delta_T )
    {

        /* Returns the temperature increase from cruise altitude as a function of time */
        if ( ( time >= t_Vortex_0 ) )
            return std::min( delta_T, delta_T * ( time - t_Vortex_0 ) / ( t_Vortex_1 - t_Vortex_0 ));
        else
            return 0.0;

    } /* End of Vortex */

}

/* End of Integrate.cpp */
