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

#include "EPM/Integrate.hpp"

namespace EPM
{
    int Integrate( RealDouble temperature_K, RealDouble pressure_Pa, RealDouble relHumidity_w, RealDouble varArray[], RealDouble fixArray[], const Aircraft &AC, const Emission &EI )
    {

        /* Get mean vortex displacement in [m] */
        RealDouble delta_z;
        if ( VORTEX_SINKING )
            delta_z = AC.getVortexdeltaz1(); 
        else
            delta_z = 0.0;

        /* Compute adiabatic and real temperature changes */
        RealDouble delta_T_ad = - physConst::GAMMA_AD * delta_z * 1.0E-03;
        RealDouble delta_T    = -            GAMMA    * delta_z * 1.0E-03;
        /*        [ K ]   =              [ K/km ] *  [ m ]  * [ km/m ] 
         * The minus sign is because delta_z is the distance pointing down */
        std::cout << "\nTemperature lapse rate: " << GAMMA << " [K/km]";
        std::cout << "\nVortex sinking: " << delta_z << " [m]";
        std::cout << "\nTemperature increase: " << delta_T << " [K]\n";

        /************************************************************************************/
        temperature_K = 210.0;
        relHumidity_w = 60.0;
        varArray[ind_H2O] *= 2;
        /************************************************************************************/

        RunMicrophysics( temperature_K, pressure_Pa, relHumidity_w, varArray, fixArray, AC, EI, delta_T_ad, delta_T );

        return EPM_SUCCESS;

    } /* End of Integrate */

    int RunMicrophysics( RealDouble temperature_K, RealDouble pressure_Pa, RealDouble relHumidity_w, RealDouble varArray[], RealDouble fixArray[], const Aircraft &AC, const Emission &EI, RealDouble delta_T_ad, RealDouble delta_T )
    {
    
        RealDouble relHumidity_i_Amb, relHumidity_i_postVortex, relHumidity_i_Final;
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

        RealDouble partPHNO3_Hom, partPHNO3_Het, satHNO3_Hom, satHNO3_Het;
        RealDouble final_Temp = temperature_K + delta_T;
        RealDouble offset_Temp = final_Temp + T_NAT_SUPERCOOL;
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
        std::cout << "\nHeterogeneous HNO3 saturation: " << 100.0 * satHNO3_Het << " %\n";

        const RealDouble SO4_NBIN = std::floor( 1 + log( pow( (SO4_R_HIG/SO4_R_LOW), 3.0 ) ) / log( SO4_VRAT ) );


        if ( SO4_VRAT <= 1.0 ) {
            std::cout << "\nVolume ratio of consecutive bins for SO4 has to be greater than 1.0 ( SO4_VRAT = " << SO4_VRAT << " )";
        }

        /* SO4 bin radius and volume centers */
        Vector_1D SO4_rJ( SO4_NBIN, 0.0 );
        Vector_1D SO4_rE( SO4_NBIN+1, 0.0 );
        Vector_1D SO4_vJ( SO4_NBIN, 0.0 );

        for ( UInt iBin = 0; iBin < SO4_NBIN+1; iBin++ ) {
            SO4_rE[iBin] = SO4_R_LOW * pow( SO4_VRAT, iBin / RealDouble(3.0) );                                /* [m] */
        }

        for ( UInt iBin = 0; iBin < SO4_NBIN; iBin++ ) {
            SO4_rJ[iBin] = 0.5 * ( SO4_rE[iBin] + SO4_rE[iBin+1] );                                            /* [m] */
            SO4_vJ[iBin] = 4.0 / RealDouble(3.0) * physConst::PI * SO4_rJ[iBin] * SO4_rJ[iBin] * SO4_rJ[iBin]; /* [m^3] */
        }

        const AIM::Coagulation Kernel( "liquid", SO4_rJ, SO4_vJ, physConst::RHO_SO4, temperature_K, pressure_Pa );
        //const AIM::Coagulation Kernel( "liquid", SO4_rJ, SO4_vJ, physConst::RHO_SO4, EI.getSootRad(), physConst::RHO_SOOT, temperature_K, pressure_Pa );

        /* For debug purposes */
        if ( DEBUG_COAGKERNEL )
            Kernel.printKernel_1D( "Kernel_Debug_File.out" );

        /* Create SO4 aerosol number distribution.
         * We allocate the PDF with a very small number of existing particle */

        RealDouble nSO4 = 1.00E-10; /* [#/cm^3] */
        RealDouble rSO4 = 5.00E-10; /* [m] */
        RealDouble sSO4 = 1.4; /* For lognormal distributions, sSO4 */
//        RealDouble sSO4 = rSO4/10; /* For normal distributions, sSO4 */

        /* Type of sulfate distribution */
        /* lognormal distribution */
        AIM::Aerosol nPDF_SO4( SO4_rJ, SO4_rE, nSO4, rSO4, sSO4, "lognormal" );
        /* normal distribution */
        // AIM::Aerosol nPDF_SO4( SO4_rJ, SO4_rE, nSO4, rSO4, sSO4, "normal" );
        /* power distribution */
        // AIM::Aerosol nPDF_SO4( SO4_rJ, SO4_rE, nSO4, 0.0, 0.0, "power", 2.0 );
        /* generalized gamma distribution */
        // AIM::Aerosol nPDF_SO4( SO4_rJ, SO4_rE, nSO4, 0.0, 0.0, "gamma", 2.0, 1.0, 1.00E+09 );
       

        /* Store ambient concentrations */
        RealDouble H2O_amb  = varArray[ind_H2O]; /* [molec/cm^3] */
        RealDouble SO4_amb  = varArray[ind_SO4]; /* [molec/cm^3] */
        RealDouble HNO3_amb = varArray[ind_HNO3]; /* [molec/cm^3] */
        RealDouble Soot_amb = 0.05; /* [#/cm^3] */

        /* Add emissions of one engine to concentration arrays */
        varArray[ind_H2O] += EI.getH2O() / ( MW_H2O  * 1.0E+03 ) * AC.getFuelFlow() / RealDouble(AC.getEngNumber()) / AC.getVFlight() * physConst::Na / Ac0     * 1.00E-06;
        /* [ molec/cm^3 ] += [ g/kgf ]   / [ kg/mol ] * [ g/kg ] * [ kgf/s ]                                        / [ m/s ]         * [ molec/mol ] / [ m^2 ] * [ m^3/cm^3 ]
         *                += [ molec/cm^3 ] */
        varArray[ind_SO4] += 0.005 * EI.getSO2() / ( MW_SO4  * 1.0E+03 ) * AC.getFuelFlow() / RealDouble(AC.getEngNumber()) / AC.getVFlight() * physConst::Na / Ac0 * 1.00E-06;


        RealDouble varSoot = Soot_amb + EI.getSoot() / ( 4.0 / RealDouble(3.0) * physConst::PI * physConst::RHO_SOOT * 1.00E+03 * EI.getSootRad() * EI.getSootRad() * EI.getSootRad() ) * AC.getFuelFlow() / RealDouble(AC.getEngNumber()) / AC.getVFlight() * 1.00E-06; /* [ #/cm^3 ] */

        /* Array storing the current value of the variable */
        RealDouble vars[N_MICROVAR];
        
        /* Adaptive time stepping? */
        bool adaptiveStep = 1;
        
        UInt nTime = 201;
        UInt iTime;
        Vector_1D timeArray(nTime);
        RealDouble timeInitial = 1.0E-04;
        RealDouble log10_timeInitial = log10(timeInitial);
        RealDouble timeFinal   = 5.0E+03;
        RealDouble log10_timeFinal = log10(timeFinal);
        for ( iTime = 0; iTime < nTime; iTime++ ) {
            timeArray[iTime] = pow( 10.0, log10_timeInitial + iTime * ( log10_timeFinal - log10_timeInitial ) / RealDouble( nTime - 1.0 ) );
        }

        /* vars[0] : Temperature,                   Unit K 
         * vars[1] : Water molecular concentration, Unit molecules/cm^3 
         * vars[2] : Gaseous SO4 concentration,     Unit molecules/cm^3
         * vars[3] : Liquid SO4 concentration,      Unit molecules/cm^3
         * vars[4] : Gaseous HNO3 concentration,    Unit molecules/cm^3
         * vars[5] : Liquid HNO3 concentration,     Unit molecules/cm^3
         * vars[6] : Soot ambient concentration,    Unit #/cm^3
         * vars[7] : Ice particle radius,           Unit m
         */
            
        static const UInt EPM_ind_T    = 0;
        static const UInt EPM_ind_P    = 1;
        static const UInt EPM_ind_H2O  = 2;
        static const UInt EPM_ind_SO4  = 3;
        static const UInt EPM_ind_HNO3 = 4;
        static const UInt EPM_ind_Part = 5;
        static const UInt EPM_ind_ParR = 6;

        struct gas_aerosol_rhs
        {
          
            const RealDouble m_temperature_K;
            const RealDouble m_pressure_Pa;
            const RealDouble m_delta_T;
            const RealDouble m_H2O_molcm3;
            const RealDouble m_SO4_molcm3;
            const RealDouble m_HNO3_molcm3;
            const RealDouble m_part_cm3;
            const RealDouble m_part_r0;

            const RealDouble sticking_SO4;
            const RealDouble sigma_SO4;

            private:

//                const UInt EPM_ind_T    = 0;
//                const UInt EPM_ind_P    = 1;
//                const UInt EPM_ind_H2O  = 2;
//                const UInt EPM_ind_SO4  = 3;
//                const UInt EPM_ind_HNO3 = 4;
//                const UInt EPM_ind_Part = 5;
//                const UInt EPM_ind_ParR = 6;


            public:
            
                gas_aerosol_rhs( RealDouble temperature_K, RealDouble pressure_Pa, RealDouble delta_T, \
                                 RealDouble H2O_molcm3, RealDouble SO4_molcm3, RealDouble HNO3_molcm3, \
                                 RealDouble part_cm3, RealDouble part_r0 ):
                    m_temperature_K( temperature_K ),
                    m_pressure_Pa( pressure_Pa ),
                    m_delta_T( delta_T ),
                    m_H2O_molcm3 ( H2O_molcm3 ),
                    m_SO4_molcm3 ( SO4_molcm3 ),
                    m_HNO3_molcm3( HNO3_molcm3 ),
                    m_part_cm3( part_cm3 ),
                    m_part_r0( part_r0 ),
                    sticking_SO4( 1.0 ),
                    sigma_SO4( 5.0E+14 )
                {

                    /* Default Constructor */

                } /* End of gas_aerosol_rhs::gas_aerosol_rhs */
        
                void operator()( const Vector_1D &x, Vector_1D &dxdt, const RealDouble t = 0 ) const
                {

                    /* DESCRIPTION:
                     * Returns the right hand side of the system of ordinary differential equations:
                     * y' = f(t,y)
                     * where y belongs to an n-dimensional Euclidean space and t is time.
                     * If f doesn't explicitly depend on time, the system is autonomous and the parameter t can be dropped
                     * Default t value is 0. */
            
                    /* INPUT:
                     * - Vector_1D x    :: vector of variables
                     * - Vector_1D dxdt :: vector of rate of change
                     *
                     * (optional)
                     * - RealDouble t   :: time expressed in s
                     *
                     * OUTPUT:
                     */

                    RealDouble dilRatio = dilutionRatio( t );

                    /* Ensure that particle radius is larger than bare soot radius */
                    RealDouble m_part_r = x[EPM_ind_ParR] > m_part_r0 ? x[EPM_ind_ParR] : m_part_r0;
                        // std::max( x[EPM_ind_ParR], m_part_r0 );

                    /* Temperature differential equation
                     * \frac{dT}{dt} = v_z(t) * \frac{dT}{dz} - w_T(t) * ( T - T_amb(t) ) */
                    dxdt[EPM_ind_T] = dT_Vortex( t, m_delta_T, 1 ) \
                                    - dilRatio * ( x[EPM_ind_T] - ( m_temperature_K + dT_Vortex ( t, m_delta_T ) ) );
                    /* Unit check:
                     * dT_Vortex(.,.,1): [K/s]
                     * dT_Vortex(.,.,0): [K]
                     * dilRatio * ...  : [1/s] * ( [K] - ( [K] + [K] ) ) = [K/s] */

                    /* Pressure differential equation
                     * \frac{dP}{dt} = 0 */
                    dxdt[EPM_ind_P] = 0.0;

                    /* Gaseous water differential equation 
                     * \frac{d[H2O]}{dt} = - w_T(t) * ( [H2O] - [H2O]_amb ) */
                    dxdt[EPM_ind_H2O] = - dilRatio * ( x[EPM_ind_H2O] - m_H2O_molcm3 ) \
                                      - ( isFreezable( m_part_r, x[EPM_ind_T], x[EPM_ind_P], x[EPM_ind_H2O], m_part_r0 ) \
                                        * depositionRate( m_part_r, x[EPM_ind_T], x[EPM_ind_P], x[EPM_ind_H2O] ) \
                                        + condensationRate( m_part_r, x[EPM_ind_T], x[EPM_ind_P], x[EPM_ind_H2O], x[9] + x[10] ) ) \
                                        * m_part_r * physConst::Na / MW_H2O;
                    /* Unit check:
                     * [1/s] * ( [molec/cm3] - [molec/cm^3] ) = [molec/cm^3/s]
                     * ([-] * [kg/s] + [kg/s]) * [#/cm^3] * [molec/kg] = [molec/cm^3/s]
                     * */

                    /* Gaseous SO4 differential equation 
                     * \frac{d[SO4]}{dt} = - w_T(t) * ( [SO4] - [SO4]_amb ) */
                    dxdt[EPM_ind_SO4] = - dilRatio * ( x[EPM_ind_SO4] - m_SO4_molcm3 );
                    
                    /* Gaseous HNO3 differential equation 
                     * \frac{d[HNO3]}{dt} = - w_T(t) * ( [HNO3] - [HNO3]_amb ) */
                    dxdt[EPM_ind_HNO3] = - dilRatio * ( x[EPM_ind_HNO3] - m_HNO3_molcm3 );

                    /* Particle differential equation 
                     * \frac{d[part]}{dt} = - w_T(t) * ( [part] - [part]_amb )*/
                    dxdt[EPM_ind_Part] = - dilRatio * ( x[EPM_ind_Part] - m_part_cm3 );

                    /* Particle radius differential equation
                     * \frac{dr}{dt} = \frac{dm}{dt}/(\rho*4*PI*r^2) */
                    dxdt[EPM_ind_ParR] = ( m_part_r >= 0.99 * m_part_r0 ) \
                                         * ( isFreezable( m_part_r, x[EPM_ind_T], x[EPM_ind_P], x[EPM_ind_H2O], m_part_r0 ) \
                                           * depositionRate( m_part_r, x[EPM_ind_T], x[EPM_ind_P], x[EPM_ind_H2O] ) \
                                           + condensationRate( m_part_r, x[EPM_ind_T], x[EPM_ind_P], x[EPM_ind_H2O], x[9] + x[10] ) ) \
                                         / ( physConst::RHO_ICE * 4.0 * physConst::PI * m_part_r * m_part_r );
                    /* Unit check: 
                     * [kg/s] / ( [kg/m^3] * [m^2]) = [m/s] */

                    dxdt[7] = 0.0;

                    dxdt[8] = 0.0;

                    /* */
                    dxdt[9] = sticking_SO4 * physFunc::thermalSpeed( x[0], MW_SO4 / physConst::Na ) * 1.0E+02 / 4 * x[3] / sigma_SO4 * ( 1.0 - x[9] - x[10] );

                    dxdt[10] = 0.0;

//                    std::cout << std::scientific << std::setprecision(6);
//                    std::cout << t << ", ";
//                    std::cout << x[2] * physConst::kB * x[0] * 1.0E+06 / physFunc::pSat_H2Os( x[0] ) << ", ";
//                    std::cout << dxdt[6] << ", ";
//                    std::cout << m_H2O_molcm3 * physConst::kB * m_temperature_K * 1.0E+06 / physFunc::pSat_H2Os( m_temperature_K ) << "\n";

                } /* End of gas_aerosol_rhs::operator() */

        };


        Vector_1D x(11);
        x[EPM_ind_T]    = Tc0;
        x[EPM_ind_P]    = pressure_Pa;
        x[EPM_ind_H2O]  = varArray[ind_H2O];
        x[EPM_ind_SO4]  = varArray[ind_SO4];
        x[EPM_ind_HNO3] = varArray[ind_HNO3];
        x[EPM_ind_Part] = varSoot;
        x[EPM_ind_ParR] = EI.getSootRad();
        x[9] = 0.0;

        Vector_2D obs_Var;
        Vector_1D obs_Time;
        EPM::streamingObserver observer( obs_Var, obs_Time, "/home/fritzt/CAPCEMM/data/Micro.out", 1 );

        /* Creating ode's right hand side */
        gas_aerosol_rhs rhs( temperature_K, pressure_Pa, delta_T, H2O_amb, SO4_amb, HNO3_amb, Soot_amb, EI.getSootRad() );
        
        /* Total number of steps */
        UInt totSteps = 0;

        RealDouble attempt_dt;

        iTime = 0;
        while ( iTime + 1 < nTime ) {
                
            attempt_dt = ( timeArray[iTime+1] - timeArray[iTime] ) / 100.0;
            if ( adaptiveStep == 1 ) {
                totSteps += boost::numeric::odeint::integrate_adaptive( boost::numeric::odeint::make_controlled< error_stepper_type >( EPM_ATOLS, EPM_RTOLS ), rhs, x, timeArray[iTime], timeArray[iTime+1], attempt_dt, observer );
            }
            else {
                totSteps += boost::numeric::odeint::integrate( rhs, x, timeArray[iTime], timeArray[iTime+1], attempt_dt, observer );
            }

            totSteps += 1;
            iTime++;
        }
        
        std::cout << "\nTotal number of steps: " << totSteps << "\n";

        observer.print2File();

        return EPM_SUCCESS;

    } /* End of RunMicrophysics */

    RealDouble dT_Vortex( const RealDouble time, const RealDouble delta_T, bool deriv )
    {

        /* DESCRIPTION:
         * Returns the temperature increase / rate of change from cruise altitude as a
         * function of time */

        /* INPUT:
         * - RealDouble time    :: time since engine exit plane expressed in s
         * - RealDouble delta_T :: temperature increase in K after vortex sinking in K
         * - bool deriv         :: rate of change or delta?
         *
         * OUTPUT:
         * - RealDouble :: temperature increase / rate of change from cruise altitude in K */

        if ( deriv == 0 ) {
            if ( ( time >= t_Vortex_0 ) )
                return std::min( delta_T, delta_T * ( time - t_Vortex_0 ) / ( t_Vortex_1 - t_Vortex_0 ) );
            else
                return 0.0;
        } else {
            if ( ( time >= t_Vortex_0 ) && ( time <= t_Vortex_1 ) )
                return delta_T / ( t_Vortex_1 - t_Vortex_0 );
            else
                return 0.0;
        }

    } /* End of dT_Vortex */
    
    RealDouble dilutionRatio( const RealDouble t )
    {

        /* DESCRIPTION:
         * Returns the plume dilution ratio in the plume early stages */

        /* INPUT:
         * - RealDouble t :: time since engine exit plane expressed in s
         *
         * OUTPUT:
         * - RealDouble :: plume dilution ratio [-] between 0 and 1 */

        RealDouble dRat = 0;

        if ( ( t >= t_0 ) && ( t <= t_j ) ) {
            if ( t < 1.0E-03 ) {
                /* Obtained from best fit */
                dRat = exp( 12.175448395306601 * log(t) + 83.748342982715030 );
            } else if ( t >= 1.0E-03 ) { 
                dRat = std::max( 0.8 / t_j - 0.8 / ( t_j * t_j ) * ( t - t_j ) + ( 1 - 0.8 / t_j + 0.8 / ( t_j * t_j ) * ( 1.0E-03 - t_j ) ) / ( ( 1.0E-03 - t_j ) * ( 1.0E-03 - t_j ) ) * ( t - t_j ) * ( t - t_j ), 1.0E-10 );
            }
        } else if ( ( t > t_j ) && ( t <= t_1 ) ) {
            dRat = 0.8 / t ;
        } else if ( ( t > t_1 ) && ( t <= t_2 ) ) {
            dRat = ( 0.8 / t_1 ) / ( 1.0 + Cv * ( 0.8 / t_1 ) * ( t - t_1 ) );
        } else if ( ( t > t_2 ) ) {
            dRat = ( 0.8 / t_1 ) / ( 1.0 + Cv * ( 0.8 / t_1 ) * ( t_2 - t_1 ) ) / ( 1.0 + 1.0 / ( 1.0 + m * exp( - ( t - t_2 ) / ( n * t_2 ) ) ) * ( 0.8 / t_1 ) / ( 1.0 + Cv * ( 0.8 / t_1 ) * ( t_2 - t_1 ) ) * ( t - t_2 ) );
        } else {
            dRat = 0;
        }

        return dRat;

    } /* End of dilutionRatio */

    RealDouble depositionRate( const RealDouble r, const RealDouble T, const RealDouble P, const RealDouble H2O )
    {

        /* DESCRIPTION:
         * Returns the water deposition rate in kg/s on a single spherical particle */

        /* INPUT:
         * - RealDouble r   :: radius in m
         * - RealDouble T   :: temperature in K
         * - RealDouble P   :: pressure in Pa
         * - RealDouble H2O :: gaseous water concentration in molec/cm^3
         *
         * OUTPUT:
         * - RealDouble :: deposition rate */

        /* Capacitance factor for spherical particles is C = r; */

        return 4.0 * physConst::PI * r * physFunc::CorrDiffCoef_H2O( r, T, P ) * \
               ( H2O * physConst::kB * T * 1.0E+06 - physFunc::pSat_H2Os( T ) * physFunc::Kelvin( r ) ) / \
               ( physConst::R / MW_H2O * T + physFunc::LHeatSubl_H2O( T ) * physFunc::CorrDiffCoef_H2O( r, T, P ) * \
                 physFunc::pSat_H2Os( T ) * ( physFunc::LHeatSubl_H2O( T ) * MW_H2O / ( physConst::R * T ) - 1.0 ) / \
                   ( physFunc::ThermalCond( r, T, P ) * T ) );


    } /* End of depositionRate */

    bool isFreezable( const RealDouble r, const RealDouble T, const RealDouble P, const RealDouble H2O, const RealDouble r0 )
    {

        /* DESCRIPTION:
         * Returns whether or not water can freeze on a particle */

        /* INPUT:
         * - RealDouble r   :: radius in m
         * - RealDouble T   :: temperature in K
         * - RealDouble P   :: pressure in Pa
         * - RealDouble H2O :: gaseous water concentration in molec/cm^3
         * - RealDouble r0  :: bare soot radius in m
         *
         * OUTPUT:
         * - bool :: isFreezable? */

        /* Returns 1 if: 
         * - Current conditions are supersaturated with respect to liquid water (freezing)
         * or 
         * - Current radius is greater than bare soot radius (melting) */

        return ( ( H2O * physConst::kB * T * 1.0E+06 >= physFunc::Kelvin( r ) * physFunc::pSat_H2Ol( T ) ) || ( r > 1.005 * r0 ) );

    } /* End of isFreezable */

    RealDouble condensationRate( const RealDouble r, const RealDouble T, const RealDouble P, const RealDouble H2O, const RealDouble theta )
    {

        /* DESCRIPTION:
         * Returns the gaseous water condensation rate in kg/s on a signe spherical particle */

        /* INPUT:
         * - RealDouble r     :: radius in m
         * - RealDouble T     :: temperature in K
         * - RealDouble P     :: pressure in Pa
         * - RealDouble H2O   :: gaseous water concentration in molec/cm^3
         * - RealDouble theta :: fractional coating of the particle
         *
         * OUTPUT:
         * - RealDouble :: condensation rate */

        if ( H2O * physConst::kB * T * 1.0E+06 >= physFunc::Kelvin( r ) * physFunc::pSat_H2Ol( T ) ) {
            return 4.0 * physConst::PI * MW_H2O * physFunc::CorrDiffCoef_H2O( r, T, P ) * r * theta \
                       * ( H2O / physConst::Na * 1.0E+06 \
                         - physFunc::Kelvin( r ) * physFunc::pSat_H2Ol( T ) / ( physConst::R * T ) );
        } else {
            return 0.0;
        }

    } /* End of condensationRate */

}

/* End of Integrate.cpp */
