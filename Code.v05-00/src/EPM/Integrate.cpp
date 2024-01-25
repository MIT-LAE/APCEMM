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
    int Integrate( double &temperature_K, double pressure_Pa, double relHumidity_w, double varArray[], \
                   const Vector_2D& aerArray, const Aircraft &AC, const Emission &EI, \
                   double &Ice_rad, double &Ice_den, double &Soot_den, double &H2O_mol, \
                   double &SO4g_mol, double &SO4l_mol, AIM::Aerosol &SO4Aer, AIM::Aerosol &IceAer, \
                   double &Area, double &Ab0, double &Tc0, const bool CHEMISTRY, std::string micro_data_out )
    {

        /* Get mean vortex displacement in [m] */
        double delta_z;
        if ( VORTEX_SINKING )
            delta_z = AC.deltaz1(); 
        else
            delta_z = 0.0;

        /* Compute adiabatic and real temperature changes */
        double delta_T_ad = - physConst::GAMMA_AD * delta_z * 1.0E-03;
        double delta_T    = -            GAMMA    * delta_z * 1.0E-03;
        /*          [ K ]     =              [ K/km ] *  [ m ]  * [ km/m ] 
         * The minus sign is because delta_z is the distance pointing down */

        int EPM_RC = RunMicrophysics( temperature_K, pressure_Pa, relHumidity_w, varArray, aerArray, AC, EI, delta_T_ad, delta_T, \
                                      Ice_rad, Ice_den, Soot_den, H2O_mol, SO4g_mol, SO4l_mol, SO4Aer, IceAer, Area, Ab0, Tc0, CHEMISTRY, micro_data_out );

        return EPM_RC;

    } /* End of Integrate */

    /* TODO: Make the original integrate function work with the new EPMOutput struct directly, and then delete this function.*/
    std::pair<EPMOutput, int> Integrate(double tempInit_K, double pressure_Pa, double rhw, double bypassArea, double coreExitTemp, double varArray[], 
                            const Vector_2D& aerArray, const Aircraft& AC,const Emission& EI, bool CHEMISTRY, 
                            std::string micro_data_out) 
    {
        EPMOutput out;
        out.finalTemp = tempInit_K;
        out.bypassArea = bypassArea;
        out.coreExitTemp = coreExitTemp;
        int returnCode = Integrate(out.finalTemp, pressure_Pa, rhw, varArray, aerArray, AC, EI, out.iceRadius,
                                    out.iceDensity, out.sootDensity, out.H2O_mol, out.SO4g_mol, out.SO4l_mol,
                                    out.SO4Aer, out.IceAer, out.area, out.bypassArea, out.coreExitTemp, CHEMISTRY, micro_data_out);
        return std::make_pair(out, returnCode);
    }

    /* FIXME: See above comment on integrate() */
    int RunMicrophysics( double &temperature_K, double pressure_Pa, double relHumidity_w, double varArray[], \
                        const Vector_2D& aerArray, const Aircraft &AC, const Emission &EI, \
                         double delta_T_ad, double delta_T, double &Ice_rad, double &Ice_den, \
                         double &Soot_den, double &H2O_mol, double &SO4g_mol, double &SO4l_mol, \
                         AIM::Aerosol &SO4Aer, AIM::Aerosol &IceAer, double &Area, double &Ab0, double &Tc0, const bool CHEMISTRY, std::string micro_data_out )
    {
    
        double relHumidity_i_Amb, relHumidity_i_postVortex, relHumidity_i_Final;

        /* relHumidity_w is in % */
        relHumidity_i_Amb        = relHumidity_w * physFunc::pSat_H2Ol( temperature_K ) \
                                                 / physFunc::pSat_H2Os( temperature_K ) / 100.0;
        relHumidity_i_postVortex = relHumidity_w * physFunc::pSat_H2Ol( temperature_K + delta_T_ad ) \
                                                 / physFunc::pSat_H2Os( temperature_K + delta_T_ad ) / 100.0;
        relHumidity_i_Final      = relHumidity_w * physFunc::pSat_H2Ol( temperature_K + delta_T ) \
                                                 / physFunc::pSat_H2Os( temperature_K + delta_T ) / 100.0;

        
        double partPHNO3_Hom, partPHNO3_Het, satHNO3_Hom, satHNO3_Het;
        double final_Temp = temperature_K + delta_T;
        double offset_Temp = final_Temp + T_NAT_SUPERCOOL; // NAT = Nitric acid trihydrate

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
        
        /* Number of SO4 size distribution bins based on specified min and max radii *
        * and volume ratio between two adjacent bins */
        const UInt SO4_NBIN = std::floor( 1 + log( pow( (LA_R_HIG/LA_R_LOW), 3.0 ) ) / log( LA_VRAT ) );

        if ( LA_VRAT <= 1.0 ) {
            std::cout << "\nVolume ratio of consecutive bins for SO4 has to be greater than 1.0 ( LA_VRAT = " << LA_VRAT << " )";
        }

        /* SO4 bin radius and volume centers */
        Vector_1D SO4_rJ( SO4_NBIN    , 0.0 ); // bin centers, radius
        Vector_1D SO4_rE( SO4_NBIN + 1, 0.0 ); // bin edges, radius
        Vector_1D SO4_vJ( SO4_NBIN    , 0.0 ); // bin centers, volume 

        /* Adjacent bin radius ratio */ 
        const double LA_RRAT = pow( LA_VRAT, 1.0 / double(3.0) );

        /* Initialize bin center and edge radii, as well as volume */
        SO4_rE[0] = LA_R_LOW;
        for ( UInt iBin = 1; iBin < SO4_NBIN + 1; iBin++ )
            SO4_rE[iBin] = SO4_rE[iBin-1] * LA_RRAT;                                                           /* [m] */

        for ( UInt iBin = 0; iBin < SO4_NBIN; iBin++ ) {
            SO4_rJ[iBin] = 0.5 * ( SO4_rE[iBin] + SO4_rE[iBin+1] );                                            /* [m] */
            SO4_vJ[iBin] = 4.0 / double(3.0) * physConst::PI * SO4_rJ[iBin] * SO4_rJ[iBin] * SO4_rJ[iBin]; /* [m^3] */
        }
        
        /* Number of ice size distribution bins based on specified min and max radii *
        * and volume ratio between two adjacent bins */
        const UInt Ice_NBIN = std::floor( 1 + log( pow( (PA_R_HIG/PA_R_LOW), 3.0 ) ) / log( PA_VRAT ) );
        
        /* Adjacent bin radius ratio */ 
        const double PA_RRAT = pow( PA_VRAT, 1.0 / double(3.0) );

        /* Ice bin center and edge radii */
        Vector_1D Ice_rJ( Ice_NBIN    , 0.0 );
        Vector_1D Ice_rE( Ice_NBIN + 1, 0.0 );
        Ice_rE[0] = PA_R_LOW;
        for ( UInt iBin = 1; iBin < Ice_NBIN + 1; iBin++ )
            Ice_rE[iBin] = Ice_rE[iBin-1] * PA_RRAT;                                                           /* [m] */

        for ( UInt iBin = 0; iBin < Ice_NBIN; iBin++ )
            Ice_rJ[iBin] = 0.5 * ( Ice_rE[iBin] + Ice_rE[iBin+1] );                                            /* [m] */


        /* Create coagulation kernels */ 
        const AIM::Coagulation Kernel( "liquid", SO4_rJ, SO4_vJ, physConst::RHO_SULF, temperature_K, pressure_Pa );
        const AIM::Coagulation Kernel_SO4_Soot( "liquid", SO4_rJ, physConst::RHO_SULF, EI.getSootRad(), physConst::RHO_SOOT, temperature_K, pressure_Pa );
       
        const Vector_1D KernelSO4Soot = Kernel_SO4_Soot.getKernel_1D();

        /* Create SO4 aerosol number distribution.
         * We allocate the PDF with a very small number of existing particles */
        
        double nSO4 = 1.00E-10; /* [#/cm^3] */
        double rSO4 = 5.00E-10; /* [m] */
        double sSO4 = 1.4; /* For lognormal distributions, sSO4 = sigma */
//        double sSO4 = rSO4/10; /* For normal distributions, sSO4 */

        /* Type of sulfate distribution */
        /* lognormal distribution */

        /* FIXME: Why do we need a thread local static variable for this? Makes literally 0 sense.
        There is nothing within this function to parallelize. - Michael */
        thread_local static AIM::Aerosol nPDF_SO4( SO4_rJ, SO4_rE, nSO4, rSO4, sSO4, "lognormal" );

        /* normal distribution */
        // AIM::Aerosol nPDF_SO4( SO4_rJ, SO4_rE, nSO4, rSO4, sSO4, "normal" );
        /* power distribution */
        // AIM::Aerosol nPDF_SO4( SO4_rJ, SO4_rE, nSO4, 0.0, 0.0, "power", 2.0 );
        /* generalized gamma distribution */
        // AIM::Aerosol nPDF_SO4( SO4_rJ, SO4_rE, nSO4, 0.0, 0.0, "gamma", 2.0, 1.0, 1.00E+09 );
       

        double n_air_amb = physConst::Na * pressure_Pa / (physConst::R * temperature_K * 1.0e6) ;
        double n_air_eng = physConst::Na * pressure_Pa / (physConst::R * Tc0 * 1.0e6) ;


        /* Store ambient concentrations */
        double H2O_amb  = varArray[ind_H2O]/n_air_amb;     /* [molec/cm^3] */
        double SO4_amb  = varArray[ind_SO4]/n_air_amb;     /* [molec/cm^3] */
        double SO4l_amb, SO4g_amb;               /* [molec/cm^3] */
        double HNO3_amb = varArray[ind_HNO3]/n_air_amb;    /* [molec/cm^3] */
        double Soot_amb = aerArray[ind_SOOT][0]/n_air_amb; ///n_air_amb; /* [#/cm^3] */

        /* Add emissions of one engine to concentration arrays */
        varArray[ind_H2O] += EI.getH2O() / ( MW_H2O  * 1.0E+03 ) * AC.FuelFlow() / double(AC.EngNumber()) / AC.VFlight() * physConst::Na / Ab0     * 1.00E-06 ;
        /* [ molec/cm^3 ] += [ g/kgf ]   / [ kg/mol ] * [ g/kg ] * [ kgf/s ]                                  / [ m/s ]      * [ molec/mol ] / [ m^2 ] * [ m^3/cm^3 ]
         *                += [ molec/cm^3 ] */

        varArray[ind_H2O] = varArray[ind_H2O] / n_air_eng ;

        /* Fixed SO2 */
        // varArray[ind_SO4] += SO2TOSO4 * 0.8 / ( MW_H2SO4  * 1.0E+03 ) * AC.FuelFlow() / double(AC.EngNumber()) / AC.VFlight() * physConst::Na / Ab0 * 1.00E-06;
        /* Variable SO2 */
        varArray[ind_SO4] += SO2TOSO4 * EI.getSO2() / ( MW_H2SO4  * 1.0E+03 ) * AC.FuelFlow() / double(AC.EngNumber()) / AC.VFlight() * physConst::Na / Ab0 * 1.00E-06;
        varArray[ind_SO4] = varArray[ind_SO4] / n_air_eng;

        double varSoot = Soot_amb + EI.getSoot() / ( 4.0 / double(3.0) * physConst::PI * physConst::RHO_SOOT * 1.00E+03 * EI.getSootRad() * EI.getSootRad() * EI.getSootRad() ) * AC.FuelFlow() / double(AC.EngNumber()) / AC.VFlight() / Ab0 * 1.00E-06 ; /* [ #/cm^3 ] */
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
        for ( iTime = 0; iTime < nTime; iTime++ ) {
            timeArray[iTime] = pow( 10.0, log10_timeInitial + iTime * ( log10_timeFinal - log10_timeInitial ) / double( nTime - 1.0 ) );
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
        
        /* TODO: Declare these in a different file (and make them constexpr). -Michael */
        static const UInt EPM_ind_Trac = 0;
        static const UInt EPM_ind_T    = 1;
        static const UInt EPM_ind_P    = 2;
        static const UInt EPM_ind_H2O  = 3;
        static const UInt EPM_ind_SO4  = 4;
        static const UInt EPM_ind_SO4l = 5;
        static const UInt EPM_ind_SO4g = 6;
        static const UInt EPM_ind_SO4s = 7;
        static const UInt EPM_ind_HNO3 = 8;
        static const UInt EPM_ind_Part = 9;
        static const UInt EPM_ind_ParR = 10;
        static const UInt EPM_ind_the1 = 11;
        static const UInt EPM_ind_the2 = 12;

        const std::vector<UInt> EPM_ind = {EPM_ind_Trac, EPM_ind_T, EPM_ind_P, EPM_ind_H2O, EPM_ind_SO4, EPM_ind_SO4l, EPM_ind_SO4g, EPM_ind_SO4s, EPM_ind_HNO3, EPM_ind_Part, EPM_ind_ParR, EPM_ind_the1, EPM_ind_the2};

        /* FIXME:  Declaring a huge struct like this in the middle of a function is not ideal. nuff said.
        Move to this a different file. -Michael*/
        struct gas_aerosol_rhs
        {
          
            const double m_temperature_K;
            const double m_pressure_Pa;
            const double m_delta_T;
            const double m_H2O_mixingratio;
            const double m_SO4_mixingratio;
            const double m_SO4l_mixingratio;
            const double m_SO4g_mixingratio;
            const double m_HNO3_mixingratio;
            const double m_part_mixingratio;
            const double m_part_r0;

            const double sticking_SO4;
            const double sigma_SO4;

            const Vector_1D KernelSO4Soot;


            private:

            public:
            
                gas_aerosol_rhs( double temperature_K, double pressure_Pa, double delta_T, \
                                 double H2O_mixingratio, double SO4_mixingratio, double SO4l_mixingratio, \
                                 double SO4g_mixingratio, double HNO3_mixingratio, double part_mixingratio, \
                                 double part_r0, Vector_1D Kernel_ ):
                    m_temperature_K( temperature_K ),
                    m_pressure_Pa( pressure_Pa ),
                    m_delta_T( delta_T ),
                    m_H2O_mixingratio ( H2O_mixingratio ),
                    m_SO4_mixingratio ( SO4_mixingratio ),
                    m_SO4l_mixingratio ( SO4l_mixingratio ),
                    m_SO4g_mixingratio ( SO4g_mixingratio ),
                    m_HNO3_mixingratio ( HNO3_mixingratio ),
                    m_part_mixingratio( part_mixingratio ),
                    m_part_r0( part_r0 ),
                    sticking_SO4( 1.0 ),
                    sigma_SO4( 5.0E+14 ),
                    KernelSO4Soot( Kernel_ )
                {

                    /* Default Constructor */

                } /* End of gas_aerosol_rhs::gas_aerosol_rhs */
        
                void operator()( const Vector_1D &x, Vector_1D &dxdt, const double t = 0 ) const
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
                     * - double t   :: time expressed in s
                     *
                     * OUTPUT:
                     */

                    double entrainRate = entrainmentRate( t );

                   

                    /* Compute air number concentration to convert number
                    concentration rates to mixing ratio rates */ 
                    double n_air = physConst::Na * m_pressure_Pa/(physConst::R * x[EPM_ind_T] * 1.0e6);

                    /* Ensure that particle radius is larger than bare soot radius and smaller than some threshold */
                    double m_part_r;
                    m_part_r = (x[EPM_ind_ParR] > m_part_r0) ? x[EPM_ind_ParR] : m_part_r0;
                    m_part_r = (m_part_r        > 1.00E-04 ) ? 1.00E-04        : m_part_r;
               
                    /* Parameters for binary homogeneous nucleation */
                    double nucRate, x_SO4, nMolec;

                    /* Compute SO4 mixing ratios */
                    
                    /* Liquid and gaseous sulfate concentration */
                    double SO4_rl, SO4_g;
                    
                    /* SO4_g represents the gaseous molecular concentration of SO4
                     * SO4_rl is the "ready to be liquid" molecular concentration of SO4. This is still gaseous SO4 because it is limited by kinetics.
                     * Gaseous SO4 is in phase equilibrium. Ensure limitations! */
                    
                    SO4_g  = ( physFunc::pSat_H2SO4( x[EPM_ind_T] ) / ( physConst::kB * x[EPM_ind_T] * 1.0E+06 ) > x[EPM_ind_SO4g] * n_air ) ? x[EPM_ind_SO4g] : physFunc::pSat_H2SO4( x[EPM_ind_T] ) / ( physConst::kB * x[EPM_ind_T] * 1.0E+06 * n_air); 
                    
                    // SO4_rl = x[EPM_ind_SO4g] * x[EPM_ind_P] / ( physConst::kB * x[EPM_ind_T] * 1.0E+06 ) - SO4_g;
                    SO4_rl = x[EPM_ind_SO4g] - SO4_g; 
                    SO4_g  = ( SO4_g > 0.0 ) ? SO4_g : 0.0;
                    SO4_rl = ( SO4_rl > 0.0 ) ? SO4_rl : 0.0;
                    

                    double RH_liquid = n_air * x[EPM_ind_H2O] * 1.0e6 * physConst::kB * x[EPM_ind_T] / physFunc::pSat_H2Ol( x[EPM_ind_T] ) ; 
                    if ( ( n_air * SO4_rl >= AIM::nThresh( x[EPM_ind_T], RH_liquid))) {
                       
                        /* Mole fraction of sulfuric acid */
                        x_SO4 = AIM::x_star(x[EPM_ind_T], RH_liquid, n_air*SO4_rl); 

                        /* Number of molecules per cluster */
                        nMolec = AIM::nTot(      x[EPM_ind_T], x_SO4, RH_liquid, n_air*SO4_rl );
                        
                        /* Nucleation rate */
                        nucRate = AIM::nuclRate( x[EPM_ind_T], x_SO4, RH_liquid, n_air*SO4_rl );

                    } else {

                        x_SO4 = 0.0;
                        nMolec = 0.0;
                        nucRate = 0.0;

                    }
                    
		            /* Compute sulfate - soot coagulation rate */
                    double CoagRate = 0;
                    Vector_1D pdf = nPDF_SO4.getPDF();
                    Vector_1D binCenters = nPDF_SO4.getBinCenters();
                    for ( unsigned int iBin = 0; iBin < nPDF_SO4.getNBin(); iBin++ ) {
                        CoagRate += ( KernelSO4Soot[iBin] * pdf[iBin] * physConst::PI * binCenters[iBin] * binCenters[iBin] ) * ( 1.0 - x[EPM_ind_the1] - x[EPM_ind_the2] );
                        /* Unit check:
                         * [cm^3/s] * [#/cm3] * [m^2] = [m^2/s]
                         */
                    }

        
                    /* Tracer differential equation
                     * \frac{d[.]}{dt} = - w_T(t) * [.] */
                    dxdt[EPM_ind_Trac] = - entrainRate * x[EPM_ind_Trac];

                    /* Temperature differential equation
                     * \frac{dT}{dt} = v_z(t) * \frac{dT}{dz} - w_T(t) * ( T - T_amb(t) ) */
                    dxdt[EPM_ind_T] = dT_Vortex( t, m_delta_T, 1 ) \
                                    - entrainRate * ( x[EPM_ind_T] - ( m_temperature_K + dT_Vortex ( t, m_delta_T ) ) );
                    /* Unit check:
                     * dT_Vortex(.,.,1): [K/s]
                     * dT_Vortex(.,.,0): [K]
                     * entrainRate * ...  : [1/s] * ( [K] - ( [K] + [K] ) ) = [K/s] */

                    /* Pressure differential equation
                     * \frac{dP}{dt} = 0 */
                    dxdt[EPM_ind_P] = 0.0;

                    /* Gaseous water differential equation 
                     * \frac{d[H2O]}{dt} = - w_T(t) * ( [H2O] - [H2O]_amb ) */
                    dxdt[EPM_ind_H2O] = - entrainRate * ( x[EPM_ind_H2O] - m_H2O_mixingratio ) \
                                        - ( isFreezable( m_part_r, x[EPM_ind_T], x[EPM_ind_H2O]*n_air, m_part_r0 ) \
                                          * depositionRate( m_part_r, x[EPM_ind_T], x[EPM_ind_P], x[EPM_ind_H2O]*n_air, m_part_r0, x[EPM_ind_the1] + x[EPM_ind_the2] ) \
                                        + condensationRate( m_part_r, x[EPM_ind_T], x[EPM_ind_P], x[EPM_ind_H2O]*n_air, x[EPM_ind_the1] + x[EPM_ind_the2] ) ) \
                                          * x[EPM_ind_Part] * physConst::Na / MW_H2O \
                                        -  nucRate * ( 1.0 - x_SO4 ) * nMolec / n_air;
                    /* Unit check:
                     * [1/s] * ( [molec/cm3] - [molec/cm^3] ) = [molec/cm^3/s]
                     * ([-] * [kg/s] + [kg/s]) * [#/cm^3] * [molec/kg] = [molec/cm^3/s]
                     * */

                    /* Liquid SO4 differential equation */
                    dxdt[EPM_ind_SO4l] = - entrainRate * ( x[EPM_ind_SO4l] - m_SO4l_mixingratio ) \
                                         + nucRate * x_SO4 * nMolec / n_air;

                    /* Gaseous SO4 differential equation */
                    dxdt[EPM_ind_SO4g] = - entrainRate * ( x[EPM_ind_SO4g] - m_SO4g_mixingratio ) \
                                         - sticking_SO4 * physFunc::thermalSpeed( x[EPM_ind_T], MW_H2SO4 / physConst::Na ) \
                                           * x[EPM_ind_Part] * n_air * 1.0E+06 * physConst::PI * m_part_r * m_part_r \
                                           * x[EPM_ind_SO4g] * ( 1.0 - x[EPM_ind_the1] - x[EPM_ind_the2] ) \
                                         - nucRate * x_SO4 * nMolec / n_air;
                    
                    /* On soot SO4 differential equation 
                     * \frac{d[SO4_s]}{dt} = alpha * v_th / 4.0 * n_part * 4.0 * \pi * r^2 * ( 1.0 - \theta ) * [SO4] */ 
                    dxdt[EPM_ind_SO4s] = - entrainRate * ( x[EPM_ind_SO4s] - 0.0 ) \
                                         + sticking_SO4 * physFunc::thermalSpeed( x[EPM_ind_T], MW_H2SO4 / physConst::Na ) \
                                           * x[EPM_ind_Part] * n_air * 1.0E+06 * physConst::PI * m_part_r * m_part_r \
                                           * x[EPM_ind_SO4g] * ( 1.0 - x[EPM_ind_the1] - x[EPM_ind_the2] );

                    /* Total SO4 differential equation
                     * \frac{d[SO4]}{dt} = - w_T(t) * ( [SO4] - [SO4]_amb ) */
                    dxdt[EPM_ind_SO4]  = dxdt[EPM_ind_SO4g] + dxdt[EPM_ind_SO4l] + dxdt[EPM_ind_SO4s];

                    /* Gaseous HNO3 differential equation 
                     * \frac{d[HNO3]}{dt} = - w_T(t) * ( [HNO3] - [HNO3]_amb ) */
                    dxdt[EPM_ind_HNO3] = - entrainRate * ( x[EPM_ind_HNO3] - m_HNO3_mixingratio );

                    /* Particle differential equation 
                     * \frac{d[part]}{dt} = - w_T(t) * ( [part] - [part]_amb )*/
                    dxdt[EPM_ind_Part] = - entrainRate * ( x[EPM_ind_Part] - m_part_mixingratio );

                    /* Particle radius differential equation
                     * \frac{dr}{dt} = \frac{dm}{dt}/(\rho * 4 * \pi * r^2) */
                    dxdt[EPM_ind_ParR] = ( isFreezable( m_part_r, x[EPM_ind_T], n_air*x[EPM_ind_H2O], m_part_r0 ) \
                                           * depositionRate( m_part_r, x[EPM_ind_T], x[EPM_ind_P], x[EPM_ind_H2O] * n_air , m_part_r0, x[EPM_ind_the1] + x[EPM_ind_the2] ) \
                                         + condensationRate( m_part_r, x[EPM_ind_T], x[EPM_ind_P], x[EPM_ind_H2O] * n_air , x[EPM_ind_the1] + x[EPM_ind_the2] ) ) \
                                         / ( physConst::RHO_ICE * 4.0 * physConst::PI * m_part_r * m_part_r );
                    /* Unit check: 
                     * [kg/s] / ( [kg/m^3] * [m^2]) = [m/s] */

                    /*Diff Eqs for theta1 (fractional soot coverage due to adsorption) and theta2 (due to scavenging)
                    CoagRate and d(theta1)/dt and d(theta2)/dt from Karcher (1998)*/
                    dxdt[EPM_ind_the1] = sticking_SO4 * physFunc::thermalSpeed( x[EPM_ind_T], MW_H2SO4 / physConst::Na ) * 1.0E+02 * 0.25 * n_air * ( SO4_g + SO4_rl ) / sigma_SO4 * ( 1.0 - x[EPM_ind_the1] - x[EPM_ind_the2] );

                    dxdt[EPM_ind_the2] = CoagRate / ( 4.0 * physConst::PI * m_part_r0 * m_part_r0 );

                } /* End of gas_aerosol_rhs::operator() */

        };


        Vector_1D x(13, 0.0);

        /* Initial conditions */
        x[EPM_ind_Trac] = 1.0;
        x[EPM_ind_T]    = Tc0;
        x[EPM_ind_P]    = pressure_Pa;
        x[EPM_ind_H2O]  = varArray[ind_H2O];
        x[EPM_ind_SO4]  = varArray[ind_SO4];

        /* Assume that emitted SO4 is purely gaseous as temperature is high */
        x[EPM_ind_SO4l] = 0.0;
        x[EPM_ind_SO4g] = varArray[ind_SO4];
        x[EPM_ind_SO4s] = 0.0;
        x[EPM_ind_HNO3] = varArray[ind_HNO3] / n_air_eng;
        x[EPM_ind_Part] = varSoot;
        x[EPM_ind_ParR] = EI.getSootRad();
        x[11] = 0.0;
        x[12] = 0.0; 

        Vector_2D obs_Var;
        Vector_1D obs_Time;

        EPM::streamingObserver observer( obs_Var, obs_Time, EPM_ind, micro_data_out, 2 );

        /* Creating ode's right hand side */
        gas_aerosol_rhs rhs( temperature_K, pressure_Pa, delta_T, H2O_amb, SO4_amb, SO4l_amb, SO4g_amb, HNO3_amb, Soot_amb, EI.getSootRad(), KernelSO4Soot);

        /* Total number of steps */
        UInt totSteps = 0;

        /* Current time step in s */
        double currTimeStep;

        /* Dilution factor */
        double dilFactor, dilFactor_b;

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
        AIM::Aerosol pSO4pdf_3mins( nPDF_SO4 );

        iTime = 0;
        while ( iTime + 1 < nTime ) {
   
            /* Update current time step */
            currTimeStep = timeArray[iTime+1] - timeArray[iTime];
            
            /* Coagulation of sulfate aerosols */
            nPDF_SO4.Coagulate( currTimeStep, Kernel );

            /* Get dilution factor */
            if ( iTime == 0 ) {
                dilFactor_b = 1.0;
                SO4l_b = 0.0;
            }
            else {
                dilFactor_b = observer.m_states[observer.m_states.size()-1][EPM_ind_Trac];
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
            
            dilFactor = observer.m_states[observer.m_states.size()-1][EPM_ind_Trac] / dilFactor_b;
            SO4l = observer.m_states[observer.m_states.size()-1][EPM_ind_SO4l] ;

            n_air = physConst::Na * x[EPM_ind_P]/(physConst::R * x[EPM_ind_T] * 1.0e6);
            n_air_prev = physConst::Na *  P_b/(physConst::R * T_b * 1.0e6);

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
                H2OMol_3mins   = observer.m_states[observer.m_states.size()-1][EPM_ind_H2O]; //* observer.m_states[observer.m_states.size()-1][EPM_ind_P] / ( physConst::kB * observer.m_states[observer.m_states.size()-1][EPM_ind_T] ) * 1.0E-06;
                Tracer_3mins   = observer.m_states[observer.m_states.size()-1][EPM_ind_Trac];
//                SO4pdf_3mins   = nPDF_SO4;
                SO4l_3mins     = observer.m_states[observer.m_states.size()-1][EPM_ind_SO4l]; // * observer.m_states[observer.m_states.size()-1][EPM_ind_P] / ( physConst::kB * observer.m_states[observer.m_states.size()-1][EPM_ind_T] * 1.0E+06 ) ;
                SO4g_3mins     = observer.m_states[observer.m_states.size()-1][EPM_ind_SO4g]; // * observer.m_states[observer.m_states.size()-1][EPM_ind_P] / ( physConst::kB * observer.m_states[observer.m_states.size()-1][EPM_ind_T] * 1.0E+06 ) ;
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
        if ( !CHEMISTRY && !observer.checkwatersat() ) {
            std::cout << "EndSim: Never reaches water saturation... ending simulation" << std::endl;
            //exit(0);
            return EPM_EARLY;
        }

        /* Persistent contrail */
        if ( relHumidity_i_Final >= 1.0 ) {

             Ice_rad  = PartRad_3mins;
             Ice_den  = PartDens_3mins;
             Soot_den = 0.0;
             H2O_mol  = physFunc::pSat_H2Os( final_Temp ); // / ( physConst::kB * final_Temp * 1.0E+06 )  ;

        } 
        /* No persistent contrail */
        else {

             Ice_rad  = EI.getSootRad();
             Ice_den  = 0.0;
             Soot_den = PartDens_3mins;
             H2O_mol  = H2OMol_3mins;
	         std::cout << "No persistent contrail..." << std::endl;
             if ( !CHEMISTRY ) return EPM_EARLY;

        }
	    std::cout << "Ice_den=" << Ice_den << std::endl;

        SO4l_mol = SO4l_3mins;
        SO4g_mol = SO4g_3mins;
        SO4Aer = pSO4pdf_3mins;

        const double expsIce = 1.15;
        AIM::Aerosol solidAer( Ice_rJ, Ice_rE, Ice_den, std::max(Ice_rad * exp ( -2.5 * log(expsIce) * log(expsIce) ), 1.5 * PA_R_LOW ), expsIce, "lognormal" );
        IceAer = solidAer;

        /* Compute plume area */
        Area = Ab0 / Tracer_3mins; 
    
//        std::cout << " Soot end: " << ( varSoot - Soot_amb  ) * 1E6 * Ab0 << " [#/m]\n";
//        std::cout << " Soot end: " << ( Soot_den - Soot_amb ) * 1E6 * Ab0 / Tracer_3mins << "[#/m]\n";

        /* Update temperature after vortex sinking */
        temperature_K = final_Temp;

        return EPM_SUCCESS;

    } /* End of RunMicrophysics */

    double dT_Vortex( const double time, const double delta_T, bool deriv )
    {

        /* DESCRIPTION:
         * Returns the temperature increase / rate of change from cruise altitude as a
         * function of time */

        /* INPUT:
         * - double time    :: time since engine exit plane expressed in s
         * - double delta_T :: temperature increase in K after vortex sinking in K
         * - bool deriv         :: rate of change or delta?
         *
         * OUTPUT:
         * - double :: temperature increase / rate of change from cruise altitude in K */

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
    
    double entrainmentRate( const double t )
    {

        /* DESCRIPTION:
         * Returns the plume entrainment rate in the plume early stages.
         * Based on Karcher 1995 */

        /* INPUT:
         * - double t :: time since engine exit plane expressed in s
         *
         * OUTPUT:
         * - double :: plume entrainment rate 1/s */

        double dRat = 0;

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

    } /* End of entrainmentRate */

    double depositionRate( const double r, const double T, const double P, const double H2O, const double r_0,  const double theta )
    {

        /* DESCRIPTION:
         * Returns the water deposition rate in kg/s on a single spherical particle */

        /* INPUT:
         * - double r     :: radius in m
         * - double T     :: temperature in K
         * - double P     :: pressure in Pa
         * - double H2O   :: gaseous water concentration in molec/cm^3
         * - double r_0   :: bare soot radius in m
         * - double theta :: particle surface coverage in [-]
         *
         * OUTPUT:
         * - double :: deposition rate */

        /* Capacitance factor for spherical particles is C = r; */

        double p_h2o = H2O * physConst::kB * T * 1.0E+06; // Partial pressure of water

        // If supersaturated w.r.t. liquid 
        if ( ( p_h2o / physFunc::pSat_H2Ol( T ) >= 1.0 ) && ( p_h2o / physFunc::pSat_H2Ol( T ) < 2.0 ) ) {
            if ( ( r * r - r_0 * r_0 ) / ( r_0 * r_0 ) <= 0.5 ) {
                return 4.0 * physConst::PI * r * theta * physFunc::CorrDiffCoef_H2O( r, T, P ) * \
                       ( p_h2o - physFunc::pSat_H2Os( T ) * physFunc::Kelvin( r ) ) / \
                       ( physConst::R / MW_H2O * T + physFunc::LHeatSubl_H2O( T ) * physFunc::CorrDiffCoef_H2O( r, T, P ) * \
                         physFunc::pSat_H2Os( T ) * ( physFunc::LHeatSubl_H2O( T ) * MW_H2O / ( physConst::R * T ) - 1.0 ) / \
                           ( physFunc::ThermalCond( r, T, P ) * T ) );
            } else {
                return 4.0 * physConst::PI * r * physFunc::CorrDiffCoef_H2O( r, T, P ) * \
                       ( p_h2o - physFunc::pSat_H2Os( T ) * physFunc::Kelvin( r ) ) / \
                       ( physConst::R / MW_H2O * T + physFunc::LHeatSubl_H2O( T ) * physFunc::CorrDiffCoef_H2O( r, T, P ) * \
                         physFunc::pSat_H2Os( T ) * ( physFunc::LHeatSubl_H2O( T ) * MW_H2O / ( physConst::R * T ) - 1.0 ) / \
                           ( physFunc::ThermalCond( r, T, P ) * T ) );
            }
        } else { // Melting
            return 4.0 * physConst::PI * r * physFunc::CorrDiffCoef_H2O( r, T, P ) * \
                   ( p_h2o  - physFunc::pSat_H2Os( T ) * physFunc::Kelvin( r ) ) / \
                   ( physConst::R / MW_H2O * T + physFunc::LHeatSubl_H2O( T ) * physFunc::CorrDiffCoef_H2O( r, T, P ) * \
                     physFunc::pSat_H2Os( T ) * ( physFunc::LHeatSubl_H2O( T ) * MW_H2O / ( physConst::R * T ) - 1.0 ) / \
                       ( physFunc::ThermalCond( r, T, P ) * T ) );
        }

    } /* End of depositionRate */

    bool isFreezable( const double r, const double T, const double H2O, const double r0 )
    {

        /* DESCRIPTION:
         * Returns whether or not water can freeze on a particle */

        /* INPUT:
         * - double r   :: radius in m
         * - double T   :: temperature in K
         * - double P   :: pressure in Pa
         * - double H2O :: gaseous water concentration in molec/cm^3
         * - double r0  :: bare soot radius in m
         *
         * OUTPUT:
         * - bool :: isFreezable? */

        /* Returns 1 if: 
         * - Current conditions are supersaturated with respect to liquid water (freezing)
         * or 
         * - Current radius is greater than bare soot radius (melting) */

        return ( ( H2O * physConst::kB * T * 1.0E+06 >= physFunc::Kelvin( r ) * physFunc::pSat_H2Ol( T ) ) || ( r > 1.005 * r0 ) );

    } /* End of isFreezable */

    double condensationRate( const double r, const double T, const double P, const double H2O, const double theta )
    {

        /* DESCRIPTION:
         * Returns the gaseous water condensation rate in kg/s on a single spherical particle */

        /* INPUT:
         * - double r     :: radius in m
         * - double T     :: temperature in K
         * - double P     :: pressure in Pa
         * - double H2O   :: gaseous water concentration in molec/cm^3
         * - double theta :: fractional coating of the particle
         *
         * OUTPUT:
         * - double :: condensation rate */
        
        double p_h2o = H2O * physConst::kB * T * 1.0E+06; // Partial pressure of water

        if ( p_h2o >= physFunc::Kelvin( r ) * physFunc::pSat_H2Ol( T ) ) {
            return 4.0 * physConst::PI * MW_H2O * physFunc::CorrDiffCoef_H2O( r, T, P ) * r * theta \
                       * ( H2O / physConst::Na * 1.0E+06 \
                         - physFunc::Kelvin( r ) * physFunc::pSat_H2Ol( T ) / ( physConst::R * T ) );
        } else {
            return 0.0;
        }

    } /* End of condensationRate */

}

/* End of Integrate.cpp */
