/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Input Header File                                                */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 11/19/2018                                */
/* File                 : Input.hpp                                 */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef INPUT_H_INCLUDED
#define INPUT_H_INCLUDED

#include "Util/ForwardDecl.hpp"

/* Input is the scenario a run simulates (aircraft, emissions and atmospheric
 * parameters).
 * Read from the PARAMETER MENU by YamlInputReader::readParamMenu
 * which fills the Input object through setters.
 * Everything that decides how the model runs is in OptInput. */
class Input
{

    double simulationTime_ = 0.0;

    double temperature_K_ = 0.0;
    double relHumidity_w_ = 0.0;
    double relHumidity_i_ = 0.0;
    double horizDiff_ = 0.0;
    double vertiDiff_ = 0.0;

    double longitude_deg_ = 0.0;
    double latitude_deg_ = 0.0;
    double pressure_Pa_ = 0.0;

    UInt emissionDOY_ = 0;
    double emissionTime_ = 0.0;

    double EI_SO2_ = 0.0;
    double EI_SO2TOSO4_ = 0.0;
    double EI_Soot_ = 0.0;
    double sootRad_ = 0.0;

    double fuelFlow_ = 0.0;
    double aircraftMass_ = 0.0;

    double backgNOx_ = 0.0;
    double backgHNO3_ = 0.0;
    double backgO3_ = 0.0;
    double backgCO_ = 0.0;
    double backgCH4_ = 0.0;
    double backgSO2_ = 0.0;

    double flightSpeed_ = 0.0;
    double numEngines_ = 0.0;
    double wingspan_ = 0.0;
    double coreExitTemp_ = 0.0;
    double bypassArea_ = 0.0;

    double nBV_ = 0.0;

    UInt emissionDay_ = 0;
    UInt emissionMonth_ = 0;

    private:
        void adjustLatLong();
        void calculate_emissionMonth();

    public:

        void checkInputValidity() const;

        double simulationTime() const { return simulationTime_; }
        double temperature_K() const { return temperature_K_; } // From the meteorology
        double pressure_Pa() const { return pressure_Pa_; }
        double relHumidity_w() const { return relHumidity_w_; } // From the meteorology
        double relHumidity_i() const { return relHumidity_i_; } // From the meteorology
        double horizDiff() const { return horizDiff_; }
        double vertiDiff() const { return vertiDiff_; }
        inline double nBV() const { return nBV_; }

        void set_temperature_K( double T_CA_K ) { temperature_K_ = T_CA_K; }
        void set_relHumidity_w( double RHW_CA ) { relHumidity_w_ = RHW_CA; }
        void set_relHumidity_i( double RHi_CA ) { relHumidity_i_ = RHi_CA; }

        double longitude_deg() const { return longitude_deg_; }
        double latitude_deg() const { return latitude_deg_; }

        UInt emissionDOY() const { return emissionDOY_; }
        UInt emissionDay() const { return emissionDay_; }
        UInt emissionMonth() const { return emissionMonth_; }
        double emissionTime() const { return emissionTime_; }

        double EI_SO2() const { return EI_SO2_; }
        double EI_SO2TOSO4() const { return EI_SO2TOSO4_; }
        double EI_Soot() const { return EI_Soot_; }
        double sootRad() const { return sootRad_; }

        double fuelFlow() const { return fuelFlow_; }

        double aircraftMass() const { return aircraftMass_; }
        double flightSpeed() const { return flightSpeed_; }
        double numEngines() const { return numEngines_; }
        double wingspan() const { return wingspan_; }
        double coreExitTemp() const { return coreExitTemp_; }
        double bypassArea() const { return bypassArea_; }

        double backgNOx() const { return backgNOx_; }
        double backgHNO3() const { return backgHNO3_; }
        double backgO3() const { return backgO3_; }
        double backgCO() const { return backgCO_; }
        double backgCH4() const { return backgCH4_; }
        double backgSO2() const { return backgSO2_; }

        /* Setters used by the input reader. The meteorology setters above
         * overwrite temperature and humidity later in the run. */

        void set_simulationTime( double simulationTime ) { simulationTime_ = simulationTime; }
        void set_horizDiff( double horizDiff ) { horizDiff_ = horizDiff; }
        void set_vertiDiff( double vertiDiff ) { vertiDiff_ = vertiDiff; }
        void set_nBV( double nBV ) { nBV_ = nBV; }

        void set_longitude_deg( double longitude_deg ) { longitude_deg_ = longitude_deg; }
        void set_latitude_deg( double latitude_deg ) { latitude_deg_ = latitude_deg; }
        void set_pressure_Pa( double pressure_Pa ) { pressure_Pa_ = pressure_Pa; }

        void set_emissionDOY( UInt emissionDOY ) { emissionDOY_ = emissionDOY; }
        void set_emissionTime( double emissionTime ) { emissionTime_ = emissionTime; }

        void set_EI_SO2( double EI_SO2 ) { EI_SO2_ = EI_SO2; }
        void set_EI_SO2TOSO4( double EI_SO2TOSO4 ) { EI_SO2TOSO4_ = EI_SO2TOSO4; }
        void set_EI_Soot( double EI_Soot ) { EI_Soot_ = EI_Soot; }
        void set_sootRad( double sootRad ) { sootRad_ = sootRad; }

        void set_fuelFlow( double fuelFlow ) { fuelFlow_ = fuelFlow; }
        void set_aircraftMass( double aircraftMass ) { aircraftMass_ = aircraftMass; }
        void set_flightSpeed( double flightSpeed ) { flightSpeed_ = flightSpeed; }
        void set_numEngines( double numEngines ) { numEngines_ = numEngines; }
        void set_wingspan( double wingspan ) { wingspan_ = wingspan; }
        void set_coreExitTemp( double coreExitTemp ) { coreExitTemp_ = coreExitTemp; }
        void set_bypassArea( double bypassArea ) { bypassArea_ = bypassArea; }

        void set_backgNOx( double backgNOx ) { backgNOx_ = backgNOx; }
        void set_backgHNO3( double backgHNO3 ) { backgHNO3_ = backgHNO3; }
        void set_backgO3( double backgO3 ) { backgO3_ = backgO3; }
        void set_backgCO( double backgCO ) { backgCO_ = backgCO; }
        void set_backgCH4( double backgCH4 ) { backgCH4_ = backgCH4; }
        void set_backgSO2( double backgSO2 ) { backgSO2_ = backgSO2; }

};

#endif /* INPUT_H_INCLUDED */
