/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* LiquidAer Header File                                            */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 10/22/2018                                */
/* File                 : LiquidAer.hpp                             */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef LIQUIDAER_H_INCLUDED
#define LIQUIDAER_H_INCLUDED

#include <vector>

double H2SO4_GASFRAC( const double temperature_K, const double SO4 );
unsigned int STRAT_AER( const double temperature_K     , const double pressure_Pa       , const double airDens ,               \
                        const double latitude_deg      , std::vector<double> &Data      , double Area ,                        \
                        std::vector<double> &KHETI_SLA , std::vector<double> &SOLIDFRAC ,                                      \
                        std::vector<double> &AERFRAC   , std::vector<double> &RAD_AER   ,                                      \
                        std::vector<double> &RHO_AER   , std::vector<double> &KG_AER    ,                                      \
                        std::vector<double> &NDENS_AER , std::vector<double> &SAD_AER   , 
                        double tropopausePressure, bool DBG_IN );
std::vector<double> SLA_GAMMA( const double T_K         , const double P_Pa     , \
                               const double WT_FRC      ,                         \
                               const double H2OSUM      , const double HClSUM   , \
                               const double HBrSUM      , const double HOBrSUM  , \
                               const double ClNO3SUM    , const double BrNO3SUM , \
                               const double RHO         , const double ARAD );
void TERNARY( const double TIN_K    , const double PIN_Pa  , const double H2OSUM_IN , \
              const double H2SO4SUM , const double HNO3SUM , const double HClSUM    , \
              const double HOClSUM  , const double HBrSUM  , const double HOBrSUM   , \
              double &W_H2SO4       , double &W_H2O        , double &W_HNO3         , \
              double &W_HCl         , double &W_HOCl       , double &W_HBr          , \
              double &W_HOBr        , double &HNO3GASFRAC  , double &HClGASFRAC     , \
              double &HOClGASFRAC   , double &HBrGASFRAC   , double &HOBrGASFRAC    , \
              double &SLA_VOL       , double &SLA_RHO );
double CARSLAW_DENSITY( const double CS, const double CN, const double T_K );


#endif /* LIQUIDAER_H_INCLUDED */
