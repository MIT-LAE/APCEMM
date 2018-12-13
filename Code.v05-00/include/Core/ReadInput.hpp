/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* ReadInput Header File                                            */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 12/10/2018                                */
/* File                 : ReadInput.hpp                             */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef INPUT_READER_H_INCLUDED
#define INPUT_READER_H_INCLUDED

#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>
#include <cstdlib>
#include <cmath>

#include "Core/Input_Mod.hpp"
#include "Util/ForwardDecl.hpp"

void Read_Input_File( OptInput &Input_Opt );
std::vector<std::string> Split_Line( std::string line2split, const std::string delimiter );
void Read_Simulation_Menu( OptInput &Input_Opt, bool &RC );
void Read_Parameters( OptInput &Input_Opt, bool &RC );
void Read_Transport_Menu( OptInput &Input_Opt, bool &RC );
void Read_Chemistry_Menu( OptInput &Input_Opt, bool &RC );
void Read_Aerosol_Menu( OptInput &Input_Opt, bool &RC );

Vector_2D ReadParameters( const OptInput Input_Opt );
Vector_2D Copy_blocked( Vector_2D& m, int n );
Vector_2D Copy_interleaved( Vector_2D& m, int n );
Vector_2D Reshape_Vector( Vector_2D& vector_2D, int n_x, int n_y );
Vector_2D CombVec( const Vector_1D& temperature_K, \
                   const Vector_1D& pressure_Pa,   \
                   const Vector_1D& relHumidity_w, \
                   const Vector_1D& longitude_deg, \
                   const Vector_1D& latitude_deg,  \
                   const Vector_1D& dayGMT,        \
                   const Vector_1D& emissionTime,  \
                   const Vector_1D& EI_NOx,        \
                   const Vector_1D& EI_CO,         \
                   const Vector_1D& EI_HC,         \
                   const Vector_1D& EI_SO2,        \
                   const Vector_1D& EI_Soot,       \
                   const Vector_1D& SootRad,       \
                   const Vector_1D& ff,            \
                   const Vector_1D& backgNOx,      \
                   const Vector_1D& backgHNO3,     \
                   const Vector_1D& backgO3,       \
                   const Vector_1D& backgCO,       \
                   const Vector_1D& backgCH4,      \
                   const Vector_1D& backgSO2 );

#endif /* INPUT_READER_H_INCLUDED */

