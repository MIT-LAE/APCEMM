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
#include <vector>
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
void Read_Meteorology_Menu( OptInput &Input_Opt, bool &RC );
void Read_Diagnostic_Menu( OptInput &Input_Opt, bool &RC );
void Read_Timeseries_Menu( OptInput &Input_Opt, bool &RC );
void Read_PL_Menu( OptInput &Input_Opt, bool &RC );

Vector_2D Copy_blocked( Vector_2D& m, int n );
Vector_2D Copy_interleaved( Vector_2D& m, int n );
Vector_2D Reshape_Vector( Vector_2D& vector_2D, int n_x, int n_y );
Vector_2D CombVec( OptInput &Input_Opt );

#endif /* INPUT_READER_H_INCLUDED */

