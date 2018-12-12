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

void Read_Input_File( const char* fileName );
std::vector<std::string> Split_Line( std::string line2split, const std::string delimiter );
void Read_Simulation_Menu( Option &Input_Opt, bool RC );

#endif /* INPUT_READER_H_INCLUDED */

