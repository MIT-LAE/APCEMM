/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* PhysConstant Header File                                         */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : PhysConstant.hpp                          */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */


#ifndef PHYSCONSTANT_H_INCLUDED
#define PHYSCONSTANT_H_INCLUDED

#include "ForwardsDecl.hpp"

#include "MolarWeights.hpp"

/* Physical constants */
#ifndef PHYS_CONSTANTS
#define PHYS_CONSTANTS

namespace physConst
{

    /* Defined physical parameters */


    /* ABSOLUTE PHYSICAL CONSTANTS */

    /* Double-precision value of \pi */
    static const RealDouble PI = 3.141592653589793238460;

    /* Avogadro number   , Unit : [ molecules / mol ] */
    static const RealDouble Na = 6.022140857E+23;  

    /* Boltzmann constant, Unit : [ J /  K ] */
    static const RealDouble kB = 1.380648528E-23;

    /* Ideal gas constant, Unit : [ J / ( K.mol ) ] */
    static const RealDouble R = 8.314459848E+00;

    /* Specific gas constant, Unit : [ J / (K.kg) ] */
    static const RealDouble R_Air = R / MW_Air;

    
    /* THERMODYNAMIC QUANTITIES */

    /* Specific heat capacity of air at 300K, Unit : [ kJ / ( kg K ) ] */
    static const RealDouble CP_Air = 1.005E+00;

    /* Heat capacity ratio of air, Unit : [ - ] */
    static const RealDouble GAMMA_Air = 1.4;


    /* MOLECULAR PROPERTIES */

    /* Mass of one air molecule, Unit : [ kg ] */
    static const RealDouble M_Air = MW_Air / Na;

    /* Density of water, Unit : [ kg / m ^ 3 ] */
    static const RealDouble RHO_H2O = 1.000E+03;

    /* Density of ice, Unit : [ kg / m ^ 3 ] */
    static const RealDouble RHO_ICE = 9.167E+02;

    /* Density of sulfuric acid, Unit : [ kg / m ^ 3 ] */
    static const RealDouble RHO_SO4 = 1.600E+03;
    
    /* Density of soot, Unit : [ kg / m ^ 3 ] */
    static const RealDouble RHO_SOOT = 2.000E+03;


    /* EARTH'S PHYSICAL PARAMETERS */

    /* Acceleration due to gravity at the surface, Unit : [ m / s^2 ] */
    static const RealDouble g = 9.80665E+00;

    /* Pressure at sea level, Unit : [ Pa ] */
    static const RealDouble PRES_SL = 101325.0;

    /* 1 Atmosphere, Unit : [ Pa ] */
    static const RealDouble ATM = PRES_SL;

    /* Temperature at sea level, Unit : [ K ] */
    static const RealDouble TEMP_SL = 288.15;

    /* Reference temperature, Unit : [ K ] */
    static const RealDouble TEMP_REF = 273.15;


    /* ATMOSPHERIC PARAMETERS */

    /* Adiabatic temperature lapse rate, Unit : [ K / km ] */
    static const RealDouble GAMMA_AD = -g / CP_Air;

    /* Turbulent dissipation rate, 
     * Range : ( 1.0E-08 - 1.0E-02 ), 
     * Unit  : [ m ^ 2 / s ^ 3 ]
     * Source: U. Schumann, A contrail cirrus prediction model, Geo. Sc. Dev., 2012 */
    static const RealDouble EPSILON = 1.00E-05;

}

#endif /* PHYS_CONSTANTS */

#endif /* PHYSCONSTANT_H_INCLUDED */
