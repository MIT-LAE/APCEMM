/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Structure Header File                                            */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : Structure.hpp                             */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef STRUCTURE_H_INCLUDED
#define STRUCTURE_H_INCLUDED

#include <cmath>
#include "AIM/Aerosol.hpp"
#include "KPP/KPP_Global.h"
#include "Util/ForwardDecl.hpp"
#include "Core/Input.hpp"
#include "Core/Input_Mod.hpp"
#include "Core/Emission.hpp"
#include "Core/Aircraft.hpp"
#include "Core/Input_Mod.hpp"
#include "Core/Meteorology.hpp"

class Solution
{
public:

    Solution(const OptInput& optInput);

    ~Solution();

    // NOTE: THIS IS ONLY CALLED FOR THE EPM
    void Initialize(std::string fileName, const Input &input,
                    const double airDens, const Meteorology &met,
                    const OptInput &Input_Opt, double *varSpeciesArray,
                    double *fixSpeciesArray, const bool DBG);

    void getData(double *varSpeciesArray, double *fixSpeciesArray,
                 const UInt i = 0, const UInt j = 0, const bool CHEMISTRY = 0);

    void SpinUp(Vector_1D &amb_Value, const Input &input,
                const double airDens, const double startTime,
                double *varSpeciesArray, double *fixSpeciesArray,
                const bool DGB = 0);

    Vector_2D getAerosol( ) const;

    UInt Nx() const { return size_x; };
    UInt Ny() const { return size_y; };

    /* Species */
    Vector_3D Species;

    /* Aerosols */
    Vector_2D sootDens, sootRadi, sootArea;

    AIM::Grid_Aerosol liquidAerosol, solidAerosol;

    UInt nBin_PA;
    UInt nBin_LA;

    double LA_nDens, LA_rEff, LA_SAD;
    double PA_nDens, PA_rEff, PA_SAD;

    Vector_1D KHETI_SLA;
    Vector_1D AERFRAC, SOLIDFRAC;
    UInt STATE_PSC;

private:
    void readInputBackgroundConditions(const Input &input, Vector_1D &amb_Value,
                                       Vector_2D &aer_Value,
                                       std::string filename);
    void processInputBackgroundLine(std::istream &s, Vector_1D &amb_Value,
                                    Vector_2D &aer_Value);
    void setAmbientConcentrations(const Input &input, Vector_1D &amb_Value);
    void initializeSpeciesH2O(const Input &input, const OptInput &input_Opt,
                              Vector_1D &amb_Value, const double airDens,
                              const Meteorology &met);
    void setSpeciesValues(Vector_1D &AERFRAC, Vector_1D &SOLIDFRAC,
                          const Vector_1D &stratData);

  const UInt nVariables;
  const UInt nAer;
  const UInt size_x;
  const UInt size_y;

  bool reducedSize;
};

#endif /* STRUCTURE_H_INCLUDED */
