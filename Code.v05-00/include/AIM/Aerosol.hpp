/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*                        AIrcraft Microphysics                     */
/*                              (AIM)                               */
/*                                                                  */
/* Aerosol Header File                                              */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 9/28/2018                                 */
/* File                 : Aerosol.hpp                               */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef AEROSOL_H_INCLUDED
#define AEROSOL_H_INCLUDED

#include <iostream>
#include <cmath>
#include <vector>
#include <cstring>
#include <boost/math/special_functions/gamma.hpp>
#ifdef OMP
    #include "omp.h"
#endif /* OMP */

/* Include Parameters.hpp for multithreading option */
#include "Core/Parameters.hpp"

#include "Util/ForwardDecl.hpp"
#include "Util/PhysConstant.hpp"
#include "AIM/Coagulation.hpp"
#include "Core/Mesh.hpp"

namespace AIM
{
    class Aerosol;
    class Grid_Aerosol;
    static const RealDouble DEFAULT_MIN_RADIUS = 1.0e-8; 
    static const RealDouble TINY = 1.0e-50;
}

class AIM::Aerosol
{

    public:

        /* Default constructor */
        Aerosol( );

        /* Constructor */
        Aerosol( Vector_1D bin_Centers, Vector_1D bin_Edges, RealDouble nPart, RealDouble mu, RealDouble sigma, const char* distType = "lognormal", RealDouble alpha_ = -1.0, RealDouble gamma_ = -1.0, RealDouble b_ = 0.0 );

        /* Copy */
        Aerosol( const Aerosol &rhs );

        /* Destructor */
        ~Aerosol();

        /* Operators */
        Aerosol& operator=( const Aerosol &rhs );
        Aerosol operator+=( const Aerosol &rhs );
        Aerosol operator-=( const Aerosol &rhs );
        Aerosol operator+( const Aerosol &rhs ) const;
        Aerosol operator-( const Aerosol &rhs ) const;

        /* Coagulation */
        void Coagulate( const RealDouble dt, const Coagulation &kernel );
        void Coagulate( const RealDouble dt, const Vector_2D &beta, const Vector_3D &f );
        
        /* Moments */
        RealDouble Moment( UInt n = 0 ) const;
        RealDouble Radius( ) const;
        RealDouble EffRadius( ) const;
        RealDouble StdDev( ) const;

        /* utils */
        void scalePdf( RealDouble factor );
        void updatePdf( Vector_1D pdf_ );

        /* gets */
        Vector_1D getBinCenters() const;
        Vector_1D getBinVCenters() const;
        Vector_1D getBinEdges() const;
        Vector_1D getBinSizes() const;
        UInt getNBin() const;
        RealDouble getNPart() const;
        const char* getType() const;
        RealDouble getAlpha() const;
        Vector_1D getPDF() const;

        static Vector_1D generateEndOfEPMTestPDF(const Vector_1D& bin_Centers);


    protected:

        Vector_1D bin_Centers;
        Vector_1D bin_VCenters;
        Vector_1D bin_Edges;
        Vector_1D bin_Sizes;
        UInt nBin;
        RealDouble nPart;
        const char* type;
        RealDouble mu;
        RealDouble sigma;
        RealDouble alpha;
        Vector_1D pdf;

    private:

};

class AIM::Grid_Aerosol
{

    public:

        /* Default constructor */
        Grid_Aerosol( );

        /* Constructor */
        Grid_Aerosol( UInt Nx_, UInt Ny_, Vector_1D bin_Centers, Vector_1D bin_Edges, RealDouble nPart, RealDouble mu, RealDouble sigma, const char* distType = "lognormal", RealDouble alpha_ = -1.0, RealDouble gamma_ = -1.0, RealDouble b_ = 0.0 );

        /* Copy */
        Grid_Aerosol( const Grid_Aerosol &rhs );

        /* Destructor */
        ~Grid_Aerosol();

        /* Operators */
        Grid_Aerosol& operator=( const Grid_Aerosol &rhs );
        Grid_Aerosol operator+=( const Grid_Aerosol &rhs );
        Grid_Aerosol operator-=( const Grid_Aerosol &rhs );
        Grid_Aerosol operator+( const Grid_Aerosol &rhs ) const;
        Grid_Aerosol operator-( const Grid_Aerosol &rhs ) const;

        /* Coagulation */
        void Coagulate( const RealDouble dt, Coagulation &kernel, const UInt N = 2, const UInt SYM = 0 );

        /* Ice crystal growth */
        void Grow( const RealDouble dt, Vector_2D &H2O, const Vector_2D &T, const Vector_1D &P, const UInt N = 2, const UInt SYM = 0 );
        RealDouble EffDiffCoef( const RealDouble r, const RealDouble T, const RealDouble P, const RealDouble H2O) const;
        void APC_Scheme(const UInt jNy, const UInt iNx, const double T, const double P,
                            const double dt, Vector_2D& H2O, Vector_2D& totH2O, Vector_3D& icePart, Vector_3D& iceVol);
        std::vector<int> ComputeBinParticleFlux(const int x_index, const int y_index, const Vector_3D& iceVol, const Vector_3D& icePart) const;
        void ApplyBinParticleFlux(const int x_index, const int y_index, const std::vector<int> &toBin, const Vector_3D &iceVol, const Vector_3D &icePart);
        
        /* Helper Functions for Coagulation and Ice Growth */
        bool CheckCoagAndGrowInputs(const UInt N, const UInt SYM, UInt& Nx_max, UInt& Ny_max, const std::string funcName) const;
        void CoagAndGrowApplySymmetry(const UInt N, const UInt SYM, const UInt Nx_max, const UInt Ny_max, const char* funcName, Vector_2D& H2O);
        /* Update bin centers - Used after aerosol transport */
        void UpdateCenters( const Vector_3D &iceV, const Vector_3D &PDF );

        /* Moments */
        Vector_2D Moment( UInt n ) const;
        RealDouble Moment( UInt n, Vector_1D PDF ) const;
        RealDouble Moment( UInt n, UInt iNx, UInt jNy ) const;

        /* Extra utils */
        Vector_3D Number( ) const;
        Vector_2D TotalNumber( ) const;
        RealDouble TotalNumber_sum( const Vector_2D cellAreas ) const;
        Vector_1D Overall_Size_Dist( const Vector_2D cellAreas ) const;
        Vector_3D Volume( ) const;
        Vector_2D TotalVolume( ) const;
        Vector_2D TotalArea( ) const;
        RealDouble TotalIceMass_sum( const Vector_2D cellAreas ) const;
        Vector_2D IWC( ) const;
        Vector_2D Extinction( ) const;
        Vector_1D PDF_Total( const Vector_2D &cellAreas ) const;
        Vector_1D PDF_Total( const Mesh &m ) const;
        Vector_1D xOD( const Vector_1D dx ) const;
        Vector_1D yOD( const Vector_1D dy ) const;
        Vector_2D Radius( ) const;
        RealDouble Radius( UInt iNx, UInt jNy ) const;
        Vector_2D EffRadius( ) const;
        RealDouble EffRadius( UInt iNx, UInt jNy ) const;
        Vector_2D StdDev( ) const;
        RealDouble StdDev( UInt iNx, UInt jNy ) const;

        /* utils */
        void updatePdf( Vector_3D pdf_ );
        Vector_1D Average( const Vector_2D &weights,   \
                           const RealDouble &totWeight ) const;
        void addPDF( const Aerosol &PDF, const Vector_2D &weights,  \
                     const Vector_2D &cellAreas, const RealDouble nCell );
        void addPDF( const Vector_1D &PDF, const Vector_2D &weights, \
                     const Vector_2D &cellAreas, const RealDouble nCell );

        /* gets */
        Vector_1D getBinCenters() const;
        Vector_1D binCenters() const { return bin_Centers; };
        Vector_3D getBinVCenters() const;
        Vector_1D getBinEdges() const;
        Vector_1D binEdges() const { return bin_Edges; };
        Vector_1D getBinSizes() const;
        Vector_1D binSizes() const { return bin_Sizes; };
        UInt getNBin() const;
        Vector_3D getPDF() const;

        Vector_3D pdf;
        Vector_3D bin_VCenters;

    protected:

        unsigned int Nx, Ny;
        Vector_1D bin_Centers;
        Vector_1D bin_Edges;
        Vector_1D bin_VEdges;
        Vector_1D bin_Sizes;
        UInt nBin;
        RealDouble nPart;
        const char* type;
        RealDouble mu;
        RealDouble sigma;
        RealDouble alpha;

    private:

};

#endif /* AEROSOL_H_INCLUDED */
