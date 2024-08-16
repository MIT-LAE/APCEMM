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

#include <cmath>
#include <vector>
#include <cstring>
#include <boost/math/special_functions/gamma.hpp>
#include "APCEMM.h"
#include <Util/VectorUtils.hpp>
#include <type_traits>
#ifdef OMP
    #include "omp.h"
#endif /* OMP */

#include "Util/ForwardDecl.hpp"
#include "AIM/Coagulation.hpp"
#include "Core/Mesh.hpp"

namespace AIM
{
    class Aerosol;
    class Grid_Aerosol;
    static const double DEFAULT_MIN_RADIUS = 1.0e-8; 
    static const double TINY = 1.0e-50;
}

class AIM::Aerosol
{

    public:
        
        /* Default constructor (Needed because EPM is terribly engineered. Ideally should be deleted.) */
        Aerosol( ) = default;

        /* Constructor */
        Aerosol( const Vector_1D& bin_Centers, const Vector_1D& bin_Edges, double nPart, double mu, double sigma, const char* distType = "lognormal", double alpha_ = -1.0, double gamma_ = -1.0, double b_ = 0.0 );

        /* Operators */
        void addAerosolToPDF( const Aerosol &rhs );

        /* Coagulation */
        void Coagulate( const double dt, const Coagulation &kernel );
        void Coagulate( const double dt, const Vector_2D &beta, const Vector_3D &f );
        
        /* Moments */
        //NOTE: This gives the moment in [- / cm3]. You need to multiply by factors to get it in [ / m] or something.
        inline double binMoment(int iBin, int n = 0) const {
            return (log(bin_Edges[iBin + 1]) - log(bin_Edges[iBin])) * pow(bin_Centers[iBin], n) * pdf[iBin];
        }
        //NOTE: This gives the moment in [- / cm3]. You need to multiply by factors to get it in [ / m] or something.
        double Moment( UInt n = 0 ) const;
        double Radius( ) const;
        double EffRadius( ) const;
        double StdDev( ) const;

        /* utils */
        void scalePdf( double factor );
        void updatePdf( const Vector_1D& pdf_ );

        /* gets */
        inline const Vector_1D& getBinCenters() const { return bin_Centers; };
        inline const Vector_1D& getBinVCenters() const { return bin_VCenters; };
        inline const Vector_1D& getBinEdges() const { return bin_Edges; };
        inline const Vector_1D& getBinSizes() const { return bin_Sizes; };
        inline UInt getNBin() const { return nBin; };
        inline const char* getType() const { return type; };
        inline double getAlpha() const { return alpha; };
        inline const Vector_1D& getPDF() const { return pdf; };

    protected:

        Vector_1D bin_Centers;
        Vector_1D bin_VCenters;
        Vector_1D bin_Edges;
        Vector_1D bin_Sizes;
        UInt nBin;
        const char* type;
        double mu;
        double sigma;
        double alpha;
        Vector_1D pdf; // represents dN [-/cm3] / d(ln r)

    private:

};

class AIM::Grid_Aerosol
{

    public:

        /* Default constructor (Needed because EPM is terribly engineered. Ideally should be deleted.)*/
        Grid_Aerosol( ) = default;

        /* Constructor */
        Grid_Aerosol( UInt Nx_, UInt Ny_, const Vector_1D& bin_Centers, const Vector_1D& bin_Edges, double nPart, double mu, double sigma, const char* distType = "lognormal", double alpha_ = -1.0, double gamma_ = -1.0, double b_ = 0.0 );

        /* Coagulation */
        void Coagulate( const double dt, Coagulation &kernel, const UInt N = 2, const UInt SYM = 0 );

        /* Ice crystal growth */
        void Grow( const double dt, Vector_2D &H2O, const Vector_2D &T, const Vector_1D &P, const UInt N = 2, const UInt SYM = 0 );
        double EffDiffCoef( const double r, const double T, const double P, const double H2O) const;
        void APC_Scheme(const UInt jNy, const UInt iNx, const double T, const double P,
                            const double dt, Vector_2D& H2O, Vector_2D& totH2O, Vector_3D& icePart, Vector_3D& iceVol);
        std::vector<int> ComputeBinParticleFlux(const int x_index, const int y_index, const Vector_3D& iceVol, const Vector_3D& icePart) const;
        void ApplyBinParticleFlux(const int x_index, const int y_index, const std::vector<int> &toBin, const Vector_3D &iceVol, const Vector_3D &icePart);
        
        /* Helper Functions for Coagulation and Ice Growth */
        bool CheckCoagAndGrowInputs(const UInt N, const UInt SYM, UInt& Nx_max, UInt& Ny_max, const std::string funcName) const;
        void CoagAndGrowApplySymmetry(const UInt N, const UInt SYM, const UInt Nx_max, const UInt Ny_max, const char* funcName, Vector_2D& H2O);
        /* Update bin centers - Used after aerosol transport */
        void UpdateCenters( const Vector_3D &iceV, const Vector_3D &PDF );
        inline void updateNx(int nx_new) { Nx = nx_new; };
        inline void updateNy(int ny_new) { Ny = ny_new; };

        /* Moments */
        Vector_2D Moment( UInt n ) const;
        double Moment( UInt n, const Vector_1D& PDF ) const;
        double Moment( UInt n, UInt iNx, UInt jNy ) const;

        /* Extra utils */
        Vector_3D Number( ) const;
        //Gives number concentration field in part / cm3
        Vector_2D TotalNumber( ) const;
        double TotalNumber_sum( const Vector_2D& cellAreas ) const;
        Vector_1D Overall_Size_Dist( const Vector_2D& cellAreas ) const;
        //Gives 3D volume field in m3 / cm3
        Vector_3D Volume( ) const;
        Vector_2D TotalVolume( ) const;
        Vector_2D TotalArea( ) const;
        double TotalIceMass_sum( const Vector_2D& cellAreas ) const;
        Vector_2D IWC( ) const;
        Vector_2D Extinction( ) const;
        std::tuple<double, int, int> extinctionWidthIndices(const Vector_1D& xCoord, double thres = 0.1) const;
        double extinctionWidth(const Vector_1D& xCoord, double thres = 0.1) const;
        std::tuple<double, int, int> extinctionDepthIndices(const Vector_1D& yCoord, double thres = 0.1) const;
        double extinctionDepth(const Vector_1D& yCoord, double thres = 0.1) const;
        double intYOD(const Vector_1D& dx, const Vector_1D& dy) const;
        Vector_1D PDF_Total( const Vector_2D &cellAreas ) const;
        Vector_1D PDF_Total( const Mesh &m ) const;
        Vector_1D xOD( const Vector_1D& dx ) const;
        Vector_1D yOD( const Vector_1D& dy ) const;
        Vector_2D Radius( ) const;
        double Radius( UInt iNx, UInt jNy ) const;
        Vector_2D EffRadius( ) const;
        double EffRadius( UInt iNx, UInt jNy ) const;
        Vector_2D StdDev( ) const;
        double StdDev( UInt iNx, UInt jNy ) const;

        //This template just lets us minimize copies/moves without writing a bunch of different overloads
        template< typename Vector3D_t, 
                  typename _  = std::enable_if_t<std::is_same_v<std::decay_t<Vector3D_t>, Vector_3D>> 
                >
        void updatePdf( Vector3D_t&& pdf_new ) {
            pdf = std::forward<Vector3D_t>(pdf_new);
        }
        /* utils */
        Vector_1D Average( const Vector_2D &weights,   \
                           const double &totWeight ) const;
        void addPDF( const Aerosol &PDF, const Vector_2D &weights,  \
                     const Vector_2D &cellAreas, const double nCell );
        void addPDF( const Vector_1D &PDF, const Vector_2D &weights, \
                     const Vector_2D &cellAreas, const double nCell );

        /* gets */
        inline const Vector_1D& getBinCenters() const { return bin_Centers; };
        inline const Vector_3D& getBinVCenters() const { return bin_VCenters; };
        inline Vector_3D& getBinVCenters_nonConstRef() { return bin_VCenters; };
        inline const Vector_1D& getBinEdges() const { return bin_Edges; };
        inline const Vector_1D& getBinSizes() const { return bin_Sizes; };
        inline UInt getNBin() const { return nBin; };
        inline const char* getType() const { return type; };
        inline double getAlpha() const { return alpha; };
        inline const Vector_3D& getPDF() const { return pdf; };
        inline Vector_3D& getPDF_nonConstRef() { return pdf; };
        inline int getNx() const { return Nx; }
        inline int getNy() const { return Ny; }


    protected:

        unsigned int Nx, Ny;
        Vector_3D pdf; //Everything with the pdf is implicitly in [ / cm3]
        Vector_3D bin_VCenters;
        Vector_1D bin_Centers;
        Vector_1D bin_Edges;
        Vector_1D bin_VEdges;
        Vector_1D bin_Sizes;
        UInt nBin;
        const char* type;
        double mu;
        double sigma;
        double alpha;

    private:

};

#endif /* AEROSOL_H_INCLUDED */
