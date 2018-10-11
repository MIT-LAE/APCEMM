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

#include "Util/ForwardsDecl.hpp"
#include "Util/PhysConstant.hpp"
#include "AIM/Coagulation.hpp"

namespace AIM
{

    class Aerosol;

}

class AIM::Aerosol
{

    public:

        /* Constructor */
        Aerosol( Vector_1D bin_Centers, Vector_1D bin_Edges, RealDouble nPart, RealDouble mu, RealDouble sigma, const char* distType = "lognormal", RealDouble alpha_ = -1.0, RealDouble gamma_ = -1.0, RealDouble b_ = 0.0 );

        /* Copy */
        Aerosol( const Aerosol &rhs );

        /* Destructor */
        ~Aerosol();

        /* Operators */
        Aerosol operator=( const Aerosol &rhs );
        Aerosol operator+=( const Aerosol &rhs );
        Aerosol operator-=( const Aerosol &rhs );
        Aerosol operator+( const Aerosol &rhs ) const;
        Aerosol operator-( const Aerosol &rhs ) const;

        /* Coagulation */
        void Coagulate( const RealDouble dt, const Coagulation &kernel );
        void Coagulate( const RealDouble dt, const Vector_2D &beta, const Vector_3D &f );
        
        /* Moments */
        RealDouble Moment( UInt n = 0 ) const;
        RealDouble getRadius( ) const;
        RealDouble getEffRadius( ) const;
        RealDouble getStdDev( ) const;

        /* utils */
        void scalePdf( RealDouble factor );

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

#endif /* AEROSOL_H_INCLUDED */
