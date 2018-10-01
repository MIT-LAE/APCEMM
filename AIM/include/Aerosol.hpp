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

#ifndef DISTRIBUTION_H_INCLUDED
#define DISTRIBUTION_H_INCLUDED

#include <iostream>
#include <cmath>
#include <vector>
#include <cstring>

#include "../../Headers/ForwardsDecl.hpp"
#include "../../Headers/PhysConstant.hpp"

namespace AIM
{

    class Aerosol;

}

class AIM::Aerosol
{

    public:

        Aerosol( Vector_1D bin_Centers, Vector_1D bin_Edges, RealDouble nPart, RealDouble mu, RealDouble sigma, const char* distType = "lognormal", RealDouble alpha_ = -1.0, RealDouble gamma_ = -1.0, RealDouble b_ = 0.0 );
        Aerosol( const Aerosol &rhs );
        ~Aerosol();
        Aerosol operator=( const Aerosol &rhs );
        Aerosol operator+=( const Aerosol &rhs );
        Aerosol operator-=( const Aerosol &rhs );
        Aerosol operator+( const Aerosol &rhs ) const;
        Aerosol operator-( const Aerosol &rhs ) const;
        RealDouble Moment( unsigned int n = 0 ) const;
        RealDouble getRadius( ) const;
        RealDouble getEffRadius( ) const;
        RealDouble getStdDev( ) const;

        Vector_1D getBinCenters() const;
        Vector_1D getBinEdges() const;
        Vector_1D getBinSizes() const;
        unsigned int getNBin() const;
        RealDouble getNPart() const;
        const char* getType() const;
        RealDouble getAlpha() const;
        Vector_1D getPDF() const;

    protected:

        Vector_1D bin_Centers;
        Vector_1D bin_Edges;
        Vector_1D bin_Sizes;
        unsigned int nBin;
        RealDouble nPart;
        const char* type;
        RealDouble mu;
        RealDouble sigma;
        RealDouble alpha;
        Vector_1D pdf;

    private:

};

#endif /* DISTRIBUTION_H_INCLUDED */
