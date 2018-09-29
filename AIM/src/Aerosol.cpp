/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*                        AIrcraft Microphysics                     */
/*                              (AIM)                               */
/*                                                                  */
/* Aerosol Program File                                             */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 9/28/2018                                 */
/* File                 : Aerosol.cpp                               */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "Aerosol.hpp"

namespace AIM
{
    
    Aerosol::Aerosol( Vector_1D bin_Centers_, Vector_1D bin_Edges_, RealDouble nPart_, RealDouble mu_, RealDouble sigma_, const char* distType, RealDouble alpha_, RealDouble gamma_, RealDouble b_ ): 
        type( distType ),
        bin_Sizes( bin_Centers_.size() ),
        pdf( bin_Centers_.size() )
    {

        /* Constructor */

        /* In the following cases, the aerosol pdf represents the following quantity: dn/d(ln(r))
         * The following identity can be used:
         * dn/dr = 1/r * dn/d(ln(r)) */
        
        /* Allocate bins centers, edges and number of bins */
        bin_Centers = bin_Centers_;
        bin_Edges = bin_Edges_;
        nBin = bin_Centers.size();

        if ( nBin <= 0 ) {
            std::cout << "\nIn Aerosol::Aerosol: distribution has " << nBin << " bins!";
        }

        if ( nBin + 1 != bin_Edges.size() ) {
            std::cout << "\nIn Aerosol::Aerosol: bin centers and/or edges are misshaped: " << nBin << " != " << bin_Edges.size();
        }

        /* Compute size of each bin */
        for ( unsigned int iBin = 0; iBin < nBin; iBin++ ) {
            bin_Sizes[iBin] = bin_Edges[iBin+1] - bin_Edges[iBin];
        }

        /* Allocate number of particles */
        nPart = nPart_;

        /* Allocate mean and standard deviation */
        if ( mu_ <= 0 ) {
            std::cout << "\nIn Aerosol::Aerosol: mean is negative: mean = " << mu_;
        }
        mu    = mu_;
        if ( sigma_ <= 0 ) {
            std::cout << "\nIn Aerosol::Aerosol: standard deviation is negative: stddev = " << sigma_;
        }
        sigma = sigma_;

        /* Allocate alpha (only used in power law and generalized gamma distributions) */
        alpha = alpha_;

        /* Initialize distribution */
        if ( ( strcmp( type, "log" ) == 0 ) || ( strcmp( type, "lognormal" ) == 0 ) ) {
            /* Log-normal distribution:
             * dn/d(ln(r)) = N / ( sqrt(2*\pi) * ln(sigma) ) * exp( - ( ln(r) - ln(r_m) ) ^ 2 / ( 2 ln(sigma) ^2 )) */

            if ( sigma <= 1.0 ) {
                std::cout << "\nIn Aerosol::Aerosol: log-normal distribution requires that stddev > 1.0 ( stddev = " << sigma << " )";
            }

            for ( unsigned int iBin = 0; iBin < bin_Centers.size(); iBin++ ) {
                pdf[iBin] = nPart * exp( - 0.5 * ( ( log( bin_Centers[iBin] ) - log( mu ) ) / log( sigma ) ) * ( ( log( bin_Centers[iBin] ) - log( mu ) ) / log( sigma ) ) ) / ( sqrt( 2.0 * physConst::PI ) * log( sigma ) );
            }
            
        } else if ( ( strcmp( type, "norm" ) == 0 ) || ( strcmp( type, "normal" ) == 0 ) ) {
            /* Normal distribution:
             * dn/d(ln(r)) = N / ( sqrt(2*\pi) * sigma ) * exp( - ( r - r_m ) ^ 2 / ( 2 * sigma ^2 )) */

            for ( unsigned int iBin = 0; iBin < bin_Centers.size(); iBin++ ) {
                pdf[iBin] = nPart * exp( - 0.5 * ( ( bin_Centers[iBin] - mu ) / sigma ) * ( ( bin_Centers[iBin] - mu ) / sigma ) ) / ( sqrt( 2.0 * physConst::PI ) * sigma );
            }

        } else if ( ( strcmp( type, "pow" ) == 0 ) || ( strcmp( type, "power" ) == 0 ) || ( strcmp( type, "junge" ) == 0 ) ) {
            /* Power distribution:
             * dn/d(ln(r)) = N * alpha * ( r / r_min ) ^ (-alpha) */

            if ( alpha <= 0.0 ) {
                std::cout << "\nIn Aerosol::Aerosol: power law requires that alpha > 0 ( alpha = " << alpha_ << " )";
            }
            
            for ( unsigned int iBin = 0; iBin < bin_Centers.size(); iBin++ ) {
                pdf[iBin] = nPart * alpha * pow( bin_Centers[iBin] / bin_Centers[0], -alpha );
            }

        } else if ( ( strcmp( type, "gam" ) == 0 ) || ( strcmp( type, "gamma" ) == 0 ) || ( strcmp( type, "generalized gamma" ) == 0 ) ) {
            /* Gamma distribution:
             * dn/d(ln(r)) = gamma * b ^ ((alpha + 1)/gamma) / Gamma((alpha+1)/gamma) * r ^ (alpha + 1) * exp( - b * r ^ (gamma) ) 
             * alpha and gamma are of the order of unity (1 - 10)
             * b must be of the order of r ^ (-gamma) >> 1 */

            if ( ( alpha <= 0 ) || (std::fmod( alpha, 1.0 ) != 0 ) ) {
                std::cout << "\nIn Aerosol::Aerosol: (generalized) gamma distribution requires that alpha is a positive integer ( alpha = " << alpha << " )";
            }
            if ( gamma_ <= 0 ) {
                std::cout << "\nIn Aerosol::Aerosol: (generalized) gamma distribution requires that gamma is positive ( gamma = " << gamma_ << " )";
            }
            if ( b_ <= 0 ) {
                std::cout << "\nIn Aerosol::Aerosol: (generalized) gamma distribution requires that b is positive ( b = " << b_ << " )";
            }
            
            for ( unsigned int iBin = 0; iBin < bin_Centers.size(); iBin++ ) {
                pdf[iBin] = nPart * gamma_ * pow( b_, (alpha+1) / gamma_ ) / std::tgamma( (alpha+1) / gamma_ ) * pow( bin_Centers[iBin], alpha + 1 ) * exp( - b_ * pow( bin_Centers[iBin], gamma_ ) );
            }

        } else {
            std::cout << "\nIn Aerosol::Aerosol: distribution type must be either lognormal, normal, power or (generalized) gamma";
            std::cout << "\nCurrent type is " << type;
        }
      
        /* Check that we get the right number of particles */
        if ( std::abs( Moment() - nPart ) / nPart > 0.10 ) {
            std::cout << "\nIn Aerosol::Aerosol: the size range doesn't cover the full distribution";
            std::cout << "\nPrescribed number: " << nPart << ", Number covered: " << Moment() << " [#/cm^3]";
        }
        
        /* For lognormal distributions:
         * rMode = rMean * exp( - ln(sigma)^2 );
         * rMedi = rMean * exp( 0.5 * ln(sigma)^2 ); */

    } /* End of Aerosol::Aerosol */

    Aerosol::Aerosol( const Aerosol &rhs )
    {

        bin_Centers = rhs.bin_Centers;
        bin_Edges = rhs.bin_Edges;
        bin_Sizes = rhs.bin_Sizes;
        nBin = rhs.nBin;
        nPart = rhs.nPart;
        type = rhs.type;
        mu = rhs.mu;
        sigma = rhs.sigma;
        alpha = rhs.alpha;
        pdf = rhs.pdf;

    } /* End of Aerosol::Aerosol */
    
    Aerosol::~Aerosol( )
    {

        /* Destructor */

    } /* End of Aerosol::~Aerosol */

    Aerosol Aerosol::operator=( const Aerosol &rhs )
    {

        if ( &rhs == this )
            return *this;

        bin_Centers = rhs.bin_Centers;
        bin_Edges = rhs.bin_Edges;
        bin_Sizes = rhs.bin_Sizes;
        nBin = rhs.nBin;
        nPart = rhs.nPart;
        type = rhs.type;
        mu = rhs.mu;
        sigma = rhs.sigma;
        alpha = rhs.alpha;
        pdf = rhs.pdf;
        return *this;

    } /* End of Aerosol::operator= */

    Aerosol Aerosol::operator+=( const Aerosol &rhs )
    {

        if ( nBin != rhs.getNBin() ) {
            std::cout << "\nIn Aerosol::operator+=: aerosol distributions do not have the same number of bins: " << nBin << " != " << rhs.getNBin();
            return *this;
        }

        Vector_1D bin_Centers_rhs = rhs.getBinCenters();
        for ( unsigned int iBin = 0; iBin < nBin; iBin++ ) {
            if ( std::abs( bin_Centers[iBin] - bin_Centers_rhs[iBin] ) > 1.0E-10 ) {
                std::cout << "\nIn Aerosol::operator+=: aerosol distributions do not have the same bin centers!";
                return *this;
            }
        }

        Vector_1D pdf_rhs = rhs.getPDF();
        for ( unsigned int iBin = 0; iBin < nBin; iBin++ ) {
            pdf[iBin] += pdf_rhs[iBin];
        }

        return *this;

    } /* End of Aerosol::operator+= */
    
    Aerosol Aerosol::operator-=( const Aerosol &rhs )
    {

        if ( nBin != rhs.getNBin() ) {
            std::cout << "\nIn Aerosol::operator+=: aerosol distributions do not have the same number of bins: " << nBin << " != " << rhs.getNBin();
            return *this;
        }

        Vector_1D bin_Centers_rhs = rhs.getBinCenters();
        for ( unsigned int iBin = 0; iBin < nBin; iBin++ ) {
            if ( std::abs( bin_Centers[iBin] - bin_Centers_rhs[iBin] ) > 1.0E-10 ) {
                std::cout << "\nIn Aerosol::operator+=: aerosol distributions do not have the same bin centers!";
                return *this;
            }
        }

        Vector_1D pdf_rhs = rhs.getPDF();
        for ( unsigned int iBin = 0; iBin < nBin; iBin++ ) {
            pdf[iBin] -= pdf_rhs[iBin];
        }

        return *this;

    } /* End of Aerosol::operator-= */
    
    Aerosol Aerosol::operator+( const Aerosol &rhs ) const
    {

        Aerosol result = *this;
        
        if ( nBin != rhs.getNBin() ) {
            std::cout << "\nIn Aerosol::operator+=: aerosol distributions do not have the same number of bins: " << nBin << " != " << rhs.getNBin();
            return *this;
        }

        Vector_1D bin_Centers_rhs = rhs.getBinCenters();
        for ( unsigned int iBin = 0; iBin < nBin; iBin++ ) {
            if ( std::abs( bin_Centers[iBin] - bin_Centers_rhs[iBin] ) > 1.0E-10 ) {
                std::cout << "\nIn Aerosol::operator+=: aerosol distributions do not have the same bin centers!";
                return *this;
            }
        }

        Vector_1D pdf_rhs = rhs.getPDF();
        for ( unsigned int iBin = 0; iBin < nBin; iBin++ ) {
            result.pdf[iBin] += pdf_rhs[iBin];
        }
        
        return result;

    } /* End of Aerosol::operator+ */
    
    Aerosol Aerosol::operator-( const Aerosol &rhs ) const
    {

        Aerosol result = *this;
        
        if ( nBin != rhs.getNBin() ) {
            std::cout << "\nIn Aerosol::operator+=: aerosol distributions do not have the same number of bins: " << nBin << " != " << rhs.getNBin();
            return *this;
        }

        Vector_1D bin_Centers_rhs = rhs.getBinCenters();
        for ( unsigned int iBin = 0; iBin < nBin; iBin++ ) {
            if ( std::abs( bin_Centers[iBin] - bin_Centers_rhs[iBin] ) > 1.0E-10 ) {
                std::cout << "\nIn Aerosol::operator+=: aerosol distributions do not have the same bin centers!";
                return *this;
            }
        }

        Vector_1D pdf_rhs = rhs.getPDF();
        for ( unsigned int iBin = 0; iBin < nBin; iBin++ ) {
            result.pdf[iBin] -= pdf_rhs[iBin];
        }
        
        return result;

    } /* End of Aerosol::operator- */
        
    RealDouble Aerosol::Moment( unsigned int n ) const
    {

        double moment = 0;
        
        for ( unsigned int iBin = 0; iBin < nBin; iBin++ ) {
            moment += ( log( bin_Edges[iBin+1] ) - log( bin_Edges[iBin] ) ) * pow( bin_Centers[iBin], n ) * pdf[iBin];
        }

        return moment;

    } /* End of Aerosol::Moment */

    RealDouble Aerosol::getRadius( ) const
    {

        double N = Moment( 0 );

        if ( N > 0 ) {
            return Moment( 1 ) / N;
        } else {
            std::cout << "\nIn Aerosol::getRadius: Number of particles is " << N << " <= 0";
            return 0;
        }

    } /* End of Aerosol::getRadius */
    
    RealDouble Aerosol::getEffRadius( ) const
    {

        double m2 = Moment(2);

        if ( m2 > 0 ) {
            return Moment(3) / m2;
        } else {
            std::cout << "\nIn Aerosol::getEffRadius: Second moment is " << m2 << " <= 0";
            return 0;
        }

    } /* End of Aerosol::getEffRadius */

    RealDouble Aerosol::getStdDev( ) const
    {

        double N = Moment( 0 );

        if ( N > 0 ) {
            return sqrt( Moment( 2 ) / N - pow( Moment( 1 ) / N, 2.0 ) );
        } else {
            std::cout << "\nIn Aerosol::getRadius: Number of particles is " << N << " <= 0";
            return 0;
        }

    } /* End of Aerosol::getStdDev */

    Vector_1D Aerosol::getBinCenters() const
    {

        return bin_Centers;

    } /* End of Aerosol::getBinCenters */
        
    Vector_1D Aerosol::getBinEdges() const
    {
    
        return bin_Edges;

    } /* End of Aerosol::getBinEdges */
    
    Vector_1D Aerosol::getBinSizes() const
    {
    
        return bin_Sizes;

    } /* End of Aerosol::getBinSizes */
    
    unsigned int Aerosol::getNBin() const
    {

        return nBin;

    } /* End of Aerosol::getNBin */

    RealDouble Aerosol::getNPart() const
    {

        return nPart;

    } /* End of Aerosol::getNPart */
    
    const char* Aerosol::getType() const
    {
    
        return type;

    } /* End of Aerosol::getType */

    RealDouble Aerosol::getAlpha() const
    {

        return alpha;

    } /* End of Aerosol::getAlpha */

    Vector_1D Aerosol::getPDF() const
    {

        return pdf;

    } /* End of Aerosol::getPDF */

}

/* End of Aerosol.cpp */
