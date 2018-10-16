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

#include "AIM/Aerosol.hpp"

namespace AIM
{

    Aerosol::Aerosol( ):
        bin_Centers( 24 ),
        bin_VCenters( 24 ),
        bin_Edges( 25 ),
        bin_Sizes( 24 ),
        nBin( 24 ),
        pdf( 24 )
    {

        /* Default constructor */

        for ( UInt iBin = 0; iBin < nBin + 1; iBin++ ) {
            bin_Edges[iBin] = 0.1E-07 * pow( 1.2, iBin / RealDouble(3.0) );
        }

        for ( UInt iBin = 0; iBin < nBin; iBin++ ) {
            bin_Centers[iBin] = 0.5 * ( bin_Edges[iBin] + bin_Edges[iBin+1] );
            bin_VCenters[iBin] = 4.0 / RealDouble(3.0) * physConst::PI * bin_Centers[iBin] * bin_Centers[iBin] * bin_Centers[iBin];
            bin_Sizes[iBin] = bin_Edges[iBin+1] - bin_Edges[iBin];
        }

        nPart = 1.0;

        type = "lognormal";

        mu = 0;
        sigma = 0;
        alpha = 0;

        for ( UInt iBin = 0; iBin < nBin; iBin++ )
            pdf[iBin] = 0.0;
        

    } /* End of Aerosol::Aerosol */

    Aerosol::Aerosol( Vector_1D bin_Centers_, Vector_1D bin_Edges_, RealDouble nPart_, RealDouble mu_, RealDouble sigma_, const char* distType, RealDouble alpha_, RealDouble gamma_, RealDouble b_ ): 
        bin_VCenters( bin_Centers_.size() ),
        bin_Sizes( bin_Centers_.size() ),
        type( distType ),
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
            std::cout << "\nIn Aerosol::Aerosol: distribution has " << nBin << " bins!\n";
        }

        if ( nBin + 1 != bin_Edges.size() ) {
            std::cout << "\nIn Aerosol::Aerosol: bin centers and/or edges are misshaped: " << nBin + 1 << " != " << bin_Edges.size() << "\n";
        }

        /* Compute size of each bin */
        for ( UInt iBin = 0; iBin < nBin; iBin++ ) {
            bin_Sizes[iBin] = bin_Edges[iBin+1] - bin_Edges[iBin];
            bin_VCenters[iBin] = 4.0 / RealDouble(3.0) * physConst::PI * bin_Centers_[iBin] * bin_Centers_[iBin] * bin_Centers_[iBin];
        }

        /* Allocate number of particles */
        nPart = nPart_;

        /* Allocate mean and standard deviation */
        if ( mu_ <= 0 ) {
            std::cout << "\nIn Aerosol::Aerosol: mean/mode is negative: mu = " << mu_ << "\n";
        }
        mu    = mu_;
        if ( sigma_ <= 0 ) {
            std::cout << "\nIn Aerosol::Aerosol: standard deviation is negative: stddev = " << sigma_ << "\n";
        }
        sigma = sigma_;

        /* Allocate alpha (only used in power law and generalized gamma distributions) */
        alpha = alpha_;

        /* Initialize distribution */
        if ( ( strcmp( type, "log" ) == 0 ) || ( strcmp( type, "lognormal" ) == 0 ) ) {
            /* Log-normal distribution:
             * dn/d(ln(r)) = N / ( sqrt(2*\pi) * ln(sigma) ) * exp( - ( ln(r) - ln(r_m) ) ^ 2 / ( 2 ln(sigma) ^2 )) */

            if ( sigma <= 1.0 ) {
                std::cout << "\nIn Aerosol::Aerosol: log-normal distribution requires that stddev > 1.0 ( stddev = " << sigma << " )\n";
            }

            for ( UInt iBin = 0; iBin < bin_Centers.size(); iBin++ ) {
                pdf[iBin] = nPart * exp( - 0.5 * ( ( log( bin_Centers[iBin] ) - log( mu ) ) / log( sigma ) ) * ( ( log( bin_Centers[iBin] ) - log( mu ) ) / log( sigma ) ) ) / ( sqrt( 2.0 * physConst::PI ) * log( sigma ) );
            }
            
        } else if ( ( strcmp( type, "norm" ) == 0 ) || ( strcmp( type, "normal" ) == 0 ) ) {
            /* Normal distribution:
             * dn/d(ln(r)) = N / ( sqrt(2*\pi) * sigma ) * exp( - ( r - r_m ) ^ 2 / ( 2 * sigma ^2 )) */

            for ( UInt iBin = 0; iBin < bin_Centers.size(); iBin++ ) {
                pdf[iBin] = nPart * exp( - 0.5 * ( ( bin_Centers[iBin] - mu ) / sigma ) * ( ( bin_Centers[iBin] - mu ) / sigma ) ) / ( sqrt( 2.0 * physConst::PI ) * sigma );
            }

        } else if ( ( strcmp( type, "pow" ) == 0 ) || ( strcmp( type, "power" ) == 0 ) || ( strcmp( type, "junge" ) == 0 ) ) {
            /* Power distribution:
             * dn/d(ln(r)) = N * alpha * ( r / r_min ) ^ (-alpha) */

            if ( alpha <= 0.0 ) {
                std::cout << "\nIn Aerosol::Aerosol: power law requires that alpha > 0 ( alpha = " << alpha_ << " )\n";
            }
            
            for ( UInt iBin = 0; iBin < bin_Centers.size(); iBin++ ) {
                pdf[iBin] = nPart * alpha * pow( bin_Centers[iBin] / bin_Centers[0], -alpha );
            }

        } else if ( ( strcmp( type, "gam" ) == 0 ) || ( strcmp( type, "gamma" ) == 0 ) || ( strcmp( type, "generalized gamma" ) == 0 ) ) {
            /* Gamma distribution:
             * dn/d(ln(r)) = gamma * b ^ ((alpha + 1)/gamma) / Gamma((alpha+1)/gamma) * r ^ (alpha + 1) * exp( - b * r ^ (gamma) ) 
             * alpha and gamma are of the order of unity (1 - 10)
             * b must be of the order of r ^ (-gamma) >> 1 */

            if ( ( alpha <= 0 ) || (std::fmod( alpha, 1.0 ) != 0 ) ) {
                std::cout << "\nIn Aerosol::Aerosol: (generalized) gamma distribution requires that alpha is a positive integer ( alpha = " << alpha << " )\n";
            }
            if ( gamma_ <= 0 ) {
                std::cout << "\nIn Aerosol::Aerosol: (generalized) gamma distribution requires that gamma is positive ( gamma = " << gamma_ << " )\n";
            }
            if ( b_ <= 0 ) {
                std::cout << "\nIn Aerosol::Aerosol: (generalized) gamma distribution requires that b is positive ( b = " << b_ << " )\n";
            }
            
            for ( UInt iBin = 0; iBin < bin_Centers.size(); iBin++ ) {
                pdf[iBin] = nPart * gamma_ * pow( b_, (alpha+1) / gamma_ ) / boost::math::tgamma( (alpha+1) / gamma_ ) * pow( bin_Centers[iBin], alpha + 1 ) * exp( - b_ * pow( bin_Centers[iBin], gamma_ ) );
            }

        } else {
            std::cout << "\nIn Aerosol::Aerosol: distribution type must be either lognormal, normal, power or (generalized) gamma\n";
            std::cout << "\nCurrent type is " << type << "\n";
        }
      
        /* Check that we get the right number of particles */
        if ( ( std::abs( Moment() - nPart ) / nPart > 0.10 ) && ( nPart > 1.0E-10 )) {
            std::cout << "\nIn Aerosol::Aerosol: the size range doesn't cover the full distribution";
            std::cout << "\nPrescribed number: " << nPart << ", Number covered: " << Moment() << " [#/cm^3]\n";
            if ( ( strcmp( type, "log" ) == 0 ) || ( strcmp( type, "lognormal" ) == 0 ) ) {
                std::cout << "For lognormal distribution, prescribed mode is: " << mu << " [m]\n";
            } else if ( ( strcmp( type, "norm" ) == 0 ) || ( strcmp( type, "normal" ) == 0 ) ) {
                std::cout << "For normal distribution, prescribed mode is: " << mu << " [m]\n";
            }
        }
        
        /* For lognormal distributions:
         * rMode = rMean * exp( - ln(sigma)^2 );
         * rMedi = rMean * exp( 0.5 * ln(sigma)^2 ); */

    } /* End of Aerosol::Aerosol */

    Aerosol::Aerosol( const Aerosol &rhs )
    {

        bin_Centers = rhs.bin_Centers;
        bin_VCenters = rhs.bin_VCenters;
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

    Aerosol& Aerosol::operator=( const Aerosol &rhs )
    {

        if ( &rhs == this )
            return *this;

        bin_Centers = rhs.bin_Centers;
        bin_VCenters = rhs.bin_VCenters;
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
        for ( UInt iBin = 0; iBin < nBin; iBin++ ) {
            if ( std::abs( bin_Centers[iBin] - bin_Centers_rhs[iBin] ) > 1.0E-10 ) {
                std::cout << "\nIn Aerosol::operator+=: aerosol distributions do not have the same bin centers!";
                return *this;
            }
        }

        Vector_1D pdf_rhs = rhs.getPDF();
        for ( UInt iBin = 0; iBin < nBin; iBin++ ) {
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
        for ( UInt iBin = 0; iBin < nBin; iBin++ ) {
            if ( std::abs( bin_Centers[iBin] - bin_Centers_rhs[iBin] ) > 1.0E-10 ) {
                std::cout << "\nIn Aerosol::operator+=: aerosol distributions do not have the same bin centers!";
                return *this;
            }
        }

        Vector_1D pdf_rhs = rhs.getPDF();
        for ( UInt iBin = 0; iBin < nBin; iBin++ ) {
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
        for ( UInt iBin = 0; iBin < nBin; iBin++ ) {
            if ( std::abs( bin_Centers[iBin] - bin_Centers_rhs[iBin] ) > 1.0E-10 ) {
                std::cout << "\nIn Aerosol::operator+=: aerosol distributions do not have the same bin centers!";
                return *this;
            }
        }

        Vector_1D pdf_rhs = rhs.getPDF();
        for ( UInt iBin = 0; iBin < nBin; iBin++ ) {
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
        for ( UInt iBin = 0; iBin < nBin; iBin++ ) {
            if ( std::abs( bin_Centers[iBin] - bin_Centers_rhs[iBin] ) > 1.0E-10 ) {
                std::cout << "\nIn Aerosol::operator+=: aerosol distributions do not have the same bin centers!";
                return *this;
            }
        }

        Vector_1D pdf_rhs = rhs.getPDF();
        for ( UInt iBin = 0; iBin < nBin; iBin++ ) {
            result.pdf[iBin] -= pdf_rhs[iBin];
        }
        
        return result;

    } /* End of Aerosol::operator- */
       
    void Aerosol::Coagulate( const RealDouble dt, const Coagulation &kernel )
    {

        /* DESCRIPTION:
         * Performs self-coagulation. Updates Aerosol.pdf
         * The numerical scheme is taken from: 
         * M.Z. Jacobson, Fundamentals of atmospheric modeling. Cambridge university press, 2005.*/
         
        /* INPUT:
         * - RealDouble dt      :: Timestep in s
         * - Coagulation kernel :: Coagulation structure containing the coagulation kernels
         *
         * OUTPUT:
         *
         */

        /* Allocate variables */
        Vector_1D P( nBin, 0.0 );
        Vector_1D L( nBin, 0.0 );
        Vector_1D v( nBin );
        UInt iBin, jBin, kBin;

        for ( iBin = 0; iBin < nBin; iBin++ ) {
            v[iBin] = pdf[iBin] * bin_VCenters[iBin] * 1.0E+06;
            /* Unit check:
             * [#/cm^3] * [m^3] * [cm^3/m^3] = [cm^3/cm^3]*/
        }
        
        /* Copy v into v_new */
        Vector_1D v_new = v;

        /* \frac{dv}{dt}[iBin] = P - L * v[iBin] 
         * Production     P = sum of all the bins (smaller than iBin) that coagulate into a particle of size iBin
         * Loss L * v[iBin] = sum of all the bins that coagulate with iBin 
         *
         * Scheme 1:
         * v_new - v = ( P - L * v ) * dt 
         * Scheme 2:
         * v_new - v = ( P - L * v_new ) * dt
         * v_new = ( v + P * dt ) / ( 1.0 + L ) 
         * The latter is mass-conserving */

        /* Can this be improved? 
         * Can for loops be removed? */
        for ( iBin = 0; iBin < nBin; iBin++ ) {

            /* Build production and loss terms */
            for ( jBin = 0; jBin < nBin; jBin++ ) {
                if ( jBin <= iBin ) {
            
                    for ( kBin = 0; kBin < iBin ; kBin++ ) {
                        /* k coagulating with j to form i */
                        if ( kernel.f[kBin][jBin][iBin] != 0.0 )
                            P[iBin] += kernel.f[kBin][jBin][iBin] * kernel.beta[kBin][jBin] * v_new[kBin] * pdf[jBin];
                    }

                }
                /* i coagulating with j to deplete i */
                if ( kernel.f[iBin][jBin][iBin] != 1.0 ) 
                    L[iBin] += ( 1.0 - kernel.f[iBin][jBin][iBin] ) * kernel.beta[iBin][jBin] * pdf[jBin];
            }

            /* Non-mass conserving scheme: */
//            v_new[iBin] = v[iBin] + ( P[iBin] - L[iBin] * v[iBin] ) * dt;

            /* Mass conserving scheme: */
            v_new[iBin] = ( v[iBin] + dt * P[iBin] ) / ( 1.0 + dt * L[iBin] );
        }

        /* For debug purposes */
        bool checkMass = 0;
        if ( checkMass ) {
            std::cout << "Coagulation is a volume-conserving process. The two following quantities should be identical\n";
            std::cout << "At t     : " << Moment(3) << "\n";
        }

        for ( iBin = 0; iBin < nBin; iBin++ ) {
            pdf[iBin] *= v_new[iBin] / v[iBin]; 
        }
        
        if ( checkMass )
            std::cout << "At t + dt: " << Moment(3) << "\n";


    } /* End of Aerosol::Coagulate */
    
    void Aerosol::Coagulate( const RealDouble dt, const Vector_2D &beta, const Vector_3D &f )
    {

        /* DESCRIPTION:
         * Performs self-coagulation. Updates Aerosol.pdf
         * The numerical scheme is taken from: 
         * M.Z. Jacobson, Fundamentals of atmospheric modeling. Cambridge university press, 2005.*/
         
        /* INPUT:
         * - RealDouble dt   :: Timestep in s
         * - Vector_2D  beta :: Coagulation kernel
         * - Vector_3D  f    :: Coagulation vector: "In which bin does the result of the coagulation of two other bins end up?"
         *
         * OUTPUT:
         *
         */

        /* Allocate variables */
        Vector_1D P( nBin );
        Vector_1D L( nBin );
        Vector_1D v( nBin );
        UInt iBin, jBin, kBin;

        for ( iBin = 0; iBin < nBin; iBin++ ) {
            v[iBin] = pdf[iBin] * bin_VCenters[iBin] * 1.0E+06;
            /* Unit check:
             * [#/cm^3] * [m^3] * [cm^3/m^3] = [cm^3/cm^3]*/
        }
        
        /* Copy v into v_new */
        Vector_1D v_new = v;

        /* \frac{dv}{dt}[iBin] = P - L * v[iBin] 
         * Production     P = sum of all the bins (smaller than iBin) that coagulate into a particle of size iBin
         * Loss L * v[iBin] = sum of all the bins that coagulate with iBin 
         *
         * Scheme 1:
         * v_new - v = ( P - L * v ) * dt 
         * Scheme 2:
         * v_new - v = ( P - L * v_new ) * dt
         * v_new = ( v + P * dt ) / ( 1.0 + L ) 
         * The latter is mass-conserving */

        /* Can this be improved? 
         * Can for loops be removed? */
        for ( iBin = 0; iBin < nBin; iBin++ ) {
            for ( jBin = 0; jBin < nBin; jBin++ ) {
                if ( jBin <= iBin ) {
                    for ( kBin = 0; kBin < iBin; kBin++ ) {
                        P[iBin] += f[kBin][jBin][iBin] * beta[kBin][jBin] * v_new[kBin] * pdf[jBin];
                    }
                }
                L[iBin] += ( 1.0 - f[iBin][jBin][iBin] ) * beta[iBin][jBin] * pdf[jBin];
            }
            v_new[iBin] = ( v[iBin] + dt * P[iBin] ) / ( 1.0 + dt * L[iBin] );
            pdf[iBin] *= v_new[iBin] / v[iBin]; 
        }

    } /* End of Aerosol::Coagulate */

    RealDouble Aerosol::Moment( UInt n ) const
    {

        RealDouble moment = 0;

        for ( UInt iBin = 0; iBin < nBin; iBin++ ) {
            moment += ( log( bin_Edges[iBin+1] ) - log( bin_Edges[iBin] ) ) * pow( bin_Centers[iBin], n ) * pdf[iBin];
        }

        return moment;

    } /* End of Aerosol::Moment */

    RealDouble Aerosol::getRadius( ) const
    {

        RealDouble N = Moment( 0 );

        if ( N > 0 ) {
            return Moment( 1 ) / N;
        } else {
            std::cout << "\nIn Aerosol::getRadius: Number of particles is " << N << " <= 0\n";
            return 0;
        }

    } /* End of Aerosol::getRadius */
    
    RealDouble Aerosol::getEffRadius( ) const
    {

        RealDouble m2 = Moment( 2 );

        if ( m2 > 0 ) {
            return Moment(3) / m2;
        } else {
            std::cout << "\nIn Aerosol::getEffRadius: Second moment is " << m2 << " <= 0\n";
            return 0;
        }

    } /* End of Aerosol::getEffRadius */

    RealDouble Aerosol::getStdDev( ) const
    {

        RealDouble N = Moment( 0 );

        if ( N > 0 ) {
            return sqrt( Moment( 2 ) / N - pow( Moment( 1 ) / N, 2.0 ) );
        } else {
            std::cout << "\nIn Aerosol::getRadius: Number of particles is " << N << " <= 0\n";
            return 0;
        }

    } /* End of Aerosol::getStdDev */

    void Aerosol::scalePdf( RealDouble scalFactor )
    {

        if ( scalFactor <= 0 ) {
            std::cout << "\nIn Aerosol::scalePdf: scaling factor is negative ( " << scalFactor << " <= 0 )\n";
        }

        for ( UInt iBin = 0; iBin < nBin; iBin++ ) {
            pdf[iBin] *= scalFactor;
        }

    } /* End of Aerosol::scalePdf */

    void Aerosol::updatePdf( Vector_1D pdf_ )
    {

        if ( pdf_.size() != pdf.size() ) {
            std::cout << "\nIn Aerosol::updatePdf:: sizes differ ( " << pdf_.size() << " != " << pdf.size() << " )\n";
        }
        else {
            pdf = pdf_;
        }

    } /* End of Aerosol::updatePdf */

    Vector_1D Aerosol::getBinCenters() const
    {

        return bin_Centers;

    } /* End of Aerosol::getBinCenters */
    
    Vector_1D Aerosol::getBinVCenters() const
    {

        return bin_VCenters;

    } /* End of Aerosol::getBinVCenters */
        
    Vector_1D Aerosol::getBinEdges() const
    {
    
        return bin_Edges;

    } /* End of Aerosol::getBinEdges */
    
    Vector_1D Aerosol::getBinSizes() const
    {
    
        return bin_Sizes;

    } /* End of Aerosol::getBinSizes */
    
    UInt Aerosol::getNBin() const
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
