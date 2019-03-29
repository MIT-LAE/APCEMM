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
        bin_Centers( 2 ),
        bin_VCenters( 2 ),
        bin_Edges( 3 ),
        bin_Sizes( 2 ),
        nBin( 2 ),
        pdf( 2 )
    {

        /* Default constructor */

        for ( UInt iBin = 0; iBin < nBin + 1; iBin++ ) {
            bin_Edges[iBin] = 0.1E-07 * pow( 1.2, iBin / RealDouble(3.0) );
        }

        for ( UInt iBin = 0; iBin < nBin; iBin++ ) {
            bin_Centers[iBin] = 0.5 * ( bin_Edges[iBin] + bin_Edges[iBin+1] );
            bin_VCenters[iBin] = 4.0 / RealDouble(3.0) * physConst::PI * ( bin_Edges[iBin] * bin_Edges[iBin] * bin_Edges[iBin] + bin_Edges[iBin+1] * bin_Edges[iBin+1] * bin_Edges[iBin+1] ) * 0.5;
            bin_Sizes[iBin] = bin_Edges[iBin+1] - bin_Edges[iBin];
        }

        nPart = 1.0;

        type = "lognormal";

        mu = 0;
        sigma = 0;
        alpha = 0;

        for ( UInt iBin = 0; iBin < nBin; iBin++ )
            pdf[iBin] = 0.0E+00;
        

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
            bin_VCenters[iBin] = 4.0 / RealDouble(3.0) * physConst::PI * ( bin_Edges_[iBin] * bin_Edges_[iBin] * bin_Edges_[iBin] + bin_Edges_[iBin+1] * bin_Edges_[iBin+1] * bin_Edges_[iBin+1] ) * 0.5;
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
        Vector_1D P( nBin, 0.0E+00 );
        Vector_1D L( nBin, 0.0E+00 );
        Vector_1D v( nBin );
        UInt iBin, jBin, kBin;

        for ( iBin = 0; iBin < nBin; iBin++ ) {
            v[iBin] = pdf[iBin] * bin_VCenters[iBin] * 1.0E+06;
            /* Unit check:
             * [#/cm^3] * [m^3] * [cm^3/m^3] = [cm^3/cm^3] */
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
        const bool checkMass = 0;
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
             * [#/cm^3] * [m^3] * [cm^3/m^3] = [cm^3/cm^3] */
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

    RealDouble Aerosol::Radius( ) const
    {

        const RealDouble N = Moment( 0 );

        if ( N > 0 ) {
            return Moment( 1 ) / N;
        } else {
            std::cout << "\nIn Aerosol::Radius: Number of particles is " << N << " <= 0\n";
            return 0;
        }

    } /* End of Aerosol::Radius */
    
    RealDouble Aerosol::EffRadius( ) const
    {

        const RealDouble m2 = Moment( 2 );

        if ( m2 > 0 ) {
            return Moment(3) / m2;
        } else {
            std::cout << "\nIn Aerosol::EffRadius: Second moment is " << m2 << " <= 0\n";
            return 0;
        }

    } /* End of Aerosol::EffRadius */

    RealDouble Aerosol::StdDev( ) const
    {

        const RealDouble N = Moment( 0 );

        if ( N > 0 ) {
            return sqrt( Moment( 2 ) / N - pow( Moment( 1 ) / N, 2.0 ) );
        } else {
            std::cout << "\nIn Aerosol::StdDev: Number of particles is " << N << " <= 0\n";
            return 0;
        }

    } /* End of Aerosol::StdDev */

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


    
    
    
    Grid_Aerosol::Grid_Aerosol( ):
        Nx( 2 ),
        Ny( 2 ),
        bin_Centers( 2 ),
        bin_Edges( 3 ),
        bin_VEdges( 3 ),
        bin_Sizes( 2 ),
        nBin( 2 )
    {

        /* Default constructor */

        for ( UInt iBin = 0; iBin < nBin + 1; iBin++ ) {
            bin_Edges[iBin]  = 0.1E-07 * pow( 1.2, iBin / RealDouble(3.0) );
            bin_VEdges[iBin] = 4.0 / (double) 3.0 * physConst::PI * bin_Edges[iBin] * \
                               bin_Edges[iBin] * bin_Edges[iBin];
        }

        for ( UInt iBin = 0; iBin < nBin; iBin++ ) {
            bin_Centers[iBin] = 0.5 * ( bin_Edges[iBin] + bin_Edges[iBin+1] );
            bin_Sizes[iBin] = bin_Edges[iBin+1] - bin_Edges[iBin];
        }
        
        bin_VCenters.resize( nBin, Vector_2D( Ny, Vector_1D( Nx, 0.0 ) ) );
       
        RealDouble vol = 0.0E+00;
        for ( UInt iBin = 0; iBin < nBin; iBin++ ) {
            vol = 0.5 * ( bin_VEdges[iBin] + bin_VEdges[iBin+1] );
            for ( UInt jNy = 0; jNy < Ny; jNy++ ) {
                for ( UInt iNx = 0; iNx < Nx; iNx++ ) {
                    bin_VCenters[iBin][jNy][iNx] = vol;
                }
            }
        }

        nPart = 1.0;

        type = "lognormal";

        mu = 0;
        sigma = 0;
        alpha = 0;
        
        pdf.resize( nBin, Vector_2D( Ny, Vector_1D( Nx, 0.0E+00 ) ) );

        for ( UInt iBin = 0; iBin < nBin; iBin++ ) {
            for ( UInt jNy = 0; jNy < Ny; jNy++ ) {
                for ( UInt iNx = 0; iNx < Nx; iNx++ ) {
                    pdf[iBin][jNy][iNx] = 0.0E+00;
                }
            }
        }
        

    } /* End of Grid_Aerosol::Grid_Aerosol */

    Grid_Aerosol::Grid_Aerosol( UInt Nx_, UInt Ny_, Vector_1D bin_Centers_, Vector_1D bin_Edges_, RealDouble nPart_, RealDouble mu_, RealDouble sigma_, const char* distType, RealDouble alpha_, RealDouble gamma_, RealDouble b_ ): 
        Nx( Nx_ ),
        Ny( Ny_ ),
        bin_VEdges( bin_Centers_.size() + 1 ),
        bin_Sizes( bin_Centers_.size() ),
        type( distType )
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
            std::cout << "\nIn Grid_Aerosol::Grid_Aerosol: distribution has " << nBin << " bins!\n";
        }

        if ( nBin + 1 != bin_Edges.size() ) {
            std::cout << "\nIn Grid_Aerosol::Grid_Aerosol: bin centers and/or edges are misshaped: " << nBin + 1 << " != " << bin_Edges.size() << "\n";
        }

        /* Compute size of each bin */
        for ( UInt iBin = 0; iBin < nBin; iBin++ ) {
            bin_VEdges[iBin] = 4.0 / (RealDouble) 3.0 * physConst::PI * bin_Edges[iBin] * \
                               bin_Edges[iBin] * bin_Edges[iBin];
            bin_Sizes[iBin]  = bin_Edges[iBin+1] - bin_Edges[iBin];
        }
        bin_VEdges[nBin] = 4.0 / (RealDouble) 3.0 * physConst::PI * bin_Edges[nBin] * \
                           bin_Edges[nBin] * bin_Edges[nBin];

        bin_VCenters.resize( nBin, Vector_2D( Ny, Vector_1D( Nx, 0.0E+00 ) ) );

        RealDouble vol = 0.0E+00;
        for ( UInt iBin = 0; iBin < nBin; iBin++ ) {
            vol = 0.5 * ( bin_VEdges[iBin] + bin_VEdges[iBin+1] );
            for ( UInt jNy = 0; jNy < Ny; jNy++ ) {
                for ( UInt iNx = 0; iNx < Nx; iNx++ )
                    bin_VCenters[iBin][jNy][iNx] = vol;
            }
        }
        
        pdf.resize( nBin, Vector_2D( Ny, Vector_1D( Nx, 0.0E+00 ) ) );

        /* Allocate number of particles */
        nPart = nPart_;

        /* Allocate mean and standard deviation */
        if ( mu_ <= 0 ) {
            std::cout << "\nIn Grid_Aerosol::Grid_Aerosol: mean/mode is negative: mu = " << mu_ << "\n";
        }
        mu    = mu_;
        if ( sigma_ <= 0 ) {
            std::cout << "\nIn Grid_Aerosol::Grid_Aerosol: standard deviation is negative: stddev = " << sigma_ << "\n";
        }
        sigma = sigma_;

        /* Allocate alpha (only used in power law and generalized gamma distributions) */
        alpha = alpha_;

        /* Initialize distribution */
        if ( ( strcmp( type, "log" ) == 0 ) || ( strcmp( type, "lognormal" ) == 0 ) ) {
            /* Log-normal distribution:
             * dn/d(ln(r)) = N / ( sqrt(2*\pi) * ln(sigma) ) * exp( - ( ln(r) - ln(r_m) ) ^ 2 / ( 2 ln(sigma) ^2 )) */

            if ( sigma <= 1.0 ) {
                std::cout << "\nIn Grid_Aerosol::Grid_Aerosol: log-normal distribution requires that stddev > 1.0 ( stddev = " << sigma << " )\n";
            }

            for ( UInt iBin = 0; iBin < bin_Centers.size(); iBin++ ) {
                pdf[iBin][0][0] = nPart * exp( - 0.5 * ( ( log( bin_Centers[iBin] ) - log( mu ) ) / log( sigma ) ) * ( ( log( bin_Centers[iBin] ) - log( mu ) ) / log( sigma ) ) ) / ( sqrt( 2.0 * physConst::PI ) * log( sigma ) );
                for ( UInt jNy = 0; jNy < Ny; jNy++ ) {
                    for ( UInt iNx = 0; iNx < Nx; iNx++ ) {
                        pdf[iBin][jNy][iNx] = pdf[iBin][0][0];
                    }
                }
            }
            
        } else if ( ( strcmp( type, "norm" ) == 0 ) || ( strcmp( type, "normal" ) == 0 ) ) {
            /* Normal distribution:
             * dn/d(ln(r)) = N / ( sqrt(2*\pi) * sigma ) * exp( - ( r - r_m ) ^ 2 / ( 2 * sigma ^2 )) */

            for ( UInt iBin = 0; iBin < bin_Centers.size(); iBin++ ) {
                pdf[iBin][0][0] = nPart * exp( - 0.5 * ( ( bin_Centers[iBin] - mu ) / sigma ) * ( ( bin_Centers[iBin] - mu ) / sigma ) ) / ( sqrt( 2.0 * physConst::PI ) * sigma );
                for ( UInt jNy = 0; jNy < Ny; jNy++ ) {
                    for ( UInt iNx = 0; iNx < Nx; iNx++ ) {
                        pdf[iBin][jNy][iNx] = pdf[iBin][0][0];
                    }
                }
            }

        } else if ( ( strcmp( type, "pow" ) == 0 ) || ( strcmp( type, "power" ) == 0 ) || ( strcmp( type, "junge" ) == 0 ) ) {
            /* Power distribution:
             * dn/d(ln(r)) = N * alpha * ( r / r_min ) ^ (-alpha) */

            if ( alpha <= 0.0 ) {
                std::cout << "\nIn Grid_Aerosol::Grid_Aerosol: power law requires that alpha > 0 ( alpha = " << alpha_ << " )\n";
            }
            
            for ( UInt iBin = 0; iBin < bin_Centers.size(); iBin++ ) {
                pdf[iBin][0][0] = nPart * alpha * pow( bin_Centers[iBin] / bin_Centers[0], -alpha );
                for ( UInt jNy = 0; jNy < Ny; jNy++ ) {
                    for ( UInt iNx = 0; iNx < Nx; iNx++ ) {
                        pdf[iBin][jNy][iNx] = pdf[iBin][0][0];
                    }
                }
            }

        } else if ( ( strcmp( type, "gam" ) == 0 ) || ( strcmp( type, "gamma" ) == 0 ) || ( strcmp( type, "generalized gamma" ) == 0 ) ) {
            /* Gamma distribution:
             * dn/d(ln(r)) = gamma * b ^ ((alpha + 1)/gamma) / Gamma((alpha+1)/gamma) * r ^ (alpha + 1) * exp( - b * r ^ (gamma) ) 
             * alpha and gamma are of the order of unity (1 - 10)
             * b must be of the order of r ^ (-gamma) >> 1 */

            if ( ( alpha <= 0 ) || (std::fmod( alpha, 1.0 ) != 0 ) ) {
                std::cout << "\nIn Grid_Aerosol::Grid_Aerosol: (generalized) gamma distribution requires that alpha is a positive integer ( alpha = " << alpha << " )\n";
            }
            if ( gamma_ <= 0 ) {
                std::cout << "\nIn Grid_Aerosol::Grid_Aerosol: (generalized) gamma distribution requires that gamma is positive ( gamma = " << gamma_ << " )\n";
            }
            if ( b_ <= 0 ) {
                std::cout << "\nIn Grid_Aerosol::Grid_Aerosol: (generalized) gamma distribution requires that b is positive ( b = " << b_ << " )\n";
            }
            
            for ( UInt iBin = 0; iBin < bin_Centers.size(); iBin++ ) {
                pdf[iBin][0][0] = nPart * gamma_ * pow( b_, (alpha+1) / gamma_ ) / boost::math::tgamma( (alpha+1) / gamma_ ) * pow( bin_Centers[iBin], alpha + 1 ) * exp( - b_ * pow( bin_Centers[iBin], gamma_ ) );
                for ( UInt jNy = 0; jNy < Ny; jNy++ ) {
                    for ( UInt iNx = 0; iNx < Nx; iNx++ ) {
                        pdf[iBin][jNy][iNx] = pdf[iBin][0][0];
                    }
                }
            }

        } else {
            std::cout << "\nIn Grid_Aerosol::Grid_Aerosol: distribution type must be either lognormal, normal, power or (generalized) gamma\n";
            std::cout << "\nCurrent type is " << type << "\n";
        }
      
//        /* Check that we get the right number of particles */
//        if ( ( std::abs( Moment() - nPart ) / nPart > 0.10 ) && ( nPart > 1.0E-10 )) {
//            std::cout << "\nIn Grid_Aerosol::Grid_Aerosol: the size range doesn't cover the full distribution";
//            std::cout << "\nPrescribed number: " << nPart << ", Number covered: " << Moment() << " [#/cm^3]\n";
//            if ( ( strcmp( type, "log" ) == 0 ) || ( strcmp( type, "lognormal" ) == 0 ) ) {
//                std::cout << "For lognormal distribution, prescribed mode is: " << mu << " [m]\n";
//            } else if ( ( strcmp( type, "norm" ) == 0 ) || ( strcmp( type, "normal" ) == 0 ) ) {
//                std::cout << "For normal distribution, prescribed mode is: " << mu << " [m]\n";
//            }
//        }
        
        /* For lognormal distributions:
         * rMode = rMean * exp( - ln(sigma)^2 );
         * rMedi = rMean * exp( 0.5 * ln(sigma)^2 ); */

    } /* End of Grid_Aerosol::Grid_Aerosol */

    Grid_Aerosol::Grid_Aerosol( const Grid_Aerosol &rhs )
    {

        bin_Centers = rhs.bin_Centers;
        bin_VCenters = rhs.bin_VCenters;
        bin_Edges = rhs.bin_Edges;
        bin_VEdges = rhs.bin_VEdges;
        bin_Sizes = rhs.bin_Sizes;
        nBin = rhs.nBin;
        nPart = rhs.nPart;
        type = rhs.type;
        mu = rhs.mu;
        sigma = rhs.sigma;
        alpha = rhs.alpha;
        pdf = rhs.pdf;

    } /* End of Grid_Aerosol::Grid_Aerosol */
    
    Grid_Aerosol::~Grid_Aerosol( )
    {

        /* Destructor */

    } /* End of Grid_Aerosol::~Grid_Aerosol */

    Grid_Aerosol& Grid_Aerosol::operator=( const Grid_Aerosol &rhs )
    {

        if ( &rhs == this )
            return *this;

        Nx = rhs.Nx;
        Ny = rhs.Ny;
        bin_Centers = rhs.bin_Centers;
        bin_VCenters = rhs.bin_VCenters;
        bin_Edges = rhs.bin_Edges;
        bin_VEdges = rhs.bin_VEdges;
        bin_Sizes = rhs.bin_Sizes;
        nBin = rhs.nBin;
        nPart = rhs.nPart;
        type = rhs.type;
        mu = rhs.mu;
        sigma = rhs.sigma;
        alpha = rhs.alpha;
        pdf = rhs.pdf;
        return *this;

    } /* End of Grid_Aerosol::operator= */

    Grid_Aerosol Grid_Aerosol::operator+=( const Grid_Aerosol &rhs )
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

        Vector_3D pdf_rhs = rhs.getPDF();
        for ( UInt iBin = 0; iBin < nBin; iBin++ ) {
            for ( UInt jNy = 0; jNy < Ny; jNy++ ) {
                for ( UInt iNx = 0; iNx < Nx; iNx++ ) {
                    pdf[iBin][jNy][iNx] += pdf_rhs[iBin][jNy][iNx];
                }
            }
        }

        return *this;

    } /* End of Grid_Aerosol::operator+= */
    
    Grid_Aerosol Grid_Aerosol::operator-=( const Grid_Aerosol &rhs )
    {

        if ( nBin != rhs.getNBin() ) {
            std::cout << "\nIn Grid_Aerosol::operator+=: aerosol distributions do not have the same number of bins: " << nBin << " != " << rhs.getNBin();
            return *this;
        }

        Vector_1D bin_Centers_rhs = rhs.getBinCenters();
        for ( UInt iBin = 0; iBin < nBin; iBin++ ) {
            if ( std::abs( bin_Centers[iBin] - bin_Centers_rhs[iBin] ) > 1.0E-10 ) {
                std::cout << "\nIn Grid_Aerosol::operator+=: aerosol distributions do not have the same bin centers!";
                return *this;
            }
        }

        Vector_3D pdf_rhs = rhs.getPDF();
        for ( UInt iBin = 0; iBin < nBin; iBin++ ) {
            for ( UInt jNy = 0; jNy < Ny; jNy++ ) {
                for ( UInt iNx = 0; iNx < Nx; iNx++ ) {
                    pdf[iBin][jNy][iNx] -= pdf_rhs[iBin][jNy][iNx];
                }
            }
        }

        return *this;

    } /* End of Grid_Aerosol::operator-= */
    
    Grid_Aerosol Grid_Aerosol::operator+( const Grid_Aerosol &rhs ) const
    {

        Grid_Aerosol result = *this;
        
        if ( nBin != rhs.getNBin() ) {
            std::cout << "\nIn Grid_Aerosol::operator+=: aerosol distributions do not have the same number of bins: " << nBin << " != " << rhs.getNBin();
            return *this;
        }

        Vector_1D bin_Centers_rhs = rhs.getBinCenters();
        for ( UInt iBin = 0; iBin < nBin; iBin++ ) {
            if ( std::abs( bin_Centers[iBin] - bin_Centers_rhs[iBin] ) > 1.0E-10 ) {
                std::cout << "\nIn Grid_Aerosol::operator+=: aerosol distributions do not have the same bin centers!";
                return *this;
            }
        }

        Vector_3D pdf_rhs = rhs.getPDF();
        for ( UInt iBin = 0; iBin < nBin; iBin++ ) {
            for ( UInt jNy = 0; jNy < Ny; jNy++ ) {
                for ( UInt iNx = 0; iNx < Nx; iNx++ ) {
                    result.pdf[iBin][jNy][iNx] += pdf_rhs[iBin][jNy][iNx];
                }
            }
        }
        
        return result;

    } /* End of Grid_Aerosol::operator+ */
    
    Grid_Aerosol Grid_Aerosol::operator-( const Grid_Aerosol &rhs ) const
    {

        Grid_Aerosol result = *this;
        
        if ( nBin != rhs.getNBin() ) {
            std::cout << "\nIn Grid_Aerosol::operator+=: aerosol distributions do not have the same number of bins: " << nBin << " != " << rhs.getNBin();
            return *this;
        }

        Vector_1D bin_Centers_rhs = rhs.getBinCenters();
        for ( UInt iBin = 0; iBin < nBin; iBin++ ) {
            if ( std::abs( bin_Centers[iBin] - bin_Centers_rhs[iBin] ) > 1.0E-10 ) {
                std::cout << "\nIn Grid_Aerosol::operator+=: aerosol distributions do not have the same bin centers!";
                return *this;
            }
        }

        Vector_3D pdf_rhs = rhs.getPDF();
        for ( UInt iBin = 0; iBin < nBin; iBin++ ) {
            for ( UInt jNy = 0; jNy < Ny; jNy++ ) {
                for ( UInt iNx = 0; iNx < Nx; iNx++ ) {
                    result.pdf[iBin][jNy][iNx] -= pdf_rhs[iBin][jNy][iNx];
                }
            }
        }
        
        return result;

    } /* End of Grid_Aerosol::operator- */
       
    void Grid_Aerosol::Coagulate( const RealDouble dt, Coagulation &kernel, const UInt N, const UInt SYM )
    {

        /* DESCRIPTION:
         * Performs self-coagulation. Updates Aerosol.pdf
         * The numerical scheme is taken from: 
         * M.Z. Jacobson, Fundamentals of atmospheric modeling. Cambridge university press, 2005.*/
         
        /* INPUT:
         * - RealDouble dt      :: Timestep in s
         * - Coagulation kernel :: Coagulation structure containing the coagulation kernels
         * - UInt N             :: Coagulation scenarios ( 0, 1 or 2 )
         * - UInt SYM           :: Symmetry?
         *
         * OUTPUT:
         *
         */
 
        /* For debug purposes */
        const bool checkMass = 0;
        if ( checkMass ) {
            std::cout << "Coagulation is a volume-conserving process. The two following quantities should be identical\n";
            std::cout << "At t     : " << Moment( 3, Nx/2, Ny/2 ) * 1.0E+18 << "[um^3/cm^3]" << std::endl;
        }

        UInt Nx_max, Ny_max;

        if ( N == 0 ) {
            /* No coagulation is performed */
            return;
        } else if ( N == 1 ) {
            /* No emitted aerosols -> Aerosol is a uniform field */
            /* Perform coagulation only once */
            Nx_max = 1;
            Ny_max = 1;
        } else if ( N == 2 ) {
            /* Perform coagulation for all */
            if ( SYM == 2 ) {
                /* Both X and Y symmetry */
                Nx_max = Nx/2;
                Ny_max = Ny/2;
            } else if ( SYM == 1 ) {
                /* Only symmetry around the Y axis */
                Nx_max = Nx/2;
                Ny_max = Ny;
            } else if ( SYM == 0 ) {
                /* No symmetry */
                Nx_max = Nx;
                Ny_max = Ny;
            } else {
                std::cout << " In Grid_Aerosol::Coagulate: Wrong input for SYM\n";
                std::cout << " SYM = " << SYM << "\n";
                return;
            }
        } else {
            std::cout << " In Grid_Aerosol::Coagulate: Wrong input for N\n";
            std::cout << " N = " << N << "\n";
            return;
        }

        /* Allocate variables */
        Vector_1D P( nBin, 0.0E+00 );
        Vector_1D L( nBin, 0.0E+00 );
        UInt iBin, jBin, kBin, kBin_;
        UInt jNy, iNx; /* Grid indices */

        /* Particle volume in each bin */
        Vector_3D v = Volume( ); /* Expressed in [m^3/cm^3] */
        /* Copy v into v_new */
        Vector_3D v_new = v;

        /* Total volume and number per grid cell */
        RealDouble totVol, nPart;

        /* Description of the algorithm:
         * \frac{dv}{dt}[iBin] = P - L * v[iBin] 
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
        for ( jNy = 0; jNy < Ny_max; jNy++ ) {
            for ( iNx = 0; iNx < Nx_max; iNx++ ) {

                /* Total aerosol volume */
                totVol = 0.0E+00;
                for ( iBin = 0; iBin < nBin; iBin++ )
                    totVol += v[iBin][jNy][iNx]; /* [m^3/cm^3] */

                if ( totVol * 1E18 > 0.1 ) {
                    /* Only run coagulation where aerosol volume is greater
                     * than 0.1 um^3/cm^3 */

                    /* Update kernel (updates kernel.f and kernel.indices) */
                    kernel.buildF( bin_VCenters, jNy, iNx );

                    for ( iBin = 0; iBin < nBin; iBin++ ) {

                        /* Reset P and L values */
                        P[iBin] = 0.0E+00;
                        L[iBin] = 0.0E+00;

                        /* Build production and loss terms */
                        for ( jBin = 0; jBin < nBin; jBin++ ) {

                            nPart = pdf[jBin][jNy][iNx] * log( bin_Edges[jBin+1] / bin_Edges[jBin] );

                            if ( jBin <= iBin ) {
                                for ( kBin_ = 0; kBin_ < kernel.indices[jBin][iBin].size(); kBin_++ ) {
                                    kBin = kernel.indices[jBin][iBin][kBin_];
                                    /* k coagulating with j to form i */
                                    if ( kBin < iBin ) {
                                        P[iBin] += kernel.f[kBin][jBin][iBin] * kernel.beta[kBin][jBin] * v_new[kBin][jNy][iNx] * nPart;
                                        /* [cm^3/#/s] * [m^3/cm^3] * [#/cm^3] = [m^3/cm^3/s] */
                                    }
                                }

                                /* The following lines perform slightly slower than the for loop above */
                                //for ( kBin = 0; kBin < iBin ; kBin++ ) {
                                //    /* k coagulating with j to form i */
                                //    if ( kernel.f[kBin][jBin][iBin] != 0.0 ) { /* f is a somewhat sparse 3D tensor */
                                //        P[iBin] += kernel.f[kBin][jBin][iBin] * kernel.beta[kBin][jBin] * v_new[kBin][jNy][iNx] * nPart;
                                //        /* [cm^3/#/s] * [m^3/cm^3] * [#/cm^3] = [m^3/cm^3/s] */
                                //    }
                                //}
                            }

                            /* i coagulating with j to deplete i */
                            if ( kernel.f[iBin][jBin][iBin] != 1.0 )
                                L[iBin] += ( 1.0 - kernel.f[iBin][jBin][iBin] ) * kernel.beta[iBin][jBin] * nPart;

                        }

                        /* Non-mass conserving scheme: */
                        //  v_new[iBin][jNy][iNx] = v[iBin][jNy][iNx] + ( P[iBin] - L[iBin] * v[iBin][jNy][iNx] ) * dt;

                        /* Mass conserving scheme: */
                        v_new[iBin][jNy][iNx] = ( v[iBin][jNy][iNx] + dt * P[iBin] ) / ( 1.0 + dt * L[iBin] );

                        if ( v[iBin][jNy][iNx] > 0.0E+00 )
                            pdf[iBin][jNy][iNx] *= v_new[iBin][jNy][iNx] / v[iBin][jNy][iNx];

                    }
                }
            }
        }

        /* Update bin centers */
        UpdateCenters( v_new, pdf );

        if ( checkMass )
            std::cout << "At t + dt: " << Moment( 3, Nx/2, Ny/2 ) * 1.0E+18 << "[um^3/cm^3]" << std::endl;

        if ( N == 1 ) {
            /* Allocate uniform results to the grid */
            for ( jNy = 0; jNy < Ny; jNy++ ) {
                for ( iNx = 0; iNx < Nx; iNx++ ) {
                    for ( iBin = 0; iBin < nBin; iBin++ ) {
                        pdf[iBin][jNy][iNx] = pdf[iBin][0][0];
                        bin_VCenters[iBin][jNy][iNx] = bin_VCenters[iBin][0][0];
                    }
                }
            }
        } else if ( N == 2 ) {
            /* Apply symmetry */
            if ( SYM == 2 ) {
                for ( iBin = 0; iBin < nBin; iBin++ ) {
                    for ( jNy = 0; jNy < Ny; jNy++ ) {
                        for ( iNx = Nx_max; iNx < Nx; iNx++ ) {
                            pdf[iBin][jNy][iNx] = pdf[iBin][jNy][Nx-1-iNx];
                            bin_VCenters[iBin][jNy][iNx] = bin_VCenters[iBin][jNy][Nx-1-iNx];
                        }
                    }
                    for ( jNy = Ny_max; jNy < Ny; jNy++ ) {
                        for ( iNx = 0; iNx < Nx; iNx++ ) {
                            pdf[iBin][jNy][iNx] = pdf[iBin][Ny-1-jNy][iNx];
                            bin_VCenters[iBin][jNy][iNx] = bin_VCenters[iBin][Ny-1-jNy][iNx];
                        }
                    }
                    for ( jNy = Ny_max; jNy < Ny; jNy++ ) {
                        for ( iNx = Nx_max; iNx < Nx; iNx++ ) {
                            pdf[iBin][jNy][iNx] = pdf[iBin][Ny-1-jNy][Nx-1-iNx];
                            bin_VCenters[iBin][jNy][iNx] = bin_VCenters[iBin][Ny-1-jNy][Nx-1-iNx];
                        }
                    }
                }
            } else if ( SYM == 1 ) {
                for ( iBin = 0; iBin < nBin; iBin++ ) {
                    for ( jNy = 0; jNy < Ny; jNy++ ) {
                        for ( iNx = Nx_max; iNx < Nx; iNx++ ) {
                            pdf[iBin][jNy][iNx] = pdf[iBin][jNy][Nx-1-iNx];
                            bin_VCenters[iBin][jNy][iNx] = bin_VCenters[iBin][jNy][Nx-1-iNx];
                        }
                    }
                }
            }
        }

    } /* End of Grid_Aerosol::Coagulate */

    void Grid_Aerosol::Grow( const RealDouble dt, Vector_2D &H2O, const Vector_2D &T, const Vector_1D &P, const UInt N, const UInt SYM )
    {

        /* DESCRIPTION:
         * Computes growth of ice crystals through direct ice deposition. Updates Aerosol.pdf and H2O
         * The numerical scheme is taken from:
         * M.Z. Jacobson, Fundamentals of atmospheric modeling. Cambridge university press, 2005.*/

        /* INPUT:
         * - RealDouble dt :: Timestep in s
         * - Vector_2D H2O :: Vector containing water vapor molecular concentrations [molec/cm^3]
         *    -> ( Ny x Nx )
         * - Vector_2D T   :: Vector containing temperature values [K]
         *    -> ( Ny x Nx )
         * - Vector_1D P   :: Vector containing pressure values [Pa]
         *    -> ( Ny )
         * - UInt N        :: Growth scenarios ( 0, 1 or 2 )
         * - UInt SYM      :: Symmetry?
         *
         * OUTPUT:
         *
         */

        UInt Nx_max, Ny_max;

        if ( N == 0 ) {
            /* No growth is performed */
            return;
        } else if ( N == 1 ) {
            /* No emitted aerosols -> Aerosol is a uniform field */
            /* Perform growth only once */
            Nx_max = 1;
            Ny_max = 1;
        } else if ( N == 2 ) {
            /* Perform growth for all */
            if ( SYM == 2 ) {
                /* Both X and Y symmetry. Never happens if settling is on!! */
                Nx_max = Nx/2;
                Ny_max = Ny/2;
            } else if ( SYM == 1 ) {
                /* Only symmetry around the Y axis. Never happens if shear is on!! */
                Nx_max = Nx/2;
                Ny_max = Ny;
            } else if ( SYM == 0 ) {
                /* No symmetry */
                Nx_max = Nx;
                Ny_max = Ny;
            } else {
                std::cout << " In Grid_Aerosol::Grow: Wrong input for SYM\n";
                std::cout << " SYM = " << SYM << "\n";
                return;
            }
        } else {
            std::cout << " In Grid_Aerosol::Grow: Wrong input for N\n";
            std::cout << " N = " << N << "\n";
            return;
        }

        /* Minimum, maximum particle volumes */
        const RealDouble MINVOL = bin_VEdges[0];
        const RealDouble MAXVOL = bin_VEdges[nBin];

        /* Conversion factor from ice volume [m^3] to [molecules] */ 
        const RealDouble UNITCONVERSION = physConst::RHO_ICE / MW_H2O * physConst::Na;
        /* Unit check: [kg/m^3] / [kg/mol] * [molec/mol] = [molec/m^3] */

        /* Scaled Boltzmann constant */
        const RealDouble kB_ = physConst::kB * 1.00E+06;

        /* Declare and initialize particle totals and water vapor array */
        Vector_3D icePart = Number( );
        Vector_3D iceVol  = Volume( );
        Vector_2D totH2O  = H2O;

        RealDouble partVol  = 0.0E+00;
        RealDouble icePart_ = 0.0E+00;
        RealDouble iceVol_  = 0.0E+00;

        int jBin = -1;
        std::vector<int> toBin( nBin, 0 );
        std::vector<int>::iterator iterBegin, iterCurr, iterEnd;

        /* Declare and initialize variable to store saturation quantities,
         * pressure and temperature */
        RealDouble pSat = 0.0E+00;
        RealDouble nSat = 0.0E+00;
        RealDouble locT = 0.0E+00;
        RealDouble locP = 0.0E+00;
        /* Declare and initialize total ice concentration */
        RealDouble totH2Oi = 0.0E+00;
        /* Vector containing Kelvin factors evaluated at each bin center */
        Vector_1D kFactor( nBin, 0.0E+00 );

        /* Declare and initialize growth rates per bin */
        Vector_1D kGrowth( nBin, 0.0E+00 );
        /* Declare and initialize aggregated growth rates */
        RealDouble totkGrowth_1 = 0.0E+00;
        RealDouble totkGrowth_2 = 0.0E+00;

        /* Compute Kelvin factor */
        for ( UInt iBin = 0; iBin < nBin; iBin++ )
            kFactor[iBin] = physFunc::Kelvin( bin_Centers[iBin] );

        for ( UInt jNy = 0; jNy < Ny_max; jNy++ ) {
            for ( UInt iNx = 0; iNx < Nx_max; iNx++ ) {
                for ( UInt iBin = 0; iBin < nBin; iBin++ ) {
                    totH2O[jNy][iNx] += iceVol[iBin][jNy][iNx] * UNITCONVERSION;
                    /* Unit check:
                     * [ molec/cm^3 ] = [ m^3 ice/cm^3 air ]   * [ molec/m^3 ice ] */
                }
            }
        }

        for ( UInt jNy = 0; jNy < Ny_max; jNy++ ) {

            /* Store local pressure.
             * TODO: 
             * That might be moved into the loop over iNx eventually to 
             * account for 2D pressure met-fields?? */
            locP = P[jNy];

            for ( UInt iNx = 0; iNx < Nx; iNx++ ) {

                /* Reinitialize total rate and concentrations */
                totkGrowth_1 = 0.0E+00;
                totkGrowth_2 = 0.0E+00;
                totH2Oi      = 0.0E+00;

                /* Store local temperature */
                locT = T[jNy][iNx];

                /* Store local saturation pressure w.r.t ice */
                pSat = physFunc::pSat_H2Os( locT );

                if ( H2O[jNy][iNx] * kB_ * locT / pSat > 0.0 ) {

                    /* ================================================= */
                    /* ================================================= */
                    /* -------------------------------------------------
                     * Analytical predictor of condensation (APC) scheme
                     * -------------------------------------------------
                     * Mark Z. Jacobson, (1997), Numerical Techniques to
                     * Solve Condensational and Dissolutional Growth
                     * Equations When Growth is Coupled to Reversible
                     * Reactions, Aerosol Science and Technology,
                     * 27:4, 491-498, DOI: 10.1080/02786829708965489     */
                    /* ================================================= */
                    /* ================================================= */

                    /* APC scheme:
                     * dc_{i}(t)/dt = k_{i}(t) * (C(t) - S'_{i}(t) * C_{s,i}(t))       (1)
                     * dC(t)/dt     = -\sum k_{i} * (C(t) - S'_{i}(t) * C_{s,i}(t))    (2)
                     *
                     * The noniterative solution to the growth equation is
                     * obtained by integrating (1) for a final aerosol
                     * concentration.
                     * c_{i}(t) = c_{i}(t-dt) + ...
                     *          dt * k_{i}(t-h) * (C(t) - S'_{i}(t-dt) * C_{s,i}(t-dt) (3)
                     * where the final gas molar concentration C(t) is 
                     * currently unknown.
                     *
                     * Final aerosol and gas concentrations are constrained
                     * by the mass-balance equation:
                     * C(t) + \sum c_{i}(t) = C(t-dt) + \sum c_{i}(t-dt} = C_{tot}
                     *
                     * Solving for the gas concentration give
                     *          C(t-dt) + dt \sum k_{i}(t) S'_{i}(t) C_{s,i}(t)
                     * C(t) = --------------------------------------------------
                     *             1.0  + dt \sum k_{i}(t)
                     *
                     * The concentration from this equation cannot fall 
                     * below zero, but can increase above the total mass
                     * of the species in the system. In such cases, gas 
                     * concentration, C(t), is limited by 
                     * C(t) = min(C(t),C_{tot})
                     *
                     * Molar aerosol concentrations are determined by 
                     * plugging the obtained C(t) back into Equation (3).
                     *
                     * The gaseous molar concentration is then updated
                     * according to
                     *
                     * C(t) = C_{tot} - \sum c_{i}
                     *
                     * The APC scheme is unconditionally stable, since all
                     * final concentrations are bounded between 0 and 
                     * C_{tot}, independently of the time step */

                    /* Compute particle growth rates through ice deposition
                     * We here assume that C_{s,i} is independent of the
                     * bin and thus the particle size and only depends
                     * on meteorological parameters. */
                    for ( UInt iBin = 0; iBin < nBin; iBin++ ) {

                        /* kGrowth is expressed in [cm^3 ice/s/part] */
                        kGrowth[iBin] = physFunc::growthRate( bin_Centers[iBin], locT, locP, H2O[jNy][iNx] );

                        /* kGrowth_* are thus in 
                         * [(cm^3 ice/s)/cm^3 air] = [1/s] */
                        totkGrowth_1 += kGrowth[iBin] * icePart[iBin][jNy][iNx] \
                                        * kFactor[iBin];
                        totkGrowth_2 += kGrowth[iBin] * icePart[iBin][jNy][iNx];
                    }

                    /* Compute the molecular saturation concentration 
                     * C_{s,i} in [molec/cm^3] */
                    nSat = pSat / ( kB_ * locT );

                    /* Update gaseous molecular concentration */
                    H2O[jNy][iNx] = ( H2O[jNy][iNx] + dt * totkGrowth_1 * nSat ) \
                                  / (   1.00E+00    + dt * totkGrowth_2        );

                    /* Make sure that molecular water does not go over 
                     * total water (gaseous + solid) concentrations */
                    H2O[jNy][iNx] = std::min( H2O[jNy][iNx], totH2O[jNy][iNx] );

                    for ( UInt iBin = 0; iBin < nBin; iBin++ ) {
                        iceVol[iBin][jNy][iNx] += dt * kGrowth[iBin] * icePart[iBin][jNy][iNx] \
                                                * ( H2O[jNy][iNx] -  kFactor[iBin] * nSat ) / UNITCONVERSION;
                        /* Unit check:
                         * [m^3 ice/cm^3 air]   = [s] * [cm^3 ice/s/part] * [part/cm^3 air] \
                         *                      * [molec/cm^3 air] * [m^3 ice/molec] 
                         *                      = [cm^3 ice/cm^3 air] * [m^3 ice/cm^3 air] 
                         *                      = [m^3 ice/cm^3 air] */

                        iceVol[iBin][jNy][iNx] = \
                                std::min( std::max( iceVol[iBin][jNy][iNx], 0.0E+00 ), icePart[iBin][jNy][iNx] * MAXVOL );

                        /* Compute total water taken up on particles */
                        totH2Oi += iceVol[iBin][jNy][iNx] * UNITCONVERSION;
                        /* Unit check:
                         * [molec/cm^3 air] = [m^3 ice/cm^3 air] * [molec/m^3 ice] */
                    }

                    H2O[jNy][iNx] = totH2O[jNy][iNx] - totH2Oi; 

                }

                /* ======================================================= */
                /* ======================================================= */
                /* ============== Moving-center structure ================ */
                /* ======================================================= */
                /* ============= Update bin center average =============== */
                /* ======================================================= */
                /* ======================================================= */

                /* 1. Compute bin particle flux */

                for ( UInt iBin = 0; iBin < nBin; iBin++ ) {

                    /* What does bin iBin grow into? */
                    toBin[iBin] = -1;

                    /* Only perform computation if number of particle is
                     * greater than some small number to avoid division
                     * by ridiculously small numbers */

                    /* Compute particle volume */
                    partVol = iceVol[iBin][jNy][iNx] / icePart[iBin][jNy][iNx];

                    /* Find which bin corresponds to this particle 
                     * volume */
                    toBin[iBin] = std::lower_bound( bin_VEdges.begin(), bin_VEdges.end(), partVol ) \
                                  - bin_VEdges.begin() - 1;

                    if ( toBin[iBin] == 0 ) {
                        if ( partVol < bin_VEdges[0] )
                            /* Particles are reduced to their core
                             * and thus considered lost */
                            toBin[iBin] = -1;
                    }

                }

                /* 2. Attribute new particles according to fluxes */

                for ( UInt iBin = 0; iBin < nBin; iBin++ ) {

                    /* Find all bins that end up in bin iBin after growth */

                    /* Initialize total new number of particles and volume 
                     * to 0 */
                    icePart_ = 0.0E+00;
                    iceVol_  = 0.0E+00;

                    /* Bin jBin -> Bin iBin */
                    jBin = -1;

                    /* Initialize iterators */
                    iterBegin = toBin.begin();
                    iterEnd   = toBin.end();
                    iterCurr  = iterBegin;

                    while ((iterCurr = std::find(iterCurr, iterEnd, iBin)) != iterEnd) {
                        jBin = iterCurr - iterBegin; //std::distance(iterBegin, iterCurr);

                        /* If jBin -> iBin, then add particle number and 
                         * volume to sum */

                        icePart_ += icePart[jBin][jNy][iNx];
                        iceVol_  += iceVol[jBin][jNy][iNx];

                        /* Iterate */
                        iterCurr++;
                    }

                    if ( icePart_ > 0.0E+00 ) {
                        /* Bin is not empty */

                        /* Compute particle volume:
                         * [m^3] = [m^3/cm^3 air] / [#/cm^3 air] 
                         * and clip it between min and max volume allowed. */

                        bin_VCenters[iBin][jNy][iNx] = std::max( std::min( iceVol_ / icePart_, bin_VEdges[iBin+1] ), bin_VEdges[iBin] );

                        pdf[iBin][jNy][iNx] = icePart_ / ( log( bin_Edges[iBin+1] / bin_Edges[iBin] ) );

                    } else {
                        /* Bin is empty */

                        /* Set bin center to average volume of the bin.
                         * This arbitrary value should not matter because
                         * no particles are in this bin */

                        bin_VCenters[iBin][jNy][iNx] = 0.5 * ( bin_VEdges[iBin] + bin_VEdges[iBin+1] );
                        pdf[iBin][jNy][iNx] = 0.0E+00;

                    }
                }
            }
        }

        if ( N == 1 ) {
            /* Allocate uniform results to the grid */
            for ( UInt jNy = 0; jNy < Ny; jNy++ ) {
                for ( UInt iNx = 0; iNx < Nx; iNx++ ) {
                    for ( UInt iBin = 0; iBin < nBin; iBin++ ) {
                        pdf[iBin][jNy][iNx] = pdf[iBin][0][0];
                        H2O[jNy][iNx] = H2O[0][0];
                        bin_VCenters[iBin][jNy][iNx] = bin_VCenters[iBin][0][0];
                    }
                }
            }
        } else if ( N == 2 ) {
            /* Apply symmetry */
            if ( SYM == 2 ) {
                /* Symmetry around the origin */
                for ( UInt iBin = 0; iBin < nBin; iBin++ ) {
                    for ( UInt jNy = 0; jNy < Ny; jNy++ ) {
                        for ( UInt iNx = Nx_max; iNx < Nx; iNx++ ) {
                            pdf[iBin][jNy][iNx] = pdf[iBin][jNy][Nx-1-iNx];
                            H2O[jNy][iNx] = H2O[jNy][Nx-1-iNx];
                            bin_VCenters[iBin][jNy][iNx] = bin_VCenters[iBin][jNy][Nx-1-iNx];
                        }
                    }
                    for ( UInt jNy = Ny_max; jNy < Ny; jNy++ ) {
                        for ( UInt iNx = 0; iNx < Nx; iNx++ ) {
                            pdf[iBin][jNy][iNx] = pdf[iBin][Ny-1-jNy][iNx];
                            H2O[jNy][iNx] = H2O[Ny-1-jNy][iNx];
                            bin_VCenters[iBin][jNy][iNx] = bin_VCenters[iBin][Ny-1-jNy][iNx];
                        }
                    }
                    for ( UInt jNy = Ny_max; jNy < Ny; jNy++ ) {
                        for ( UInt iNx = Nx_max; iNx < Nx; iNx++ ) {
                            pdf[iBin][jNy][iNx] = pdf[iBin][Ny-1-jNy][Nx-1-iNx];
                            H2O[jNy][iNx] = H2O[Ny-1-jNy][Nx-1-iNx];
                            bin_VCenters[iBin][jNy][iNx] = bin_VCenters[iBin][Ny-1-jNy][Nx-1-iNx];
                        }
                    }
                }
            } else if ( SYM == 1 ) {
                /* Symmetry around the Y-axis */
                for ( UInt iBin = 0; iBin < nBin; iBin++ ) {
                    for ( UInt jNy = 0; jNy < Ny; jNy++ ) {
                        for ( UInt iNx = Nx_max; iNx < Nx; iNx++ ) {
                            pdf[iBin][jNy][iNx] = pdf[iBin][jNy][Nx-1-iNx];
                            H2O[jNy][iNx] = H2O[jNy][Nx-1-iNx];
                            bin_VCenters[iBin][jNy][iNx] = bin_VCenters[iBin][jNy][Nx-1-iNx];
                        }
                    }
                }
            }
        }

    } /* End of Grid::Aerosol::Grow */

    void Grid_Aerosol::UpdateCenters( const Vector_3D &iceV, const Vector_3D &PDF ) {

        const RealDouble TINY = 1.00E-50;

        RealDouble ratio = 0.0E+00;

        for ( UInt iBin = 0; iBin < nBin; iBin++ ) {
            ratio = log( bin_Edges[iBin+1] / bin_Edges[iBin] );
            for ( UInt jNy = 0; jNy < Ny; jNy++ ) {
                for ( UInt iNx = 0; iNx < Nx; iNx++ ) {
                    if ( PDF[iBin][jNy][iNx] > TINY )
                        bin_VCenters[iBin][jNy][iNx] = std::max( std::min( iceV[iBin][jNy][iNx] / PDF[iBin][jNy][iNx] / ratio, \
                                                                           0.9999 * bin_VEdges[iBin+1] ), \
                                                                 1.0001 * bin_VEdges[iBin] );
                    else
                        bin_VCenters[iBin][jNy][iNx] = 0.5 * ( bin_VEdges[iBin] + bin_VEdges[iBin+1] );
                }
            }
        }

    } /* End of Grid_Aerosol::UpdateCenters */

    Vector_2D Grid_Aerosol::Moment( UInt n ) const
    {

        Vector_2D moment( Ny, Vector_1D( Nx, 0.0E+00 ) );
        const RealDouble FACTOR = 3.0 / RealDouble( 4.0 * physConst::PI );

        for ( UInt jNy = 0; jNy < Ny; jNy++ ) {
            for ( UInt iNx = 0; iNx < Nx; iNx++ ) {
                for ( UInt iBin = 0; iBin < nBin; iBin++ ) {
                    moment[jNy][iNx] += ( log( bin_Edges[iBin+1] / bin_Edges[iBin] ) ) * pow( FACTOR * bin_VCenters[iBin][jNy][iNx], n / RealDouble( 3.0 ) ) * pdf[iBin][jNy][iNx];
                }
            }
        }

        return moment;

    } /* End of Grid_Aerosol::Moment */

    Vector_3D Grid_Aerosol::Number( ) const
    {

        Vector_3D number( nBin, Vector_2D( Ny, Vector_1D( Nx, 0.0E+00 ) ) );
        RealDouble ratio = 0.0E+00;

        for ( UInt iBin = 0; iBin < nBin; iBin++ ) {
            ratio = log( bin_Edges[iBin+1] / bin_Edges[iBin] );
            for ( UInt jNy = 0; jNy < Ny; jNy++ ) {
                for ( UInt iNx = 0; iNx < Nx; iNx++ ) {
                    number[iBin][jNy][iNx] = ratio * pdf[iBin][jNy][iNx];
                    /* Unit check: [#/cm^3] */
                }
            }
        }

        return number;

    } /* End of Grid_Aerosol::Number */

    Vector_2D Grid_Aerosol::TotalNumber( ) const
    {

        return Moment( 0 );

    } /* End of Grid_Aerosol::TotalNumber */

    Vector_3D Grid_Aerosol::Volume( ) const
    {

        Vector_3D volume( nBin, Vector_2D( Ny, Vector_1D( Nx, 0.0E+00 ) ) );
        RealDouble ratio = 0.0E+00;

        for ( UInt iBin = 0; iBin < nBin; iBin++ ) {
            ratio = log( bin_Edges[iBin+1] / bin_Edges[iBin] );
            for ( UInt jNy = 0; jNy < Ny; jNy++ ) {
                for ( UInt iNx = 0; iNx < Nx; iNx++ ) {
                    volume[iBin][jNy][iNx] = ratio * bin_VCenters[iBin][jNy][iNx] * pdf[iBin][jNy][iNx];
                    /* Unit check:                   [m^3] * [#/cm^3] = [m^3/cm^3] */
                }
            }
        }

        return volume;

    } /* End of Grid_Aerosol::Volume */
 
    Vector_2D Grid_Aerosol::TotalVolume( ) const
    {

        Vector_2D m3 = Moment( 3 );
        const RealDouble FACTOR = 4.0 / RealDouble(3.0) * physConst::PI;

        /* V = 4.0/3.0*pi*m3 */
        for ( UInt jNy = 0; jNy < Ny; jNy++ ) {
            for ( UInt iNx = 0; iNx < Nx; iNx++ )
                m3[jNy][iNx] = m3[jNy][iNx] * FACTOR;
        }

        return m3;

    } /* End of Grid_Aerosol::TotalVolume */

    Vector_2D Grid_Aerosol::IWC( ) const
    {

        Vector_2D TVol = TotalVolume();
        const RealDouble FACTOR = physConst::RHO_ICE * 1.0E+06;

        /* Unit check: [m^3/cm^3] * [kg/m^3] * [cm^3/m^3] = [kg/m^3] */
        for ( UInt jNy = 0; jNy < Ny; jNy++ ) {
            for ( UInt iNx = 0; iNx < Nx; iNx++ )
                TVol[jNy][iNx] = TVol[jNy][iNx] * FACTOR;
        }

        return TVol;

    } /* End of Grid_Aerosol::IWC */

    Vector_2D Grid_Aerosol::Extinction( ) const
    {

        Vector_2D chi = IWC();
        Vector_2D rE  = EffRadius();

        const RealDouble a = 3.448E+00; /* [m^2/kg] */
        const RealDouble b = 2.431E-03; /* [m^3/kg] */

        for ( UInt jNy = 0; jNy < Ny; jNy++ ) {
            for ( UInt iNx = 0; iNx < Nx; iNx++ ) {
                if ( rE[jNy][iNx] > 1.00E-15 ) {
                    chi[jNy][iNx] = chi[jNy][iNx] * ( a + b / rE[jNy][iNx] );
                    /* Unit check: 
                     * [1/m] = [kg/m^3] * [m^2/kg] */
                } else
                    chi[jNy][iNx] = 0.0E+00;
            }
        }

        return chi;

    } /* End of Grid_Aerosol::Extinction */

    Vector_1D Grid_Aerosol::xOD( const Vector_1D xE ) const
    {

        Vector_1D tau_x( Ny, 0.0E+00 );
        Vector_2D chi = Extinction();

        for ( UInt jNy = 0; jNy < Ny; jNy++ ) {
            for ( UInt iNx = 0; iNx < Nx; iNx++ )
                tau_x[jNy] += ( xE[iNx+1] - xE[iNx] ) * chi[jNy][iNx];
        }

        return tau_x;

    } /* End of Grid_Aerosol::tau_x */

    Vector_1D Grid_Aerosol::yOD( const Vector_1D yE ) const
    {

        Vector_1D tau_y( Nx, 0.0E+00 );
        Vector_2D chi = Extinction();

        for ( UInt iNx = 0; iNx < Nx; iNx++ ) {
            for ( UInt jNy = 0; jNy < Ny; jNy++ )
                tau_y[iNx] += ( yE[jNy+1] - yE[jNy] ) * chi[jNy][iNx];
        }

        return tau_y;

    } /* End of Grid_Aerosol::tau_y */

    Vector_2D Grid_Aerosol::Radius( ) const
    {

        Vector_2D r( Ny, Vector_1D( Nx, 0.0E+00 ) );

        const Vector_2D m0 = Moment( 0 );
        const Vector_2D m1 = Moment( 1 );

        for ( UInt jNy = 0; jNy < Ny; jNy++ ) {
            for ( UInt iNx = 0; iNx < Nx; iNx++ ) {
                if ( m0[jNy][iNx] > 1.00E-50 )
                    r[jNy][iNx] = m1[jNy][iNx]/m0[jNy][iNx];
                else
                    r[jNy][iNx] = 0.0E+00;
            }
        }

        return r;

    } /* End of Grid_Aerosol::Radius */
    
    Vector_2D Grid_Aerosol::EffRadius( ) const
    {

        Vector_2D r_eff( Ny, Vector_1D( Nx, 0.0E+00 ) );

        const Vector_2D m2 = Moment( 2 );
        const Vector_2D m3 = Moment( 3 );

        for ( UInt jNy = 0; jNy < Ny; jNy++ ) {
            for ( UInt iNx = 0; iNx < Nx; iNx++ ) {
                if ( m2[jNy][iNx] > 1.00E-50 )
                    r_eff[jNy][iNx] = m3[jNy][iNx]/m2[jNy][iNx];
                else
                    r_eff[jNy][iNx] = 0.0E+00;
            }
        }

        return r_eff;

    } /* End of Grid_Aerosol::EffRadius */

    Vector_2D Grid_Aerosol::StdDev( ) const
    {

        Vector_2D sigma( Ny, Vector_1D( Nx, 0.0E+00 ) );

        const Vector_2D m0 = Moment( 0 );
        const Vector_2D m1 = Moment( 1 ); 
        const Vector_2D m2 = Moment( 2 ); 

        for ( UInt jNy = 0; jNy < Ny; jNy++ ) {
            for ( UInt iNx = 0; iNx < Nx; iNx++ ) {
                if ( m0[jNy][iNx] > 1.00E-50 )
                    sigma[jNy][iNx] = sqrt( m2[jNy][iNx] / m0[jNy][iNx] - m1[jNy][iNx] * m1[jNy][iNx] / ( m0[jNy][iNx] * m0[jNy][iNx] ) );
                else
                    sigma[jNy][iNx] = 0.0E+00;
            }
        }

        return sigma;

    } /* End of Grid_Aerosol::StdDev */
    
    RealDouble Grid_Aerosol::Moment( UInt n, Vector_1D PDF ) const
    {

        RealDouble moment = 0;

        for ( UInt iBin = 0; iBin < nBin; iBin++ ) {
            moment += ( log( bin_Edges[iBin+1] ) - log( bin_Edges[iBin] ) ) * pow( bin_Centers[iBin], n ) * PDF[iBin];
        }

        return moment;

    } /* End of Grid_Aerosol::Moment */
    
    RealDouble Grid_Aerosol::Moment( UInt n, UInt jNy, UInt iNx ) const
    {

        RealDouble moment = 0;
        const RealDouble FACTOR = 3.0 / RealDouble( 4.0 * physConst::PI );

        for ( UInt iBin = 0; iBin < nBin; iBin++ )
            moment += ( log( bin_Edges[iBin+1] / bin_Edges[iBin] ) ) * pow( FACTOR * bin_VCenters[iBin][jNy][iNx], n / RealDouble( 3.0 ) ) * pdf[iBin][jNy][iNx];

        return moment;

    } /* End of Grid_Aerosol::Moment */

    RealDouble Grid_Aerosol::Radius( UInt jNy, UInt iNx ) const
    {

        const RealDouble N = Moment( 0, jNy, iNx );

        if ( N > 0 ) {
            return Moment( 1, jNy, iNx ) / N;
        } else {
            //std::cout << "\nIn Grid_Aerosol::Radius: Number of particles is " << N << " <= 0\n";
            return 0;
        }

    } /* End of Grid_Aerosol::Radius */
    
    RealDouble Grid_Aerosol::EffRadius( UInt jNy, UInt iNx ) const
    {

        const RealDouble m2 = Moment( 2, jNy, iNx );

        if ( m2 > 0 ) {
            return Moment( 3, jNy, iNx ) / m2;
        } else {
            //std::cout << "\nIn Grid_Aerosol::EffRadius: Second moment is " << m2 << " <= 0\n";
            return 0;
        }

    } /* End of Aerosol::EffRadius */

    RealDouble Grid_Aerosol::StdDev( UInt jNy, UInt iNx ) const
    {

        const RealDouble N = Moment( 0, jNy, iNx );

        if ( N > 0 ) {
            return sqrt( Moment( 2, jNy, iNx ) / N - pow( Moment( 1, jNy, iNx ) / N, 2.0 ) );
        } else {
            //std::cout << "\nIn Grid_Aerosol::StdDev: Number of particles is " << N << " <= 0\n";
            return 0;
        }

    } /* End of Grid_Aerosol::StdDev */

    void Grid_Aerosol::updatePdf( Vector_3D pdf_ )
    {

        if ( pdf_.size() != pdf.size() ) {
            std::cout << "\nIn Grid_Aerosol::updatePdf:: sizes differ ( " << pdf_.size() << " != " << pdf.size() << " )\n";
        }
        else {
            pdf = pdf_;
        }

    } /* End of Grid_Aerosol::updatePdf */

    Vector_1D Grid_Aerosol::Average( const Vector_2D &weights,   \
                                     const RealDouble &totalWeight ) const
    {

        Vector_1D out( 4, 0.0E+00 );

        Vector_1D PDF( nBin, 0.0E+00 );

        unsigned int iBin, iNx, jNy;
        RealDouble w = 0.0E+00;

        for ( jNy = 0; jNy < Ny; jNy++ ) {
            for ( iNx = 0; iNx < Nx; iNx++ ) {
                w = weights[jNy][iNx] / totalWeight;
                for ( iBin = 0; iBin < nBin; iBin++ )
                    PDF[iBin] += pdf[iBin][jNy][iNx] * w;
            }
        }

        out[0] = Moment( 0, PDF );
        if ( out[0] > 1.00E-50 ) {
            out[2] = 4.0 * physConst::PI * Moment( 2, PDF );
            out[3] = 4.0 / 3.0 * physConst::PI * Moment( 3, PDF );
            out[1] = 3.0 * out[3] / out[2];
        } else {
            out[1] = 1.00E-07;
            out[2] = 0.00E+00;
            out[3] = 0.00E+00;
        }

        return out;

    } /* End of Grid_Aerosol::Average */

    void Grid_Aerosol::addPDF( const Aerosol &aerosol, const Vector_2D &weights )
    {

        Vector_1D AerPDF = aerosol.getPDF();

        unsigned int iBin, iNx, jNy;

        for ( jNy = 0; jNy < Ny; jNy++ ) {
            for ( iNx = 0; iNx < Nx; iNx++ ) {
                if ( weights[jNy][iNx] != 0.0E+00 ) {
                    for ( iBin = 0; iBin < nBin; iBin++ )
                        pdf[iBin][jNy][iNx] += AerPDF[iBin];
                }
            }
        }

    } /* End of Grid_Aerosol::addPDF */

    void Grid_Aerosol::addPDF( const Vector_1D &PDF, const Vector_2D &weights ) 
    {

        unsigned int iBin, iNx, jNy;

        for ( jNy = 0; jNy < Ny; jNy++ ) {
            for ( iNx = 0; iNx < Nx; iNx++ ) {
                if ( weights[jNy][iNx] != 0.0E+00 ) {
                    for ( iBin = 0; iBin < nBin; iBin++ )
                        pdf[iBin][jNy][iNx] += PDF[iBin];
                }
            }
        }

    } /* End of Grid_Aerosol::addPDF */

    Vector_1D Grid_Aerosol::getBinCenters() const
    {

        return bin_Centers;

    } /* End of Grid_Aerosol::getBinCenters */
    
    Vector_3D Grid_Aerosol::getBinVCenters() const
    {

        return bin_VCenters;

    } /* End of Grid_Aerosol::getBinVCenters */
        
    Vector_1D Grid_Aerosol::getBinEdges() const
    {
    
        return bin_Edges;

    } /* End of Grid_Aerosol::getBinEdges */
    
    Vector_1D Grid_Aerosol::getBinSizes() const
    {
    
        return bin_Sizes;

    } /* End of Grid_Aerosol::getBinSizes */
    
    UInt Grid_Aerosol::getNBin() const
    {

        return nBin;

    } /* End of Grid_Aerosol::getNBin */

    Vector_3D Grid_Aerosol::getPDF() const
    {

        return pdf;

    } /* End of Grid_Aerosol::getPDF */

}

/* End of Aerosol.cpp */
