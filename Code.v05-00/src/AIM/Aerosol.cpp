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
#include "Core/Structure.hpp"

namespace AIM
{

    Aerosol::Aerosol() : bin_Centers(2),
                         bin_VCenters(2),
                         bin_Edges(3),
                         bin_Sizes(2),
                         nBin(2),
                         pdf(2)
    {

        /* Default constructor */

        for (UInt iBin = 0; iBin < nBin + 1; iBin++)
        {
            bin_Edges[iBin] = DEFAULT_MIN_RADIUS * pow(1.2, iBin / 3.0);
        }

        for (UInt iBin = 0; iBin < nBin; iBin++)
        {
            bin_Centers[iBin] = 0.5 * (bin_Edges[iBin] + bin_Edges[iBin + 1]);
            bin_VCenters[iBin] = 4.0 / 3.0 * physConst::PI * (pow(bin_Edges[iBin], 3) + pow(bin_Edges[iBin + 1], 3)) * 0.5;
            bin_Sizes[iBin] = bin_Edges[iBin + 1] - bin_Edges[iBin];
        }

        nPart = 1.0;

        type = "lognormal";

        mu = 0;
        sigma = 0;
        alpha = 0;

        for (UInt iBin = 0; iBin < nBin; iBin++)
            pdf[iBin] = 0.0E+00;

    } /* End of Aerosol::Aerosol */

    Aerosol::Aerosol(Vector_1D bin_Centers_, Vector_1D bin_Edges_, RealDouble nPart_, RealDouble mu_, RealDouble sigma_, const char *distType, RealDouble alpha_, RealDouble gamma_, RealDouble b_) : bin_VCenters(bin_Centers_.size()),
                                                                                                                                                                                                      bin_Sizes(bin_Centers_.size()),
                                                                                                                                                                                                      type(distType),
                                                                                                                                                                                                      pdf(bin_Centers_.size())
    {

        /* Constructor */

        /* In the following cases, the aerosol pdf represents the following quantity: dn/d(ln(r))
         * The following identity can be used:
         * dn/dr = 1/r * dn/d(ln(r)) */

        /* Allocate bins centers, edges and number of bins */
        bin_Centers = bin_Centers_;
        bin_Edges = bin_Edges_;
        nBin = bin_Centers.size();

        if (nBin <= 0) { std::cout << "\nIn Aerosol::Aerosol: distribution has " << nBin << " bins!\n"; }

        if (nBin + 1 != bin_Edges.size()) { std::cout << "\nIn Aerosol::Aerosol: bin centers and/or edges are misshaped: " << 
                                            nBin + 1 << " != " << bin_Edges.size() << "\n"; }

        /* Compute size of each bin */
        for (UInt iBin = 0; iBin < nBin; iBin++)
        {
            bin_Sizes[iBin] = bin_Edges[iBin + 1] - bin_Edges[iBin];
            bin_VCenters[iBin] = 4.0 / 3.0 * physConst::PI * (pow(bin_Edges_[iBin], 3) + pow(bin_Edges_[iBin + 1], 3)) * 0.5;
        }

        /* Allocate number of particles */
        nPart = nPart_;

        /* Allocate mean and standard deviation */
        if (mu_ <= 0) { std::cout << "\nIn Aerosol::Aerosol: mean/mode is negative: mu = " << mu_ << "\n"; }
        mu = mu_;

        if (sigma_ <= 0) { std::cout << "\nIn Aerosol::Aerosol: standard deviation is negative: stddev = " << sigma_ << "\n"; }
        sigma = sigma_;

        /* Allocate alpha (only used in power law and generalized gamma distributions) */
        alpha = alpha_;

        /* Initialize distribution */

        if ((strcmp(type, "log") == 0) || (strcmp(type, "lognormal") == 0))
        {
            /* Log-normal distribution:
             * dn/d(ln(r)) = N / ( sqrt(2*\pi) * ln(sigma) ) * exp( - ( ln(r) - ln(r_m) ) ^ 2 / ( 2 ln(sigma) ^2 )) */

            if (sigma <= 1.0) { std::cout << "\nIn Aerosol::Aerosol: log-normal distribution requires that stddev > 1.0 ( stddev = " << sigma << " )\n"; }

            for (UInt iBin = 0; iBin < bin_Centers.size(); iBin++)
            {
                pdf[iBin] = nPart * exp(-0.5 * ((log(bin_Centers[iBin]) - log(mu)) / log(sigma)) * ((log(bin_Centers[iBin]) - log(mu)) / log(sigma))) / (sqrt(2.0 * physConst::PI) * log(sigma));
            }
        }
        else if ((strcmp(type, "norm") == 0) || (strcmp(type, "normal") == 0))
        {
            /* Normal distribution:
             * dn/d(ln(r)) = N / ( sqrt(2*\pi) * sigma ) * exp( - ( r - r_m ) ^ 2 / ( 2 * sigma ^2 )) */

            for (UInt iBin = 0; iBin < bin_Centers.size(); iBin++)
            {
                pdf[iBin] = nPart * exp(-0.5 * ((bin_Centers[iBin] - mu) / sigma) * ((bin_Centers[iBin] - mu) / sigma)) / (sqrt(2.0 * physConst::PI) * sigma);
            }
        }
        else if ((strcmp(type, "pow") == 0) || (strcmp(type, "power") == 0) || (strcmp(type, "junge") == 0))
        {
            /* Power distribution:
             * dn/d(ln(r)) = N * alpha * ( r / r_min ) ^ (-alpha) */

            if (alpha <= 0.0) { std::cout << "\nIn Aerosol::Aerosol: power law requires that alpha > 0 ( alpha = " << alpha_ << " )\n"; }

            for (UInt iBin = 0; iBin < bin_Centers.size(); iBin++)
            {
                pdf[iBin] = nPart * alpha * pow(bin_Centers[iBin] / bin_Centers[0], -alpha);
            }
        }
        else if ((strcmp(type, "gam") == 0) || (strcmp(type, "gamma") == 0) || (strcmp(type, "generalized gamma") == 0))
        {
            /* Gamma distribution:
             * dn/d(ln(r)) = gamma * b ^ ((alpha + 1)/gamma) / Gamma((alpha+1)/gamma) * r ^ (alpha + 1) * exp( - b * r ^ (gamma) )
             * alpha and gamma are of the order of unity (1 - 10)
             * b must be of the order of r ^ (-gamma) >> 1 */

            if ((alpha <= 0) || (std::fmod(alpha, 1.0) != 0)) { std::cout << "\nIn Aerosol::Aerosol: (generalized) gamma distribution requires that alpha is a positive integer ( alpha = " << alpha << " )\n"; }
            if (gamma_ <= 0) { std::cout << "\nIn Aerosol::Aerosol: (generalized) gamma distribution requires that gamma is positive ( gamma = " << gamma_ << " )\n"; }
            if (b_ <= 0) { std::cout << "\nIn Aerosol::Aerosol: (generalized) gamma distribution requires that b is positive ( b = " << b_ << " )\n"; }

            for (UInt iBin = 0; iBin < bin_Centers.size(); iBin++)
            {
                pdf[iBin] = nPart * gamma_ * pow(b_, (alpha + 1) / gamma_) / boost::math::tgamma((alpha + 1) / gamma_) * pow(bin_Centers[iBin], alpha + 1) * exp(-b_ * pow(bin_Centers[iBin], gamma_));
            }
        }
        else
        {
            std::cout << "\nIn Aerosol::Aerosol: distribution type must be either lognormal, normal, power or (generalized) gamma\n";
            std::cout << "\nCurrent type is " << type << "\n";
        }

        /* Check that we get the right number of particles */
        if ((std::abs(Moment() - nPart) / nPart > 0.10) && (nPart > 1.0E-10))
        {
            std::cout << "\nIn Aerosol::Aerosol: the size range doesn't cover the full distribution";
            std::cout << "\nPrescribed number: " << nPart << ", Number covered: " << Moment() << " [#/cm^3]\n";
            if ((strcmp(type, "log") == 0) || (strcmp(type, "lognormal") == 0))
            {
                std::cout << "For lognormal distribution, prescribed mode is: " << mu << " [m]\n";
            }
            else if ((strcmp(type, "norm") == 0) || (strcmp(type, "normal") == 0))
            {
                std::cout << "For normal distribution, prescribed mode is: " << mu << " [m]\n";
            }
        }

        /* For lognormal distributions:
         * rMode = rMean * exp( - ln(sigma)^2 );
         * rMedi = rMean * exp( 0.5 * ln(sigma)^2 ); */

    } /* End of Aerosol::Aerosol */

    Aerosol::Aerosol(const Aerosol &rhs)
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

    Aerosol::~Aerosol()
    {

        /* Destructor */

    } /* End of Aerosol::~Aerosol */

    Aerosol &Aerosol::operator=(const Aerosol &rhs)
    {

        if (&rhs == this)
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

    Aerosol Aerosol::operator+=(const Aerosol &rhs)
    {

        if (nBin != rhs.getNBin())
        {
            std::cout << "\nIn Aerosol::operator+=: aerosol distributions do not have the same number of bins: " << nBin << " != " << rhs.getNBin();
            return *this;
        }

        Vector_1D bin_Centers_rhs = rhs.getBinCenters();
        for (UInt iBin = 0; iBin < nBin; iBin++)
        {
            if (std::abs(bin_Centers[iBin] - bin_Centers_rhs[iBin]) > 1.0E-10)
            {
                std::cout << "\nIn Aerosol::operator+=: aerosol distributions do not have the same bin centers!";
                return *this;
            }
        }

        Vector_1D pdf_rhs = rhs.getPDF();
        for (UInt iBin = 0; iBin < nBin; iBin++)
        {
            pdf[iBin] += pdf_rhs[iBin];
        }

        return *this;

    } /* End of Aerosol::operator+= */

    Aerosol Aerosol::operator-=(const Aerosol &rhs)
    {

        if (nBin != rhs.getNBin())
        {
            std::cout << "\nIn Aerosol::operator+=: aerosol distributions do not have the same number of bins: " << nBin << " != " << rhs.getNBin();
            return *this;
        }

        Vector_1D bin_Centers_rhs = rhs.getBinCenters();
        for (UInt iBin = 0; iBin < nBin; iBin++)
        {
            if (std::abs(bin_Centers[iBin] - bin_Centers_rhs[iBin]) > 1.0E-10)
            {
                std::cout << "\nIn Aerosol::operator+=: aerosol distributions do not have the same bin centers!";
                return *this;
            }
        }

        Vector_1D pdf_rhs = rhs.getPDF();
        for (UInt iBin = 0; iBin < nBin; iBin++)
        {
            pdf[iBin] -= pdf_rhs[iBin];
        }

        return *this;

    } /* End of Aerosol::operator-= */

    Aerosol Aerosol::operator+(const Aerosol &rhs) const
    {

        Aerosol result = *this;

        if (nBin != rhs.getNBin())
        {
            std::cout << "\nIn Aerosol::operator+=: aerosol distributions do not have the same number of bins: " << nBin << " != " << rhs.getNBin();
            return *this;
        }

        Vector_1D bin_Centers_rhs = rhs.getBinCenters();
        for (UInt iBin = 0; iBin < nBin; iBin++)
        {
            if (std::abs(bin_Centers[iBin] - bin_Centers_rhs[iBin]) > 1.0E-10)
            {
                std::cout << "\nIn Aerosol::operator+=: aerosol distributions do not have the same bin centers!";
                return *this;
            }
        }

        Vector_1D pdf_rhs = rhs.getPDF();
        for (UInt iBin = 0; iBin < nBin; iBin++)
        {
            result.pdf[iBin] += pdf_rhs[iBin];
        }

        return result;

    } /* End of Aerosol::operator+ */

    Aerosol Aerosol::operator-(const Aerosol &rhs) const
    {

        Aerosol result = *this;

        if (nBin != rhs.getNBin())
        {
            std::cout << "\nIn Aerosol::operator+=: aerosol distributions do not have the same number of bins: " << nBin << " != " << rhs.getNBin();
            return *this;
        }

        Vector_1D bin_Centers_rhs = rhs.getBinCenters();
        for (UInt iBin = 0; iBin < nBin; iBin++)
        {
            if (std::abs(bin_Centers[iBin] - bin_Centers_rhs[iBin]) > 1.0E-10)
            {
                std::cout << "\nIn Aerosol::operator+=: aerosol distributions do not have the same bin centers!";
                return *this;
            }
        }

        Vector_1D pdf_rhs = rhs.getPDF();
        for (UInt iBin = 0; iBin < nBin; iBin++)
        {
            result.pdf[iBin] -= pdf_rhs[iBin];
        }

        return result;

    } /* End of Aerosol::operator- */

    void Aerosol::Coagulate(const RealDouble dt, const Coagulation &kernel)
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
        Vector_1D P(nBin, 0.0E+00);
        Vector_1D L(nBin, 0.0E+00);
        Vector_1D v(nBin);
        UInt iBin, jBin, kBin;

        for (iBin = 0; iBin < nBin; iBin++)
        {
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
        for (iBin = 0; iBin < nBin; iBin++)
        {

            /* Build production and loss terms */
            for (jBin = 0; jBin < nBin; jBin++)
            {
                if (jBin <= iBin)
                {

                    for (kBin = 0; kBin < iBin; kBin++)
                    {
                        /* k coagulating with j to form i */
                        if (kernel.f[kBin][jBin][iBin] != 0.0)
                            P[iBin] += kernel.f[kBin][jBin][iBin] * kernel.beta[kBin][jBin] * v_new[kBin] * pdf[jBin];
                    }
                }
                /* i coagulating with j to deplete i */
                if (kernel.f[iBin][jBin][iBin] != 1.0)
                    L[iBin] += (1.0 - kernel.f[iBin][jBin][iBin]) * kernel.beta[iBin][jBin] * pdf[jBin];
            }

            /* Non-mass conserving scheme: */
            //            v_new[iBin] = v[iBin] + ( P[iBin] - L[iBin] * v[iBin] ) * dt;

            /* Mass conserving scheme: */
            v_new[iBin] = (v[iBin] + dt * P[iBin]) / (1.0 + dt * L[iBin]);
        }

        /* For debug purposes */
        const bool checkMass = 0;
        if (checkMass)
        {
            std::cout << "Coagulation is a volume-conserving process. The two following quantities should be identical\n";
            std::cout << "At t     : " << Moment(3) << "\n";
        }

        for (iBin = 0; iBin < nBin; iBin++)
        {
            pdf[iBin] = 1.0e-6 * v_new[iBin] / bin_VCenters[iBin];
        }

        if (checkMass)
            std::cout << "At t + dt: " << Moment(3) << "\n";

    } /* End of Aerosol::Coagulate */

    void Aerosol::Coagulate(const RealDouble dt, const Vector_2D &beta, const Vector_3D &f)
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
        Vector_1D P(nBin);
        Vector_1D L(nBin);
        Vector_1D v(nBin);
        UInt iBin, jBin, kBin;

        for (iBin = 0; iBin < nBin; iBin++)
        {
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
        for (iBin = 0; iBin < nBin; iBin++)
        {
            for (jBin = 0; jBin < nBin; jBin++)
            {
                if (jBin <= iBin)
                {
                    for (kBin = 0; kBin < iBin; kBin++)
                    {
                        P[iBin] += f[kBin][jBin][iBin] * beta[kBin][jBin] * v_new[kBin] * pdf[jBin];
                    }
                }
                L[iBin] += (1.0 - f[iBin][jBin][iBin]) * beta[iBin][jBin] * pdf[jBin];
            }
            v_new[iBin] = (v[iBin] + dt * P[iBin]) / (1.0 + dt * L[iBin]);
            pdf[iBin] = 1.0e-6 * v_new[iBin] / bin_VCenters[iBin]; // Convert back to m3/cm3
        }

    } /* End of Aerosol::Coagulate */

    RealDouble Aerosol::Moment(UInt n) const
    {

        RealDouble moment = 0;

        for (UInt iBin = 0; iBin < nBin; iBin++)
        {
            moment += (log(bin_Edges[iBin + 1]) - log(bin_Edges[iBin])) * pow(bin_Centers[iBin], n) * pdf[iBin];
        }

        return moment;

    } /* End of Aerosol::Moment */

    RealDouble Aerosol::Radius() const
    {

        const RealDouble N = Moment(0);

        if (N > 0)
        {
            return Moment(1) / N;
        }
        else
        {
            std::cout << "\nIn Aerosol::Radius: Number of particles is " << N << " <= 0\n";
            return 0;
        }

    } /* End of Aerosol::Radius */

    RealDouble Aerosol::EffRadius() const
    {

        const RealDouble m2 = Moment(2);

        if (m2 > 0)
        {
            return Moment(3) / m2;
        }
        else
        {
            std::cout << "\nIn Aerosol::EffRadius: Second moment is " << m2 << " <= 0\n";
            return 0;
        }

    } /* End of Aerosol::EffRadius */

    RealDouble Aerosol::StdDev() const
    {

        const RealDouble N = Moment(0);

        if (N > 0)
        {
            return sqrt(Moment(2) / N - pow(Moment(1) / N, 2.0));
        }
        else
        {
            std::cout << "\nIn Aerosol::StdDev: Number of particles is " << N << " <= 0\n";
            return 0;
        }

    } /* End of Aerosol::StdDev */

    void Aerosol::scalePdf(RealDouble scalFactor)
    {

        if (scalFactor < 0)
        {
            std::cout << "\nIn Aerosol::scalePdf: scaling factor is negative ( " << scalFactor << " < 0 )\n";
        }

        for (UInt iBin = 0; iBin < nBin; iBin++)
        {
            pdf[iBin] *= scalFactor;
        }

    } /* End of Aerosol::scalePdf */

    void Aerosol::updatePdf(Vector_1D pdf_)
    {

        if (pdf_.size() != pdf.size())
        {
            std::cout << "\nIn Aerosol::updatePdf:: sizes differ ( " << pdf_.size() << " != " << pdf.size() << " )\n";
        }
        else
        {
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

    const char *Aerosol::getType() const
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

    Grid_Aerosol::Grid_Aerosol() : Nx(2),
                                   Ny(2),
                                   bin_Centers(2),
                                   bin_Edges(3),
                                   bin_VEdges(3),
                                   bin_Sizes(2),
                                   nBin(2)
    {

        /* Default constructor */

        for (UInt iBin = 0; iBin < nBin + 1; iBin++)
        {
            bin_Edges[iBin] = DEFAULT_MIN_RADIUS * pow(1.2, iBin / 3.0);
            bin_VEdges[iBin] = 4.0 / 3.0 * physConst::PI * bin_Edges[iBin] *
                               bin_Edges[iBin] * bin_Edges[iBin];
        }

        for (UInt iBin = 0; iBin < nBin; iBin++)
        {
            bin_Centers[iBin] = 0.5 * (bin_Edges[iBin] + bin_Edges[iBin + 1]);
            bin_Sizes[iBin] = bin_Edges[iBin + 1] - bin_Edges[iBin];
        }

        bin_VCenters.resize(nBin, Vector_2D(Ny, Vector_1D(Nx, 0.0)));

        RealDouble vol = 0.0E+00;
        for (UInt iBin = 0; iBin < nBin; iBin++)
        {
            vol = 0.5 * (bin_VEdges[iBin] + bin_VEdges[iBin + 1]);
            for (UInt jNy = 0; jNy < Ny; jNy++)
            {
                for (UInt iNx = 0; iNx < Nx; iNx++)
                {
                    bin_VCenters[iBin][jNy][iNx] = vol;
                }
            }
        }

        nPart = 1.0;

        type = "lognormal";

        mu = 0;
        sigma = 0;
        alpha = 0;

        pdf.resize(nBin, Vector_2D(Ny, Vector_1D(Nx, 0.0E+00)));

        for (UInt iBin = 0; iBin < nBin; iBin++)
        {
            for (UInt jNy = 0; jNy < Ny; jNy++)
            {
                for (UInt iNx = 0; iNx < Nx; iNx++)
                {
                    pdf[iBin][jNy][iNx] = 0.0E+00;
                }
            }
        }

    } /* End of Grid_Aerosol::Grid_Aerosol */

    Grid_Aerosol::Grid_Aerosol(UInt Nx_, UInt Ny_, Vector_1D bin_Centers_, Vector_1D bin_Edges_, RealDouble nPart_, RealDouble mu_, RealDouble sigma_, const char *distType, RealDouble alpha_, RealDouble gamma_, RealDouble b_) : Nx(Nx_),
                                                                                                                                                                                                                                    Ny(Ny_),
                                                                                                                                                                                                                                    bin_VEdges(bin_Centers_.size() + 1),
                                                                                                                                                                                                                                    bin_Sizes(bin_Centers_.size()),
                                                                                                                                                                                                                                    type(distType)
    {

        /* Constructor */

        /* In the following cases, the aerosol pdf represents the following quantity: dn/d(ln(r))
         * The following identity can be used:
         * dn/dr = 1/r * dn/d(ln(r)) */

        /* Allocate bins centers, edges and number of bins */
        bin_Centers = bin_Centers_;
        bin_Edges = bin_Edges_;
        nBin = bin_Centers.size();

        if (nBin <= 0) { std::cout << "\nIn Grid_Aerosol::Grid_Aerosol: distribution has " << nBin << " bins!\n"; }

        if (nBin + 1 != bin_Edges.size()) { std::cout << "\nIn Grid_Aerosol::Grid_Aerosol: bin centers and/or edges are misshaped: " 
                                    << nBin + 1 << " != " << bin_Edges.size() << "\n"; }

        /* Compute size of each bin */
        for (UInt iBin = 0; iBin < nBin; iBin++)
        {
            bin_VEdges[iBin] = 4.0 / 3.0 * physConst::PI * pow(bin_Edges[iBin], 3);
            bin_Sizes[iBin] = bin_Edges[iBin + 1] - bin_Edges[iBin];
        }
        bin_VEdges[nBin] = 4.0 / 3.0 * physConst::PI * pow(bin_Edges[nBin], 3);

        bin_VCenters.resize(nBin, Vector_2D(Ny, Vector_1D(Nx, 0.0E+00)));

        for (UInt iBin = 0; iBin < nBin; iBin++)
        {
            for (UInt jNy = 0; jNy < Ny; jNy++)
            {
                for (UInt iNx = 0; iNx < Nx; iNx++)
                    bin_VCenters[iBin][jNy][iNx] = 0.5 * (bin_VEdges[iBin] + bin_VEdges[iBin + 1]);
            }
        }

        pdf.resize(nBin, Vector_2D(Ny, Vector_1D(Nx, 0.0E+00)));

        /* Allocate number of particles */
        nPart = nPart_;

        /* Allocate mean and standard deviation */
        if (mu_ <= 0) { std::cout << "\nIn Grid_Aerosol::Grid_Aerosol: mean/mode is negative: mu = " << mu_ << "\n"; }
        mu = mu_;

        if (sigma_ <= 0) { std::cout << "\nIn Grid_Aerosol::Grid_Aerosol: standard deviation is negative: stddev = " << sigma_ << "\n"; }
        sigma = sigma_;

        /* Allocate alpha (only used in power law and generalized gamma distributions) */
        alpha = alpha_;

        /* Initialize distribution */
        if ((strcmp(type, "log") == 0) || (strcmp(type, "lognormal") == 0))
        {
            /* Log-normal distribution:
             * dn/d(ln(r)) = N / ( sqrt(2*\pi) * ln(sigma) ) * exp( - ( ln(r) - ln(r_m) ) ^ 2 / ( 2 ln(sigma) ^2 )) */

            if (sigma <= 1.0) {
                std::cout << "\nIn Grid_Aerosol::Grid_Aerosol: log-normal distribution requires that stddev > 1.0 ( stddev = " << sigma << " )\n"; }

            for (UInt iBin = 0; iBin < bin_Centers.size(); iBin++)
            {
                pdf[iBin][0][0] = nPart * exp(-0.5 * ((log(bin_Centers[iBin]) - log(mu)) / log(sigma)) * ((log(bin_Centers[iBin]) - log(mu)) / log(sigma))) / (sqrt(2.0 * physConst::PI) * log(sigma));
                for (UInt jNy = 0; jNy < Ny; jNy++)
                {
                    for (UInt iNx = 0; iNx < Nx; iNx++)
                    {
                        pdf[iBin][jNy][iNx] = pdf[iBin][0][0];
                    }
                }
            }
        }
        else if ((strcmp(type, "norm") == 0) || (strcmp(type, "normal") == 0))
        {
            /* Normal distribution:
             * dn/d(ln(r)) = N / ( sqrt(2*\pi) * sigma ) * exp( - ( r - r_m ) ^ 2 / ( 2 * sigma ^2 )) */

            for (UInt iBin = 0; iBin < bin_Centers.size(); iBin++)
            {
                pdf[iBin][0][0] = nPart * exp(-0.5 * ((bin_Centers[iBin] - mu) / sigma) * ((bin_Centers[iBin] - mu) / sigma)) / (sqrt(2.0 * physConst::PI) * sigma);
                for (UInt jNy = 0; jNy < Ny; jNy++)
                {
                    for (UInt iNx = 0; iNx < Nx; iNx++)
                    {
                        pdf[iBin][jNy][iNx] = pdf[iBin][0][0];
                    }
                }
            }
        }
        else if ((strcmp(type, "pow") == 0) || (strcmp(type, "power") == 0) || (strcmp(type, "junge") == 0))
        {
            /* Power distribution:
             * dn/d(ln(r)) = N * alpha * ( r / r_min ) ^ (-alpha) */

            if (alpha <= 0.0) {
                std::cout << "\nIn Grid_Aerosol::Grid_Aerosol: power law requires that alpha > 0 ( alpha = " << alpha_ << " )\n"; }

            for (UInt iBin = 0; iBin < bin_Centers.size(); iBin++)
            {
                pdf[iBin][0][0] = nPart * alpha * pow(bin_Centers[iBin] / bin_Centers[0], -alpha);
                for (UInt jNy = 0; jNy < Ny; jNy++)
                {
                    for (UInt iNx = 0; iNx < Nx; iNx++)
                    {
                        pdf[iBin][jNy][iNx] = pdf[iBin][0][0];
                    }
                }
            }
        }
        else if ((strcmp(type, "gam") == 0) || (strcmp(type, "gamma") == 0) || (strcmp(type, "generalized gamma") == 0))
        {
            /* Gamma distribution:
             * dn/d(ln(r)) = gamma * b ^ ((alpha + 1)/gamma) / Gamma((alpha+1)/gamma) * r ^ (alpha + 1) * exp( - b * r ^ (gamma) )
             * alpha and gamma are of the order of unity (1 - 10)
             * b must be of the order of r ^ (-gamma) >> 1 */

            if ((alpha <= 0) || (std::fmod(alpha, 1.0) != 0)) { std::cout << "\nIn Grid_Aerosol::Grid_Aerosol: (generalized) gamma distribution requires that alpha is a positive integer ( alpha = " << alpha << " )\n"; }
            if (gamma_ <= 0) { std::cout << "\nIn Grid_Aerosol::Grid_Aerosol: (generalized) gamma distribution requires that gamma is positive ( gamma = " << gamma_ << " )\n"; }
            if (b_ <= 0) { std::cout << "\nIn Grid_Aerosol::Grid_Aerosol: (generalized) gamma distribution requires that b is positive ( b = " << b_ << " )\n"; }

            for (UInt iBin = 0; iBin < bin_Centers.size(); iBin++)
            {
                pdf[iBin][0][0] = nPart * gamma_ * pow(b_, (alpha + 1) / gamma_) / boost::math::tgamma((alpha + 1) / gamma_) * pow(bin_Centers[iBin], alpha + 1) * exp(-b_ * pow(bin_Centers[iBin], gamma_));
                for (UInt jNy = 0; jNy < Ny; jNy++)
                {
                    for (UInt iNx = 0; iNx < Nx; iNx++)
                    {
                        pdf[iBin][jNy][iNx] = pdf[iBin][0][0];
                    }
                }
            }
        }
        else
        {
            std::cout << "\nIn Grid_Aerosol::Grid_Aerosol: distribution type must be either lognormal, normal, power or (generalized) gamma\n";
            std::cout << "\nCurrent type is " << type << "\n";
        }
    } /* End of Grid_Aerosol::Grid_Aerosol */

    Grid_Aerosol::Grid_Aerosol(const Grid_Aerosol &rhs)
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
        Nx = rhs.Nx;
        Ny = rhs.Ny;

    } /* End of Grid_Aerosol::Grid_Aerosol */

    Grid_Aerosol::~Grid_Aerosol()
    {
        /* Destructor */

    } /* End of Grid_Aerosol::~Grid_Aerosol */

    Grid_Aerosol &Grid_Aerosol::operator=(const Grid_Aerosol &rhs)
    {

        if (&rhs == this)
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

    Grid_Aerosol Grid_Aerosol::operator+=(const Grid_Aerosol &rhs)
    {   
        if (nBin != rhs.getNBin())
        {
            std::cout << "\nIn Aerosol::operator+=: aerosol distributions do not have the same number of bins: " << nBin << " != " << rhs.getNBin();
            return *this;
        }

        Vector_1D bin_Centers_rhs = rhs.getBinCenters();
        for (UInt iBin = 0; iBin < nBin; iBin++)
        {
            if (std::abs(bin_Centers[iBin] - bin_Centers_rhs[iBin]) > 1.0E-10)
            {
                std::cout << "\nIn Aerosol::operator+=: aerosol distributions do not have the same bin centers!";
                return *this;
            }
        }

        Vector_3D pdf_rhs = rhs.getPDF();
        for (UInt iBin = 0; iBin < nBin; iBin++)
        {
            for (UInt jNy = 0; jNy < Ny; jNy++)
            {
                for (UInt iNx = 0; iNx < Nx; iNx++)
                {
                    pdf[iBin][jNy][iNx] += pdf_rhs[iBin][jNy][iNx];
                }
            }
        }

        return *this;

    } /* End of Grid_Aerosol::operator+= */

    Grid_Aerosol Grid_Aerosol::operator-=(const Grid_Aerosol &rhs)
    {
        if (nBin != rhs.getNBin())
        {
            std::cout << "\nIn Grid_Aerosol::operator+=: aerosol distributions do not have the same number of bins: " << nBin << " != " << rhs.getNBin();
            return *this;
        }

        Vector_1D bin_Centers_rhs = rhs.getBinCenters();
        for (UInt iBin = 0; iBin < nBin; iBin++)
        {
            if (std::abs(bin_Centers[iBin] - bin_Centers_rhs[iBin]) > 1.0E-10)
            {
                std::cout << "\nIn Grid_Aerosol::operator+=: aerosol distributions do not have the same bin centers!";
                return *this;
            }
        }

        Vector_3D pdf_rhs = rhs.getPDF();
        for (UInt iBin = 0; iBin < nBin; iBin++)
        {
            for (UInt jNy = 0; jNy < Ny; jNy++)
            {
                for (UInt iNx = 0; iNx < Nx; iNx++)
                {
                    pdf[iBin][jNy][iNx] -= pdf_rhs[iBin][jNy][iNx];
                }
            }
        }

        return *this;

    } /* End of Grid_Aerosol::operator-= */

    Grid_Aerosol Grid_Aerosol::operator+(const Grid_Aerosol &rhs) const
    {
        Grid_Aerosol result = *this;

        if (nBin != rhs.getNBin())
        {
            std::cout << "\nIn Grid_Aerosol::operator+=: aerosol distributions do not have the same number of bins: " << nBin << " != " << rhs.getNBin();
            return *this;
        }

        Vector_1D bin_Centers_rhs = rhs.getBinCenters();
        for (UInt iBin = 0; iBin < nBin; iBin++)
        {
            if (std::abs(bin_Centers[iBin] - bin_Centers_rhs[iBin]) > 1.0E-10)
            {
                std::cout << "\nIn Grid_Aerosol::operator+=: aerosol distributions do not have the same bin centers!";
                return *this;
            }
        }

        Vector_3D pdf_rhs = rhs.getPDF();
        for (UInt iBin = 0; iBin < nBin; iBin++)
        {
            for (UInt jNy = 0; jNy < Ny; jNy++)
            {
                for (UInt iNx = 0; iNx < Nx; iNx++)
                {
                    result.pdf[iBin][jNy][iNx] += pdf_rhs[iBin][jNy][iNx];
                }
            }
        }

        return result;

    } /* End of Grid_Aerosol::operator+ */

    Grid_Aerosol Grid_Aerosol::operator-(const Grid_Aerosol &rhs) const
    {
        Grid_Aerosol result = *this;

        if (nBin != rhs.getNBin())
        {
            std::cout << "\nIn Grid_Aerosol::operator+=: aerosol distributions do not have the same number of bins: " << nBin << " != " << rhs.getNBin();
            return *this;
        }

        Vector_1D bin_Centers_rhs = rhs.getBinCenters();
        for (UInt iBin = 0; iBin < nBin; iBin++)
        {
            if (std::abs(bin_Centers[iBin] - bin_Centers_rhs[iBin]) > 1.0E-10)
            {
                std::cout << "\nIn Grid_Aerosol::operator+=: aerosol distributions do not have the same bin centers!";
                return *this;
            }
        }

        Vector_3D pdf_rhs = rhs.getPDF();
        for (UInt iBin = 0; iBin < nBin; iBin++)
        {
            for (UInt jNy = 0; jNy < Ny; jNy++)
            {
                for (UInt iNx = 0; iNx < Nx; iNx++)
                {
                    result.pdf[iBin][jNy][iNx] -= pdf_rhs[iBin][jNy][iNx];
                }
            }
        }

        return result;
    } /* End of Grid_Aerosol::operator- */

    void Grid_Aerosol::Coagulate(const RealDouble dt, Coagulation &kernel, const UInt N, const UInt SYM)
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
        if (checkMass)
        {
            std::cout << "Coagulation is a volume-conserving process. The two following quantities should be identical\n";
            std::cout << "At t     : " << Moment(3, Nx / 2, Ny / 2) * 1.0E+18 << "[um^3/cm^3]" << std::endl;
        }

        UInt Nx_max, Ny_max;

        bool performCoag = CheckCoagAndGrowInputs(N, SYM, Nx_max, Ny_max, "Coagulation");
        if(performCoag == false) { return; }

        /* Grid indices */
        UInt iNx = 0;
        UInt jNy = 0;

        /* Bin indices */
        UInt iBin = 0;
        UInt jBin = 0;
        UInt kBin = 0;
        UInt kBin_ = 0;

        /* Particle volume in each bin */
        Vector_3D v = Volume(); /* Expressed in [m^3/cm^3] */
        /* Copy v into v_new */
        Vector_3D v_new = v;

        /* Allocate variables */
        RealDouble P[nBin];
        RealDouble L[nBin];

        /* Total volume and number per grid cell */
        RealDouble totVol, nPart;

        /* Description of the algorithm:
         * \frac{dv}{dt}[iBin] = P - L * v[iBin]
         * Production     P = sum of all the bins (smaller than iBin) that
         *                    coagulate into a particle of size iBin
         * Loss L * v[iBin] = sum of all the bins that coagulate with iBin
         *
         * Scheme 1:
         * v_new - v = ( P - L * v ) * dt
         * Scheme 2:
         * v_new - v = ( P - L * v_new ) * dt
         * v_new = ( v + P * dt ) / ( 1.0 + L )
         * The latter is mass-conserving */

        for (jNy = 0; jNy < Ny_max; jNy++)
        {

            for (iNx = 0; iNx < Nx_max; iNx++)
            {

                /* Total aerosol volume */
                totVol = 0.0E+00;

                for (iBin = 0; iBin < nBin; iBin++)
                    totVol += v[iBin][jNy][iNx]; /* [m^3/cm^3] */

                if (totVol * 1E18 > 0.1)
                {
                    /* Only run coagulation where aerosol volume is greater
                     * than 0.1 um^3/cm^3 */

                    #pragma omp master
                    {
                        /* Update kernel (updates kernel.f and kernel.indices).
                         * This needs to be performed in serial as the class
                         * kernel is specific to the grid cell (jNy, iNx) */
                        kernel.buildF(bin_VCenters, jNy, iNx);
                    }

                    #pragma omp barrier

                    #pragma omp parallel for default(shared) private(iBin, jBin, kBin, kBin_, nPart) \
                        schedule(dynamic, 1) if (!PARALLEL_CASES)
                    for (iBin = 0; iBin < nBin; iBin++)
                    {

                        /* Reset P and L values */
                        P[iBin] = 0.0E+00;
                        L[iBin] = 0.0E+00;

                        /* Build production and loss terms */
                        for (jBin = 0; jBin < nBin; jBin++)
                        {

                            nPart = pdf[jBin][jNy][iNx] * log(bin_Edges[jBin + 1] / bin_Edges[jBin]);

                            if (jBin <= iBin)
                            {
                                for (kBin_ = 0; kBin_ < kernel.indices[jBin][iBin].size(); kBin_++)
                                {
                                    kBin = kernel.indices[jBin][iBin][kBin_];
                                    /* k coagulating with j to form i */
                                    if (kBin < iBin)
                                    {
                                        P[iBin] += kernel.f[kBin][jBin][iBin] * kernel.beta[kBin][jBin] * v_new[kBin][jNy][iNx] * nPart;
                                        /* [cm^3/#/s] * [m^3/cm^3] * [#/cm^3] = [m^3/cm^3/s] */
                                    }
                                }
                            }

                            /* i coagulating with j to deplete i */
                            if (kernel.f[iBin][jBin][iBin] != 1.0)
                                L[iBin] += (1.0 - kernel.f[iBin][jBin][iBin]) * kernel.beta[iBin][jBin] * nPart;
                        }

                        /* Mass conserving scheme: */
                        v_new[iBin][jNy][iNx] = (v[iBin][jNy][iNx] + dt * P[iBin]) / (1.0 + dt * L[iBin]);
                        if (v[iBin][jNy][iNx] > 0.0E+00)
                            pdf[iBin][jNy][iNx] *= v_new[iBin][jNy][iNx] / v[iBin][jNy][iNx];
                    }
                }
            }
        }

        /* Update bin centers */
        UpdateCenters(v_new, pdf);

        if (checkMass) { std::cout << "At t + dt: " << Moment(3, Nx / 2, Ny / 2) * 1.0E+18 << "[um^3/cm^3]" << std::endl; }
        
        //Apply Symmetry
        Vector_2D temp = Vector_2D();
        Vector_2D& temp1 = temp;
        CoagAndGrowApplySymmetry(N, SYM, Nx_max, Ny_max, "Coagulate", temp1);


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
        bool performGrowth = CheckCoagAndGrowInputs(N, SYM, Nx_max, Ny_max, "Grow");
        if(performGrowth == false) { return; }

        UInt iNx  = 0, jNy  = 0, iBin = 0;

        /* Conversion factor from ice volume [m^3] to [molecules] */ 
        const RealDouble UNITCONVERSION = physConst::RHO_ICE / MW_H2O * physConst::Na;

        /* Scaled Boltzmann constant */
        const RealDouble kB_ = physConst::kB * 1.00E+06;

        /* Declare and initialize particle totals and water vapor array */
        Vector_3D icePart = Number( );
        Vector_3D iceVol  = Volume( );
        Vector_2D totH2O  = H2O;
        
        RealDouble pSat;

        #pragma omp parallel if( !PARALLEL_CASES ) default( shared )
        {

            /* All declarations here are enforced as thread private */

            RealDouble partVol  = 0.0E+00;
            RealDouble icePart_ = 0.0E+00;
            RealDouble iceVol_  = 0.0E+00;

            int jBin = -1;
            std::vector<int> toBin( nBin, 0 );
            std::vector<int>::iterator iterBegin, iterCurr, iterEnd;

            /* Declare and initialize variable to store saturation quantities,
            * pressure and temperature */
            RealDouble locT = 0.0E+00;
            RealDouble locP = 0.0E+00;


            #pragma omp for                                                               \
            private ( iNx, jNy, iBin                                            ) \
            schedule( dynamic, 1                                                )
            for ( jNy = 0; jNy < Ny_max; jNy++ ) {
                for ( iNx = 0; iNx < Nx_max; iNx++ ) {
                    for ( iBin = 0; iBin < nBin; iBin++ ) {
                        totH2O[jNy][iNx] += iceVol[iBin][jNy][iNx] * UNITCONVERSION;
                        /* Unit check:
                        * [ molec/cm^3 ] = [ m^3 ice/cm^3 air ]   * [ molec/m^3 ice ] */
                    }
                }
            }

            #pragma omp for                                                               \
            private ( iNx, jNy, iBin, jBin, locP, locT              ) \
            private ( partVol, icePart_, iceVol_                                ) \
            schedule( dynamic, 1                                                )
            for ( jNy = 0; jNy < Ny_max; jNy++ ) {

                /* Store local pressure.
                * TODO: 
                * That might be moved into the loop over iNx eventually to 
                * account for 2D pressure met-fields?? */
                locP = P[jNy];

                for ( iNx = 0; iNx < Nx; iNx++ ) {
                    /* Store local temperature */
                    locT = T[jNy][iNx];

                    /* Store local saturation pressure w.r.t ice */
                    pSat = physFunc::pSat_H2Os( locT );

                    if ( H2O[jNy][iNx] * kB_ * locT / pSat > 0.0 ) {
                        APC_Scheme(jNy,iNx, locT, locP, dt, H2O, totH2O, icePart, iceVol );
                    }
                    /* ============== Moving-center structure ================ */
                    /* ======================================================= */
                    /* ============= Update bin center average =============== */


                    /* 1. Compute bin particle flux */
                    toBin = ComputeBinParticleFlux(iNx, jNy, iceVol, icePart);

                    /* 2. Attribute new particles according to fluxes */
                    ApplyBinParticleFlux(iNx, jNy, toBin, iceVol, icePart);
                }
            }
        } /* pragma omp parallel */

        //Apply Symmetries if there are any
        CoagAndGrowApplySymmetry(N, SYM, Nx_max, Ny_max, "Grow", H2O);
    } /* End of Grid::Aerosol::Grow */

    void Grid_Aerosol::APC_Scheme(const UInt jNy, const UInt iNx, const double T, const double P,
                            const double dt, Vector_2D& H2O, Vector_2D& totH2O, Vector_3D& icePart, Vector_3D& iceVol){
        
        RealDouble totPart = 0.0, totalkGrowth = 0.0, totalkGrowth_kelvin = 0.0, totH2Oi = 0.0;
        RealDouble pSat = physFunc::pSat_H2Os( T );
        Vector_1D  kGrowth(nBin, 0);
        RealDouble MAXVOL = bin_VEdges[nBin]; 
        RealDouble kB_ = physConst::kB * 1.00E+06; //SCALED boltzmann constant [J cm^3/K]
        RealDouble c_qit, C_qt, C_qsi; //Quantities used in APC scheme.
        RealDouble UNITCONVERSION = physConst::RHO_ICE / MW_H2O * physConst::Na; // Conversion factor from ice volume [m^3] to [molecules]
        /* -------------------------------------------------
        * Analytical predictor of condensation (APC) scheme
        * -------------------------------------------------
        * Mark Z. Jacobson, (1997), Numerical Techniques to
        * Solve Condensational and Dissolutional Growth
        * Equations When Growth is Coupled to Reversible
        * Reactions, Aerosol Science and Technology,
        * 27:4, 491-498, DOI: 10.1080/02786829708965489     */

        /* Check if partNum greater than a limit */
        for ( UInt iBin = 0; iBin < nBin; iBin++ ) {

            totPart += icePart[iBin][jNy][iNx];

        }

        /* Compute particle growth rates through ice deposition
        * We here assume that C_{s,i} is independent of the
        * bin and thus the particle size and only depends
        * on meteorological parameters. */
        if ( totPart < 0.00 ) { return; }
        for ( UInt iBin = 0; iBin < nBin; iBin++ ) {
        
            //Factor of 1e6 for cm3 - m3 conversion. 
            kGrowth[iBin] = 1.0e6 * icePart[iBin][jNy][iNx] * 4.0 * physConst::PI * bin_Centers[iBin]\
                * EffDiffCoef( bin_Centers[iBin], T, P, H2O[jNy][iNx]);  

            totalkGrowth += kGrowth[iBin];
            totalkGrowth_kelvin += kGrowth[iBin] * physFunc::Kelvin(bin_Centers[iBin]);
        }
        
        /* Compute the molar saturation concentration 
        * C_{q,s,i} in [mol/cm^3] */
        C_qsi = pSat / ( kB_ * T  * physConst::Na);
    

        /* Update gaseous molar concentration. S' is always 1 so not included.*/
        C_qt = ((H2O[jNy][iNx]/physConst::Na) + dt * totalkGrowth_kelvin * C_qsi) / (1.0 + dt*totalkGrowth);
        H2O[jNy][iNx] = C_qt * physConst::Na;

        
        /* Make sure that molecular water does not go over 
        * total water (gaseous + solid) concentrations */
        H2O[jNy][iNx] = std::min( H2O[jNy][iNx], totH2O[jNy][iNx] );
        
        for ( UInt iBin = 0; iBin < nBin; iBin++ ) {
            //Update molar concentration of ice [mol/cm3] and convert to volumetric concentration [m3/cm3]
            c_qit = (iceVol[iBin][jNy][iNx] * physConst::RHO_ICE / MW_H2O) + dt*kGrowth[iBin]*(C_qt - physFunc::Kelvin(bin_Centers[iBin])*C_qsi);
            iceVol[iBin][jNy][iNx] = c_qit * MW_H2O / physConst::RHO_ICE;
        
            iceVol[iBin][jNy][iNx] = \
                    std::min( std::max( iceVol[iBin][jNy][iNx], 0.0E+00 ), icePart[iBin][jNy][iNx] * MAXVOL );
        
            /* Compute total water taken up on particles */
            totH2Oi += iceVol[iBin][jNy][iNx] * UNITCONVERSION;
            /* Unit check:
            * [molec/cm^3 air] = [m^3 ice/cm^3 air] * [molec/m^3 ice] */
        }
        
        H2O[jNy][iNx] = totH2O[jNy][iNx] - totH2Oi; 
    } //End of Grid_Aerosol::APC_Scheme

    RealDouble Grid_Aerosol::EffDiffCoef( const RealDouble r, const RealDouble T, const RealDouble P, const RealDouble H2O ) const
    {

        /* DESCRIPTION:
         * Returns the "Effective Diffusion Coefficient" for ice crystal condensational growth [m^2/s]
         *
         * Source: Toon et al. 1989 "Physical Processes in Polar Stratospheric Clouds,
         *     and Jacobson, M.Z and Turco, R.P. "Simulating Condensational Growth, Evaporation, and Coagulation of Aerosols..." (1995)
         * Uses Jacobson's D_{eff} expression with a Kelvin effect correction in the numerator from Toon (1989) and other literature*/

        /* INPUT PARAMETERS:
         * - RealDouble r     :: radius in m
         * - RealDouble T     :: temperature in K
         * - RealDouble P     :: pressure in Pa
         * - RealDouble H2O   :: gaseous water concentration in molec/cm^3 air
         *
         * OUTPUT PARAMETERS:
         * - RealDouble :: Effective Diffusion Coefficient in m^2/s */

         const RealDouble dCoef = physFunc::CorrDiffCoef_H2O( r, T, P ); /* [m^2/s] */
         const RealDouble p_sat = physFunc::pSat_H2Os( T );
         const RealDouble p_h2o = physConst::R * T * H2O /(physConst::Na*1e-6);
         const RealDouble latS  = physFunc::LHeatSubl_H2O( T ); /* [J/kg] */
         const RealDouble nSat  = physFunc::pSat_H2Os( T ) / ( physConst::kB * T ); /* [#/m^3] */
        
        //Base D_{eff} expression in Jacobson (1995), equation 2
        //We can NOT directly use expressions for something like dm/dt, dr/dt, dv/dt, etc because of the (S' - 1) term in the 
        //numerator. Plugging one of those expressions into APC results in double-multiplication of (S' - 1), leading to a 1000%
        //for 110% ice supersaturation.
        //TODO Currently does not take into account Knudsen number effects for dCoef & ThermalCond or radiation effects.
         RealDouble Deff = dCoef/(1 + (dCoef * latS*latS*MW_H2O*MW_H2O * (H2O*1.0e6/physConst::Na)) / (physFunc::ThermalCond(r,T,P)*  physConst::R*T*T) );

         return Deff;
    } // End of Grid_Aerosol::EffDiffCoef 

    // TODO: Decide on a better way to handle ice particles that go above max volume. Currently,
    // they just stay in the highest volume box.
    std::vector<int> Grid_Aerosol::ComputeBinParticleFlux(const int x_index, const int y_index,
                                                          const Vector_3D &iceVol, const Vector_3D &icePart) const
    {
        std::vector<int> toBin(nBin, 0);
        double partVol;
        for (int iBin = 0; iBin < nBin; iBin++)
        {

            toBin[iBin] = -1;
            partVol = iceVol[iBin][y_index][x_index] / icePart[iBin][y_index][x_index];

            toBin[iBin] = std::lower_bound(bin_VEdges.begin(), bin_VEdges.end(), partVol) - bin_VEdges.begin() - 1;

            // Handling of ice crystals that grow past the max volume
            if (partVol > bin_VEdges[nBin])
            {
                // std::cout << "WARNING: ICE CRYSTALS GROWING PAST MAX BIN VOLUME" << std::endl;

                toBin[iBin] = nBin - 1;
            }
            // Ice crystals that drop below the minimum radius are considered lost.
            else if (toBin[iBin] == 0 && partVol < bin_VEdges[0])
            {
                toBin[iBin] = -1;
            }
        }
        return toBin;
    } //End of Grid_Aerosol::ComputeBinParticleFlux

    void Grid_Aerosol::ApplyBinParticleFlux(const int x_index, const int y_index,
                                            const std::vector<int> &toBin, const Vector_3D &iceVol, const Vector_3D &icePart)
    {

        std::vector<int>::const_iterator iterBegin, iterCurr, iterEnd;
        double icePart_, iceVol_;
        int jBin;
        for (int iBin = 0; iBin < nBin; iBin++)
        {

            icePart_ = 0.0E+00;
            iceVol_ = 0.0E+00;
            jBin = -1;
            iterBegin = toBin.begin();
            iterEnd = toBin.end();
            iterCurr = iterBegin;
            // TODO Can be optimized, use unordered_map to not have to loop through the entire toBin vector every time...
            // Sums up all ice particles being assigned to a given
            while ((iterCurr = std::find(iterCurr, iterEnd, iBin)) != iterEnd)
            {
                jBin = iterCurr - iterBegin; // j is less than or more than i, and this amount of will transfer to
                icePart_ += icePart[jBin][y_index][x_index];
                iceVol_ += iceVol[jBin][y_index][x_index];
                iterCurr++;
            }

            if (icePart_ > 0.0E+00)
            {
                // Bin is not empty. Compute particle volume, and clip it between min and max volume allowed.
                bin_VCenters[iBin][y_index][x_index] = std::max(std::min(iceVol_ / icePart_, bin_VEdges[iBin + 1]), bin_VEdges[iBin]);
                pdf[iBin][y_index][x_index] = icePart_ / (log(bin_Edges[iBin + 1] / bin_Edges[iBin]));
            }
            else
            {
                // Bin is empty. Set bin center to average volume of the bin.
                // This arbitrary value should not matter because there are no particles in this bin

                bin_VCenters[iBin][y_index][x_index] = 0.5 * (bin_VEdges[iBin] + bin_VEdges[iBin + 1]);
                pdf[iBin][y_index][x_index] = 0.0E+00;
            }
        }
    } //End of Grid_Aerosol::ApplyBinParticleFlux

    bool Grid_Aerosol::CheckCoagAndGrowInputs(const UInt N, const UInt SYM, UInt& Nx_max, UInt& Ny_max, const std::string funcName) const 
    {
        if (N == 0) { return false; } // Nothing is performed
        else if (N != 1 && N != 2){
            std::cout << " In Grid_Aerosol::" + funcName + ": Wrong input for N\n";
            std::cout << " N = " << N << "\n";
            return false;
        }
        else if (N == 1)
        {
            /* No emitted aerosols -> Aerosol is a uniform field. Perform coagulation only once */
            Nx_max = 1;
            Ny_max = 1;
        }
        else if (N == 2)
        {
            if(SYM != 0 && SYM != 1 && SYM != 2) {
                std::cout << " In Grid_Aerosol::" + funcName + ": Wrong input for SYM\n";
                std::cout << " SYM = " << SYM << "\n";
                return false;
            }
            //SYM: 2 = Both X and Y symmetry, 1 = Only sym around y axis, 0 = no symmetry
            Nx_max = (SYM == 2 || SYM == 1) ? Nx/2 : Nx;
            Ny_max = (SYM == 2) ? Ny/2 : Ny;

        }
        return true;
    }

    void Grid_Aerosol::CoagAndGrowApplySymmetry(const UInt N, const UInt SYM, const UInt Nx_max, const UInt Ny_max, const char* funcName, Vector_2D& H2O) 
    {
        if(strcmp(funcName, "Coagulate") != 0 && strcmp(funcName, "Grow") != 0){
            std::cout << "In Grid_Aerosol::CoagAndGrowApplySymmetry: " << funcName << "is not a valid function name input." <<std::endl;
            return;
        }

        UInt jNy, iNx, iBin;
        if ( N == 1 ) {

            /* Allocate uniform results to the grid */

            #pragma omp parallel for                                                      \
            default ( shared                                                ) \
            private ( iNx, jNy, iBin                                        ) \
            schedule( dynamic, 1                                            ) \
            if      ( !PARALLEL_CASES                                       )
            for ( jNy = 0; jNy < Ny; jNy++ ) {
                for ( iNx = 0; iNx < Nx; iNx++ ) {
                    if(strcmp(funcName, "Grow") == 0) { H2O[jNy][iNx] = H2O[0][0]; }
                    for ( iBin = 0; iBin < nBin; iBin++ ) {
                        pdf[iBin][jNy][iNx] = pdf[iBin][0][0];
                        bin_VCenters[iBin][jNy][iNx] = bin_VCenters[iBin][0][0];
                    }
                }
            }

        } else if ( N == 2 ) {

            /* Apply symmetry */

            if ( SYM == 2 ) {

                /* Symmetry around the origin */

                #pragma omp parallel for                                                      \
                default ( shared                                            ) \
                private ( iNx, jNy, iBin                                    ) \
                schedule( dynamic, 1                                        ) \
                if      ( !PARALLEL_CASES                                   )
                for ( jNy = 0; jNy < Ny; jNy++ ) {
                    for ( iNx = Nx_max; iNx < Nx; iNx++ ) {
                        if(strcmp(funcName, "Grow") == 0) { H2O[jNy][iNx] = H2O[jNy][Nx-1-iNx]; }
                        for ( iBin = 0; iBin < nBin; iBin++ ) {
                            pdf[iBin][jNy][iNx] = pdf[iBin][jNy][Nx-1-iNx];
                            bin_VCenters[iBin][jNy][iNx] = bin_VCenters[iBin][jNy][Nx-1-iNx];
                        }
                    }
                }

                #pragma omp parallel for                                                      \
                default ( shared                                            ) \
                private ( iNx, jNy, iBin                                    ) \
                schedule( dynamic, 1                                        ) \
                if      ( !PARALLEL_CASES                                   )
                for ( jNy = Ny_max; jNy < Ny; jNy++ ) {
                    for ( iNx = 0; iNx < Nx; iNx++ ) {
                        if(strcmp(funcName, "Grow") == 0) { H2O[jNy][iNx] = H2O[Ny-1-jNy][iNx]; };
                        for ( iBin = 0; iBin < nBin; iBin++ ) {
                            pdf[iBin][jNy][iNx] = pdf[iBin][Ny-1-jNy][iNx];
                            bin_VCenters[iBin][jNy][iNx] = bin_VCenters[iBin][Ny-1-jNy][iNx];
                        }
                    }
                }

                #pragma omp parallel for                                                      \
                default ( shared                                            ) \
                private ( iNx, jNy, iBin                                    ) \
                schedule( dynamic, 1                                        ) \
                if      ( !PARALLEL_CASES                                   )
                for ( jNy = Ny_max; jNy < Ny; jNy++ ) {
                    for ( iNx = Nx_max; iNx < Nx; iNx++ ) {
                        if(strcmp(funcName, "Grow") == 0) { H2O[jNy][iNx] = H2O[Ny-1-jNy][Nx-1-iNx]; };
                        for ( iBin = 0; iBin < nBin; iBin++ ) {
                            pdf[iBin][jNy][iNx] = pdf[iBin][Ny-1-jNy][Nx-1-iNx];
                            bin_VCenters[iBin][jNy][iNx] = bin_VCenters[iBin][Ny-1-jNy][Nx-1-iNx];
                        }
                    }
                }

            } else if ( SYM == 1 ) {

                /* Symmetry around the Y-axis */

                #pragma omp parallel for                                                      \
                default ( shared                                            ) \
                private ( iNx, jNy, iBin                                    ) \
                schedule( dynamic, 1                                        ) \
                if      ( !PARALLEL_CASES                                   )
                for ( jNy = 0; jNy < Ny; jNy++ ) {
                    for ( iNx = Nx_max; iNx < Nx; iNx++ ) {
                        if(strcmp(funcName, "Grow") == 0) { H2O[jNy][iNx] = H2O[jNy][Nx-1-iNx]; };
                        for ( iBin = 0; iBin < nBin; iBin++ ) {
                            pdf[iBin][jNy][iNx] = pdf[iBin][jNy][Nx-1-iNx];
                            bin_VCenters[iBin][jNy][iNx] = bin_VCenters[iBin][jNy][Nx-1-iNx];
                        }
                    }
                }
            }
        }
    }
    
    void Grid_Aerosol::UpdateCenters(const Vector_3D &iceV, const Vector_3D &PDF)
    {

        UInt iNx = 0;
        UInt jNy = 0;
        UInt iBin = 0;

        RealDouble ratio = 0.0E+00;

        #pragma omp parallel for default(shared) private(iNx, jNy, iBin, ratio) \
            schedule(dynamic, 1) if (!PARALLEL_CASES)
        for (UInt iBin = 0; iBin < nBin; iBin++)
        {
            ratio = log(bin_Edges[iBin + 1] / bin_Edges[iBin]);
            for (UInt jNy = 0; jNy < Ny; jNy++)
            {
                for (UInt iNx = 0; iNx < Nx; iNx++)
                {
                    if (PDF[iBin][jNy][iNx] > TINY)
                        bin_VCenters[iBin][jNy][iNx] =
                            std::max(std::min(iceV[iBin][jNy][iNx] / PDF[iBin][jNy][iNx] / ratio,
                                              0.9999 * bin_VEdges[iBin + 1]),
                                     1.0001 * bin_VEdges[iBin]);
                    else
                        bin_VCenters[iBin][jNy][iNx] =
                            0.5 * (bin_VEdges[iBin] + bin_VEdges[iBin + 1]);
                }
            }
        }

    } /* End of Grid_Aerosol::UpdateCenters */

    Vector_2D Grid_Aerosol::Moment(UInt n) const
    {

        UInt jNy = 0;
        UInt iNx = 0;
        UInt iBin = 0;

        Vector_2D moment(Ny, Vector_1D(Nx, 0.0E+00));
        const RealDouble FACTOR = 3.0 / RealDouble(4.0 * physConst::PI);

        #pragma omp parallel for default(shared) private(iNx, jNy, iBin) \
            schedule(dynamic, 1) if (!PARALLEL_CASES)
        for (jNy = 0; jNy < Ny; jNy++)
        {
            for (iNx = 0; iNx < Nx; iNx++)
            {
                for (iBin = 0; iBin < nBin; iBin++)
                {
                    moment[jNy][iNx] += (log(bin_Edges[iBin + 1] / bin_Edges[iBin])) * pow(FACTOR * bin_VCenters[iBin][jNy][iNx], n / RealDouble(3.0)) * pdf[iBin][jNy][iNx];
                }
            }
        }

        return moment;

    } /* End of Grid_Aerosol::Moment */

    Vector_3D Grid_Aerosol::Number() const
    {

        UInt jNy = 0;
        UInt iNx = 0;
        UInt iBin = 0;

        Vector_3D number(nBin, Vector_2D(Ny, Vector_1D(Nx, 0.0E+00)));
        RealDouble ratio = 0.0E+00;

        #pragma omp parallel for default(shared) private(iNx, jNy, iBin, ratio) \
            schedule(dynamic, 1) if (!PARALLEL_CASES)
        for (iBin = 0; iBin < nBin; iBin++)
        {
            ratio = log(bin_Edges[iBin + 1] / bin_Edges[iBin]);
            for (jNy = 0; jNy < Ny; jNy++)
            {
                for (iNx = 0; iNx < Nx; iNx++)
                {
                    number[iBin][jNy][iNx] = ratio * pdf[iBin][jNy][iNx];
                    /* Unit check: [#/cm^3] */
                }
            }
        }

        return number;

    } /* End of Grid_Aerosol::Number */

    Vector_2D Grid_Aerosol::TotalNumber() const
    {

        return Moment(0);

    } /* End of Grid_Aerosol::TotalNumber */

    RealDouble Grid_Aerosol::TotalNumber_sum(const Vector_2D cellAreas) const
    {

        Vector_2D TotalNumber_pcell = TotalNumber();
        RealDouble totalnumber_sum = 0;
        UInt iNx = 0;
        UInt jNy = 0;

        #pragma omp parallel for default(shared) private(iNx, jNy) \
            reduction(+                                            \
                    : totalnumber_sum)                           \
                schedule(dynamic, 1) if (!PARALLEL_CASES)
        for (jNy = 0; jNy < Ny; jNy++)
        {
            for (iNx = 0; iNx < Nx; iNx++)
            {
                totalnumber_sum += TotalNumber_pcell[jNy][iNx] * cellAreas[jNy][iNx] * 1.0E+06;
                                    // part/cm3 * m2 * cm3/m3 = part/m
            }
        }

        return totalnumber_sum;
    }

    Vector_1D Grid_Aerosol::Overall_Size_Dist(const Vector_2D cellAreas) const
    {

        Vector_1D overall_size_dist(nBin, 0.00E+00);
        UInt jNy = 0;
        UInt iNx = 0;
        UInt iBin = 0;

        #pragma omp parallel for default(shared) private(iNx, jNy, iBin) \
            schedule(dynamic, 1) if (!PARALLEL_CASES)
        for (iBin = 0; iBin < nBin; iBin++)
        {
            for (jNy = 0; jNy < Ny; jNy++)
            {
                for (iNx = 0; iNx < Nx; iNx++)
                {
                    overall_size_dist[iBin] += pdf[iBin][jNy][iNx] * cellAreas[jNy][iNx] * 1.0E+06;
                }
            }
        }

        return overall_size_dist;
    }

    Vector_3D Grid_Aerosol::Volume() const
    {

        UInt jNy = 0;
        UInt iNx = 0;
        UInt iBin = 0;

        Vector_3D volume(nBin, Vector_2D(Ny, Vector_1D(Nx, 0.0E+00)));
        RealDouble ratio = 0.0E+00;

        #pragma omp parallel for default(shared) private(iNx, jNy, iBin, ratio) \
            schedule(dynamic, 1) if (!PARALLEL_CASES)
        for (iBin = 0; iBin < nBin; iBin++)
        {
            ratio = log(bin_Edges[iBin + 1] / bin_Edges[iBin]);
            for (jNy = 0; jNy < Ny; jNy++)
            {
                for (iNx = 0; iNx < Nx; iNx++)
                {
                    volume[iBin][jNy][iNx] = ratio * bin_VCenters[iBin][jNy][iNx] *
                                             pdf[iBin][jNy][iNx];
                    /* Unit check:                   [m^3] * \
                     *                               [#/cm^3] \
                     *                             = [m^3/cm^3] */
                }
            }
        }

        return volume;

    } /* End of Grid_Aerosol::Volume */

    Vector_2D Grid_Aerosol::TotalArea() const
    {

        UInt jNy = 0;
        UInt iNx = 0;

        Vector_2D m2 = Moment(2);
        const RealDouble FACTOR = 4.0 * physConst::PI;
        /* A = 4.0*pi*m2 */

        #pragma omp parallel for default(shared) private(iNx, jNy) \
            schedule(dynamic, 1) if (!PARALLEL_CASES)
        for (jNy = 0; jNy < Ny; jNy++)
        {
            for (iNx = 0; iNx < Nx; iNx++)
                m2[jNy][iNx] = m2[jNy][iNx] * FACTOR;
        }

        return m2;

    } /* End of Grid_Aerosol::TotalVolume */

    Vector_2D Grid_Aerosol::TotalVolume() const
    {

        UInt jNy = 0;
        UInt iNx = 0;

        Vector_2D m3 = Moment(3);
        const RealDouble FACTOR = 4.0 / RealDouble(3.0) * physConst::PI;
        /* V = 4.0/3.0*pi*m3 */

        #pragma omp parallel for default(shared) private(iNx, jNy) \
            schedule(dynamic, 1) if (!PARALLEL_CASES)
        for (jNy = 0; jNy < Ny; jNy++)
        {
            for (iNx = 0; iNx < Nx; iNx++)
                m3[jNy][iNx] = m3[jNy][iNx] * FACTOR;
        }

        return m3;

    } /* End of Grid_Aerosol::TotalVolume */

    RealDouble Grid_Aerosol::TotalIceMass_sum(const Vector_2D cellAreas) const
    {
        Vector_2D iwc = IWC();
        RealDouble totalicemass_sum = 0;
        UInt iNx = 0;
        UInt jNy = 0;

        #pragma omp parallel for default(shared) private(iNx, jNy) \
            reduction(+                                            \
                    : totalicemass_sum)                           \
                schedule(dynamic, 1) if (!PARALLEL_CASES)
        for (jNy = 0; jNy < Ny; jNy++)
        {
            for (iNx = 0; iNx < Nx; iNx++)
            {
                totalicemass_sum += iwc[jNy][iNx] * cellAreas[jNy][iNx];
            }
        }
        return totalicemass_sum;
    }

    Vector_2D Grid_Aerosol::IWC() const
    {

        UInt jNy = 0;
        UInt iNx = 0;

        Vector_2D TVol = TotalVolume();
        const RealDouble FACTOR = physConst::RHO_ICE * 1.0E+06;

        #pragma omp parallel for default(shared) private(iNx, jNy) \
            schedule(dynamic, 1) if (!PARALLEL_CASES)
        for (jNy = 0; jNy < Ny; jNy++)
        {
            for (iNx = 0; iNx < Nx; iNx++)
                TVol[jNy][iNx] *= FACTOR;
            /* Unit check: [m^3/cm^3] * [kg/m^3] * [cm^3/m^3] = [kg/m^3] */
        }

        return TVol;

    } /* End of Grid_Aerosol::IWC */

    Vector_2D Grid_Aerosol::Extinction() const
    {

        UInt jNy = 0;
        UInt iNx = 0;

        Vector_2D chi = IWC();
        Vector_2D rE = EffRadius();

        const RealDouble a = 3.448E+00; /* [m^2/kg] */
        const RealDouble b = 2.431E-03; /* [m^3/kg] */

        #pragma omp parallel for default(shared) private(iNx, jNy) \
            schedule(dynamic, 1) if (!PARALLEL_CASES)
        for (jNy = 0; jNy < Ny; jNy++)
        {
            for (iNx = 0; iNx < Nx; iNx++)
            {
                if (rE[jNy][iNx] > 1.00E-15)
                {
                    chi[jNy][iNx] = chi[jNy][iNx] * (a + b / rE[jNy][iNx]);
                    /* Unit check:
                     * [1/m] = [kg/m^3] * [m^2/kg] */
                }
                else
                    chi[jNy][iNx] = 0.0E+00;
            }
        }

        return chi;

    } /* End of Grid_Aerosol::Extinction */

    Vector_1D Grid_Aerosol::PDF_Total(const Vector_2D &cellAreas) const
    {

        UInt jNy = 0;
        UInt iNx = 0;
        UInt iBin = 0;

        Vector_1D PDF(nBin, 0.0E+00);

        #pragma omp parallel for default(shared) private(iNx, jNy, iBin) \
            schedule(dynamic, 1) if (!PARALLEL_CASES)
        for (iBin = 0; iBin < nBin; iBin++)
        {
            for (jNy = 0; jNy < Ny; jNy++)
            {
                for (iNx = 0; iNx < Nx; iNx++)
                {
                    PDF[iBin] = PDF[iBin] + pdf[iBin][jNy][iNx] * cellAreas[jNy][iNx] * 1.0E+06;
                    /* Unit check:
                     * [#/cm^3] * [m^2] * [cm^3/m^3] = [#/m] */
                }
            }
        }

        return PDF;

    } /* End of Grid_Aerosol::PDF_Total */

    Vector_1D Grid_Aerosol::PDF_Total(const Mesh &m) const
    {

        const Vector_2D cellAreas = m.areas();

        UInt jNy = 0;
        UInt iNx = 0;
        UInt iBin = 0;

        Vector_1D PDF(nBin, 0.0E+00);

        #pragma omp parallel for default(shared) private(iNx, jNy, iBin) \
            schedule(dynamic, 1) if (!PARALLEL_CASES)
        for (iBin = 0; iBin < nBin; iBin++)
        {
            for (jNy = 0; jNy < Ny; jNy++)
            {
                for (iNx = 0; iNx < Nx; iNx++)
                {
                    PDF[iBin] = PDF[iBin] + pdf[iBin][jNy][iNx] * cellAreas[jNy][iNx] * 1.0E+06;
                    /* Unit check:
                     * [#/cm^3] * [m^2] * [cm^3/m^3] = [#/m] */
                }
            }
        }

        return PDF;

    } /* End of Grid_Aerosol::PDF_Total */

    Vector_1D Grid_Aerosol::xOD(const Vector_1D dx) const
    {

        UInt jNy = 0;
        UInt iNx = 0;

        Vector_1D tau_x(Ny, 0.0E+00);
        Vector_2D chi = Extinction();

        #pragma omp parallel for default(shared) private(iNx, jNy) \
            schedule(dynamic, 1) if (!PARALLEL_CASES)
        for (jNy = 0; jNy < Ny; jNy++)
        {
            for (iNx = 0; iNx < Nx; iNx++)
                tau_x[jNy] += dx[iNx] * chi[jNy][iNx];
        }

        return tau_x;

    } /* End of Grid_Aerosol::tau_x */

    Vector_1D Grid_Aerosol::yOD(const Vector_1D dy) const
    {

        UInt jNy = 0;
        UInt iNx = 0;

        Vector_1D tau_y(Nx, 0.0E+00);
        Vector_2D chi = Extinction();

        #pragma omp parallel for default(shared) private(iNx, jNy) \
            schedule(dynamic, 1) if (!PARALLEL_CASES)
        for (iNx = 0; iNx < Nx; iNx++)
        {
            for (jNy = 0; jNy < Ny; jNy++)
                tau_y[iNx] += dy[jNy] * chi[jNy][iNx];
        }

        return tau_y;

    } /* End of Grid_Aerosol::tau_y */

    Vector_2D Grid_Aerosol::Radius() const
    {

        UInt jNy = 0;
        UInt iNx = 0;

        Vector_2D r(Ny, Vector_1D(Nx, 0.0E+00));

        const Vector_2D m0 = Moment(0);
        const Vector_2D m1 = Moment(1);

        #pragma omp parallel for default(shared) private(iNx, jNy) \
            schedule(dynamic, 1) if (!PARALLEL_CASES)
        for (jNy = 0; jNy < Ny; jNy++)
        {
            for (iNx = 0; iNx < Nx; iNx++)
            {
                if (m0[jNy][iNx] > TINY)
                    r[jNy][iNx] = m1[jNy][iNx] / m0[jNy][iNx];
                else
                    r[jNy][iNx] = 0.0E+00;
            }
        }

        return r;

    } /* End of Grid_Aerosol::Radius */

    Vector_2D Grid_Aerosol::EffRadius() const
    {

        UInt jNy = 0;
        UInt iNx = 0;

        Vector_2D r_eff(Ny, Vector_1D(Nx, 0.0E+00));

        const Vector_2D m2 = Moment(2);
        const Vector_2D m3 = Moment(3);

        #pragma omp parallel for default(shared) private(iNx, jNy) \
            schedule(dynamic, 1) if (!PARALLEL_CASES)
        for (jNy = 0; jNy < Ny; jNy++)
        {
            for (iNx = 0; iNx < Nx; iNx++)
            {
                if (m2[jNy][iNx] > TINY)
                    r_eff[jNy][iNx] = m3[jNy][iNx] / m2[jNy][iNx];
                else
                    r_eff[jNy][iNx] = 0.0E+00;
            }
        }

        return r_eff;

    } /* End of Grid_Aerosol::EffRadius */

    Vector_2D Grid_Aerosol::StdDev() const
    {

        UInt jNy = 0;
        UInt iNx = 0;

        Vector_2D sigma(Ny, Vector_1D(Nx, 0.0E+00));

        const Vector_2D m0 = Moment(0);
        const Vector_2D m1 = Moment(1);
        const Vector_2D m2 = Moment(2);

        #pragma omp parallel for default(shared) private(iNx, jNy) \
            schedule(dynamic, 1) if (!PARALLEL_CASES)
        for (jNy = 0; jNy < Ny; jNy++)
        {
            for (iNx = 0; iNx < Nx; iNx++)
            {
                if (m0[jNy][iNx] > TINY)
                    sigma[jNy][iNx] = sqrt(m2[jNy][iNx] / m0[jNy][iNx] - m1[jNy][iNx] * m1[jNy][iNx] / (m0[jNy][iNx] * m0[jNy][iNx]));
                else
                    sigma[jNy][iNx] = 0.0E+00;
            }
        }

        return sigma;

    } /* End of Grid_Aerosol::StdDev */

    RealDouble Grid_Aerosol::Moment(UInt n, Vector_1D PDF) const
    {

        UInt iBin = 0;

        RealDouble moment = 0.0E+00;

        #pragma omp parallel for default(shared) private(iBin) \
            reduction(+                                        \
                    : moment)                                \
                schedule(dynamic, 1) if (!PARALLEL_CASES)
        for (iBin = 0; iBin < nBin; iBin++)
        {
            moment += (log(bin_Edges[iBin + 1] / bin_Edges[iBin])) * pow(bin_Centers[iBin], n) * PDF[iBin];
        }

        return moment;

    } /* End of Grid_Aerosol::Moment */

    RealDouble Grid_Aerosol::Moment(UInt n, UInt jNy, UInt iNx) const
    {

        UInt iBin = 0;

        RealDouble moment = 0.0E+00;
        const RealDouble FACTOR = 3.0 / RealDouble(4.0 * physConst::PI);

        #pragma omp parallel for default(shared) private(iBin) \
            reduction(+                                        \
                    : moment)                                \
                schedule(dynamic, 1) if (!PARALLEL_CASES)
        for (iBin = 0; iBin < nBin; iBin++)
            moment += (log(bin_Edges[iBin + 1] / bin_Edges[iBin])) * pow(FACTOR * bin_VCenters[iBin][jNy][iNx], n / RealDouble(3.0)) * pdf[iBin][jNy][iNx];

        return moment;

    } /* End of Grid_Aerosol::Moment */

    RealDouble Grid_Aerosol::Radius(UInt jNy, UInt iNx) const
    {

        const RealDouble N = Moment(0, jNy, iNx);

        if (N > 0)
        {
            return Moment(1, jNy, iNx) / N;
        }
        else
        {
            // std::cout << "\nIn Grid_Aerosol::Radius: Number of particles is " << N << " <= 0\n";
            return 0;
        }

    } /* End of Grid_Aerosol::Radius */

    RealDouble Grid_Aerosol::EffRadius(UInt jNy, UInt iNx) const
    {

        const RealDouble m2 = Moment(2, jNy, iNx);

        if (m2 > 0)
        {
            return Moment(3, jNy, iNx) / m2;
        }
        else
        {
            // std::cout << "\nIn Grid_Aerosol::EffRadius: Second moment is " << m2 << " <= 0\n";
            return 0;
        }

    } /* End of Aerosol::EffRadius */

    RealDouble Grid_Aerosol::StdDev(UInt jNy, UInt iNx) const
    {

        const RealDouble N = Moment(0, jNy, iNx);

        if (N > 0)
        {
            return sqrt(Moment(2, jNy, iNx) / N - pow(Moment(1, jNy, iNx) / N, 2.0));
        }
        else
        {
            // std::cout << "\nIn Grid_Aerosol::StdDev: Number of particles is " << N << " <= 0\n";
            return 0;
        }

    } /* End of Grid_Aerosol::StdDev */

    void Grid_Aerosol::updatePdf(Vector_3D pdf_)
    {

        if (pdf_.size() != pdf.size())
        {
            std::cout << "\nIn Grid_Aerosol::updatePdf:: sizes differ ( " << pdf_.size() << " != " << pdf.size() << " )\n";
        }
        else
        {
            pdf = pdf_;
        }

    } /* End of Grid_Aerosol::updatePdf */

    Vector_1D Grid_Aerosol::Average(const Vector_2D &weights,
                                    const RealDouble &totalWeight) const
    {

        UInt iNx = 0;
        UInt jNy = 0;
        UInt iBin = 0;

        Vector_1D out(4, 0.0E+00);

        Vector_1D PDF(nBin, 0.0E+00);

        RealDouble w = 0.0E+00;

        #pragma omp parallel for default(shared) private(iNx, jNy, iBin) \
            schedule(dynamic, 1) if (!PARALLEL_CASES)
        for (iBin = 0; iBin < nBin; iBin++)
        {
            for (jNy = 0; jNy < Ny; jNy++)
            {
                for (iNx = 0; iNx < Nx; iNx++)
                    PDF[iBin] += pdf[iBin][jNy][iNx] * weights[jNy][iNx] / totalWeight;
            }
        }

        out[0] = Moment(0, PDF);
        if (out[0] > TINY)
        {
            out[2] = 4.0 * physConst::PI * Moment(2, PDF);
            out[3] = 4.0 / 3.0 * physConst::PI * Moment(3, PDF);
            out[1] = 3.0 * out[3] / out[2];
        }
        else
        {
            out[1] = 1.00E-07;
            out[2] = 0.00E+00;
            out[3] = 0.00E+00;
        }

        return out;

    } /* End of Grid_Aerosol::Average */

    void Grid_Aerosol::addPDF(const Aerosol &aerosol, const Vector_2D &weights,
                              const Vector_2D &cellAreas, const RealDouble nCell)
    {

        UInt iNx = 0;
        UInt jNy = 0;
        UInt iBin = 0;

        Vector_1D AerPDF = aerosol.getPDF();

        #pragma omp parallel for default(shared) private(iNx, jNy, iBin) \
            schedule(dynamic, 1) if (!PARALLEL_CASES)
        for (jNy = 0; jNy < Ny; jNy++)
        {
            for (iNx = 0; iNx < Nx; iNx++)
            {
                if (weights[jNy][iNx] != 0.0E+00)
                {
                    for (iBin = 0; iBin < nBin; iBin++)
                        pdf[iBin][jNy][iNx] += AerPDF[iBin] / (nCell * cellAreas[jNy][iNx]);
                }
            }
        }

    } /* End of Grid_Aerosol::addPDF */

    void Grid_Aerosol::addPDF(const Vector_1D &PDF, const Vector_2D &weights,
                              const Vector_2D &cellAreas, const RealDouble nCell)
    {

        UInt iNx = 0;
        UInt jNy = 0;
        UInt iBin = 0;

        #pragma omp parallel for default(shared) private(iNx, jNy, iBin) \
            schedule(dynamic, 1) if (!PARALLEL_CASES)
        for (jNy = 0; jNy < Ny; jNy++)
        {
            for (iNx = 0; iNx < Nx; iNx++)
            {
                if (weights[jNy][iNx] != 0.0E+00)
                {
                    for (iBin = 0; iBin < nBin; iBin++)
                        pdf[iBin][jNy][iNx] += PDF[iBin] / (nCell * cellAreas[jNy][iNx]);
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

Vector_1D AIM::Aerosol::generateEndOfEPMTestPDF(const Vector_1D& bin_Centers){
    const int nBin = bin_Centers.size();
    Vector_1D pdf (bin_Centers.size(), 0);
    const double R_MIN = bin_Centers[0];
    const double R_MAX = bin_Centers[nBin-1];
    double slope, radius;
    for (int i = 0; i < nBin; i++){
        radius = bin_Centers[i];
        if(radius < 0.500e-6) {
            pdf[i] = std::max(0.00, 378.65543+22.652366*log(radius));
        }
        else if (radius >= 0.500e-6 && radius <= 2.00e-6){
            pdf[i] = 465.22382 + 28.67319*log(radius); 
        }
        else if (radius > 2.00e-6 && radius <= 5.0e-6){
            pdf[i] = -979.1635 - 81.28436*log(radius);
        }
        else if (radius > 5.0e-6 && radius <= 9.0e-6){
            pdf[i] = -275.7263 - 23.81817*log(radius);
        }
        else{
            pdf[i] = std::max(0.0, -21.74413 - 1.957615 * log(radius));
        }
    }
    return pdf;
}

/* End of Aerosol.cpp */
