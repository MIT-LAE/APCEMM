#include "AIM/Aerosol.hpp"
#include "Util/ForwardDecl.hpp"
#include "Util/PhysConstant.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <fstream>
#include <iostream>

using namespace AIM;

double volume_fraction(double V, double v1, double v2) {
    return ((v2-V)/(v2-v1)) * v1/V;
}

TEST_CASE ("Aerosol", "[single-file]" ) {

    
    int nBins = 100;
    double r_min = 1e-8;
    double r_max = 10e-6;
    // double bin_size = (r_max - r_min)/nBins;

    Vector_1D bin_centers(nBins);
    Vector_1D bin_edges(nBins+1);
    Vector_1D volume_bin_centers(nBins); 

    double ratio = r_max/r_min;
    double exponent;


    for (int i = 0; i < nBins; i++) {
        exponent = double(i) / double(nBins) ; 
        bin_centers[i] = 0.5 * r_min * pow(ratio, exponent) *( 1+ pow(ratio, 1.0/nBins));
        bin_edges[i] = r_min * pow(ratio, exponent);
        volume_bin_centers[i] = pow(bin_centers[i], 3.0) * (4.0*physConst::PI)/3.0 ;
    }
    bin_edges[nBins] = r_max;

    double nPart = 1.0e12;
    double mu = 0.4e-7;
    double sigma = 1.5;

    Aerosol test_pdf(bin_centers, bin_edges, nPart, mu, sigma);

    Vector_1D pdf = test_pdf.getPDF();

    SECTION("getBinCenters") {
        Vector_1D retrieved_centers = test_pdf.getBinCenters();
        for (int i = 0; i < nBins; i++) {
            REQUIRE(retrieved_centers[i] == Catch::Approx(bin_centers[i]));
        }
    }

    SECTION("getBinVCenters") {
        Vector_1D retrieved_volumes = test_pdf.getBinVCenters();
        for (int i = 0; i < nBins; i++) {
            REQUIRE(retrieved_volumes[i] == Catch::Approx(volume_bin_centers[i]).margin(1e-12));
        }
    }

    SECTION("getBinEdges") {
        Vector_1D retrieved_edges = test_pdf.getBinEdges();
        for (int i = 0; i < nBins+1; i++) {
            REQUIRE(retrieved_edges[i] == Catch::Approx(bin_edges[i]));
        }
    }

    SECTION("getNBin") {
        int retrieved_Nbin = test_pdf.getNBin();
        REQUIRE(retrieved_Nbin == nBins);
    }

    SECTION("getType") {
        const char* retrieved_type = test_pdf.getType();
        REQUIRE(strcmp(retrieved_type, "lognormal") == 0);
        
    }

    SECTION("Moment 0") {
        double moment0 = test_pdf.Moment(0);
        REQUIRE(moment0 == Catch::Approx(nPart).epsilon(0.1));
    }

    SECTION("Moment 1") {
        // See equation 8.44 in Seinfeld and Pandis (2006)
        double moment1 = test_pdf.Moment(1);
        REQUIRE(moment1 == Catch::Approx(nPart * mu * exp(pow(log(sigma), 2)/2.0)).epsilon(0.1));
    }

    SECTION("Moment 2") {

        double moment2 = test_pdf.Moment(2);
        REQUIRE(moment2 == Catch::Approx(nPart * mu*mu * exp(2*pow(log(sigma),2) )) );
    }

    SECTION("Moment 3") {

        double moment3 = test_pdf.Moment(3);
        REQUIRE(moment3 == Catch::Approx(nPart * mu*mu *mu* exp((9.0/2.0)*pow(log(sigma),2))) );
    }
    SECTION("Moment 4") {

        double moment4 = test_pdf.Moment(4);
        REQUIRE(moment4 == Catch::Approx(nPart * pow(mu,4) * exp(8.0*pow(log(sigma),2))) );
    }

    SECTION("Radius") {
        // See equation 8.44 in Seinfeld and Pandis (2006)
        double radius = test_pdf.Radius();
        REQUIRE(radius == Catch::Approx(mu * exp(pow(log(sigma), 2)/2.0)).epsilon(0.1));
    }

    SECTION("Effective radius") {

        double m3 = test_pdf.Moment(3);
        double m2 = test_pdf.Moment(2);
        REQUIRE(test_pdf.EffRadius() == Catch::Approx(m3/m2));

    }

    SECTION("Standard deviation") {
        double variance = mu*mu*exp(pow(log(sigma),2))*(exp(pow(log(sigma),2)) -1);
        REQUIRE(test_pdf.StdDev() == Catch::Approx(pow(variance,0.5)).epsilon(0.01));

    }

    SECTION("Scale PDF") {
        Vector_1D temp_pdf = test_pdf.getPDF();
        test_pdf.scalePdf(2.0);
        Vector_1D scaled_pdf = test_pdf.getPDF();
        for (int i = 0; i < nBins ; i ++) {
            REQUIRE(scaled_pdf[i] == Catch::Approx(2.0*temp_pdf[i]));
        }
        // Ensure we can map the pdf back to its original values
        test_pdf.scalePdf(0.5);
        scaled_pdf = test_pdf.getPDF();
        for (int i = 0; i < nBins ; i ++) {
            REQUIRE(scaled_pdf[i] == Catch::Approx(temp_pdf[i]));
        }
    }

    SECTION("Update PDF") {
        Vector_1D all_ones(nBins);
        std::fill(all_ones.begin(), all_ones.end(), 1.0);
        test_pdf.updatePdf(all_ones);
        Vector_1D temp_pdf = test_pdf.getPDF(); 
        for (int i = 0; i < nBins ; i ++) {
            REQUIRE(temp_pdf[i] == Catch::Approx(1.0));
        }

    }

    SECTION("Coagulation, Smoluchowski monodisperse analytical solution") {
        // Test against figure 15.2 from Jacobson (2005)
        double T = 298.;
        // double p = 101325.;
        // double rho = 1.0e3;
        double n0 = 10.e6; 
        double r0 = 3.0e-9; 
        double V0 = (4/3.0)*physConst::PI*pow(r0, 3);
        // double h = 12.0 * 3600;
        double kb = 1.38064852e-23;
        double eta = 1.8235e-5;
        double beta = 1.0e6 * 8*kb*T/(3*eta);
        int N = 100;

        Vector_1D bin_numbers(N);
        Vector_1D volumes(N);
        Vector_1D bin_centers(N);
        Vector_1D bin_edges(N+1);
        Vector_1D pdf(N); 
        std::fill(pdf.begin(), pdf.end(), 0);
        
        double V_rat = 1.2;

        for (int i = 0; i < N ; i++) {
            bin_numbers[i] = i;
            volumes[i] = V0 * pow(V_rat, i);
            bin_centers[i] = pow(0.75*volumes[i]/ physConst::PI, 1.0/3.0);
            bin_edges[i] = pow(0.75*(volumes[i]-V0/2.0)/ physConst::PI, 1.0/3.0);
            bin_edges[i+1] = pow(0.75*(volumes[i]+V0/2.0)/ physConst::PI, 1.0/3.0);

        }


        Vector_2D beta_vec( N, Vector_1D( N, 0.0E+00 ) );
        Vector_3D f_vec(N, Vector_2D(N, Vector_1D(N, 0.0)));
        double Vsum, Vfraction;
        int low_idx;  
        for (int i = 0; i < N ; i++) {
            for (int j = 0; j < N ; j++) {
                beta_vec[i][j] = beta;
                Vsum = volumes[i] + volumes[j];
                low_idx = std::lower_bound(volumes.begin(), volumes.end(), Vsum) - volumes.begin() - 1;
            
                if (low_idx == N-1) {
                    f_vec[i][j][N-1] = 1.0;
                }
                else {
                    Vfraction = volume_fraction(Vsum, volumes[low_idx], volumes[low_idx+1]);
                    f_vec[i][j][low_idx] = Vfraction;
                    f_vec[i][j][low_idx+1] = 1-Vfraction;
                }

            }
        }

        pdf[0] = n0;

        Aerosol aerosol(bin_centers, bin_edges, 0.0, mu, sigma);
        aerosol.updatePdf(pdf);

        int N_steps = 66;

        for (int i = 0; i < N_steps*6; i++) {
            aerosol.Coagulate(100., beta_vec, f_vec);
        }
        
        // Convert result to dn/dlog(diameter)
        Vector_1D result = aerosol.getPDF();
        for (int i = 0; i < N ; i++) {
            result[i] /= log10(2*bin_edges[i+1]) - log10(2*bin_edges[i]);
        }
        
        // Check first entro to be non-zero but also smaller than the initial condition
        REQUIRE(result[0] < 1.0e5);
        REQUIRE(result[0] > 0.0);

        // Check maximum
        double max;
        max = *std::max_element(result.begin(), result.end());
        REQUIRE(max == Catch::Approx(2.0e5).epsilon(0.5));

        // Check decay beyond maximum
        low_idx = std::lower_bound(bin_centers.begin(), bin_centers.end(), 0.6*0.03e-6) - result.begin() -1;
        REQUIRE(result[low_idx] < 10.0);
    }

}