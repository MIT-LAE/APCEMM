#include "AIM/Aerosol.hpp"
#include "Util/ForwardDecl.hpp"
#include "Util/PhysConstant.hpp"
#include "Util/PhysFunction.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

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
        double n0 = 10.e6; // [#/cm^3]
        double r0 = 3.0e-9; 
        double V0 = (4/3.0)*physConst::PI*pow(r0, 3);
        double kb = 1.38064852e-23;
        double eta = 1.8235e-5;
        double beta = 1.0e6 * 8*kb*T/(3*eta); // Constant kernel [cm^3/s]
        int N = 100;
        double V_rat = 1.2;

        // Volume-ratio grid. The edge volumes are 2 V_i / (1 + V_rat), so
        // that the mean of the two edge volumes of bin i, which Aerosol uses
        // as the bin volume, is exactly V_i.
        Vector_1D volumes(N);
        Vector_1D bin_centers(N);
        Vector_1D bin_edges(N+1);
        for (int i = 0; i < N ; i++) {
            volumes[i] = V0 * pow(V_rat, i);
            bin_centers[i] = cbrt(0.75*volumes[i] / physConst::PI);
        }
        for (int i = 0; i <= N ; i++) {
            double V_edge = 2.0 * V0 * pow(V_rat, i) / (1.0 + V_rat);
            bin_edges[i] = cbrt(0.75*V_edge / physConst::PI);
        }

        Vector_2D beta_vec( N, Vector_1D( N, beta ) );
        Vector_3D f_vec(N, Vector_2D(N, Vector_1D(N, 0.0)));
        double Vsum, Vfraction;
        int low_idx;  
        for (int i = 0; i < N ; i++) {
            for (int j = 0; j < N ; j++) {
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

        // pdf is a density dN/dln(r): put n0 particles in the first bin
        const double dlnr = log(bin_edges[1]) - log(bin_edges[0]);
        Vector_1D pdf(N, 0.0);
        pdf[0] = n0 / dlnr;

        Aerosol aerosol(bin_centers, bin_edges, 0.0, mu, sigma);
        aerosol.updatePdf(pdf);
        REQUIRE(aerosol.Moment(0) == Catch::Approx(n0).epsilon(1e-12));
        const double M3_0 = aerosol.Moment(3);

        const double dt = 100.;
        const int N_steps = 66*6;
        for (int i = 0; i < N_steps; i++) {
            aerosol.Coagulate(dt, beta_vec, f_vec);
        }

        // Total number against the analytical solution for a constant
        // kernel: N(t) = N0 / (1 + beta N0 t / 2)
        const double t = N_steps * dt;
        const double N_exact = n0 / (1.0 + 0.5 * beta * n0 * t);
        REQUIRE(aerosol.Moment(0) == Catch::Approx(N_exact).epsilon(1e-2));

        // Coagulation conserves volume
        REQUIRE(aerosol.Moment(3) == Catch::Approx(M3_0).epsilon(1e-12));

        // Convert result to dN/dlog10(D): dN/dlog10(D) = ln(10) * dN/dln(r)
        Vector_1D result = aerosol.getPDF();
        for (int i = 0; i < N ; i++) {
            result[i] *= log(10.0);
        }
        
        // Check first entry to be non-zero but also smaller than the initial condition
        REQUIRE(result[0] < 1.0e5);
        REQUIRE(result[0] > 0.0);

        // Check maximum
        double max;
        max = *std::max_element(result.begin(), result.end());
        REQUIRE(max == Catch::Approx(2.0e5).epsilon(0.5));

        // Check decay beyond the maximum (D = 0.03 um): dN/dlog10(D) < 10 at D = 0.08 um
        low_idx = std::lower_bound(bin_centers.begin(), bin_centers.end(), 0.04e-6) - bin_centers.begin() -1;
        REQUIRE(result[low_idx] < 10.0);
    }

}

TEST_CASE("Aerosol Flexible Bin Grid Construction", "[aerosol][bins]") {
    const double r_low = 5.0e-8;  // 0.05 um
    const double r_hig = 8.0e-5;  // 80.0 um

    SECTION("Default Grid Parameters (VRAT = 1.50)") {
        const double vrat = 1.50;
        const UInt expected_nbins = static_cast<UInt>(std::floor(1 + 3.0 * log(r_hig / r_low) / log(vrat)));
        REQUIRE(expected_nbins == 55);

        const double rrat = cbrt(vrat);
        Vector_1D edges(expected_nbins + 1);
        Vector_1D centers(expected_nbins);
        edges[0] = r_low;
        for (UInt i = 1; i <= expected_nbins; i++) {
            edges[i] = edges[i - 1] * rrat;
        }
        for (UInt i = 0; i < expected_nbins; i++) {
            centers[i] = 0.5 * (edges[i] + edges[i + 1]);
        }

        REQUIRE(edges.front() == Catch::Approx(5.0e-8));
        REQUIRE(edges.back() >= 8.0e-5);
    }

    SECTION("VRAT Spectrum Geometric Invariance (1.2 to 3.0)") {
        const std::vector<double> vrats = {1.2, 1.35, 1.5, 1.8, 2.0, 2.1, 2.4, 2.7, 3.0};
        const std::vector<UInt> expected_counts = {
            static_cast<UInt>(std::floor(1 + 3.0 * log(r_hig / r_low) / log(1.2))),  // ~122
            static_cast<UInt>(std::floor(1 + 3.0 * log(r_hig / r_low) / log(1.35))), // 75
            static_cast<UInt>(std::floor(1 + 3.0 * log(r_hig / r_low) / log(1.5))),  // 55
            static_cast<UInt>(std::floor(1 + 3.0 * log(r_hig / r_low) / log(1.8))),  // 38
            static_cast<UInt>(std::floor(1 + 3.0 * log(r_hig / r_low) / log(2.0))),  // 32
            static_cast<UInt>(std::floor(1 + 3.0 * log(r_hig / r_low) / log(2.1))),  // 30
            static_cast<UInt>(std::floor(1 + 3.0 * log(r_hig / r_low) / log(2.4))),  // 26
            static_cast<UInt>(std::floor(1 + 3.0 * log(r_hig / r_low) / log(2.7))),  // 23
            static_cast<UInt>(std::floor(1 + 3.0 * log(r_hig / r_low) / log(3.0)))   // 21
        };

        for (size_t k = 0; k < vrats.size(); k++) {
            const double vrat = vrats[k];
            const UInt nbins = static_cast<UInt>(std::floor(1 + 3.0 * log(r_hig / r_low) / log(vrat)));
            REQUIRE(nbins == expected_counts[k]);

            const double rrat = cbrt(vrat);
            Vector_1D edges(nbins + 1);
            Vector_1D centers(nbins);
            edges[0] = r_low;
            for (UInt i = 1; i <= nbins; i++) {
                edges[i] = edges[i - 1] * rrat;
                // Ratio test for each adjacent edge
                REQUIRE(edges[i] / edges[i - 1] == Catch::Approx(rrat).epsilon(1e-12));
            }
            for (UInt i = 0; i < nbins; i++) {
                centers[i] = 0.5 * (edges[i] + edges[i + 1]);
            }

            REQUIRE(edges.front() == Catch::Approx(r_low));
            REQUIRE(edges.back() >= r_hig);
        }
    }

    SECTION("Physical Moment Conservation across VRAT Spectrum") {
        const double nPart = 1.0e12; // #/cm3
        const double mu = 1.0e-6;    // 1.0 um
        const double sigma = 1.6;

        const std::vector<double> vrats = {1.2, 1.5, 1.8, 2.0, 2.1, 2.4, 2.7, 3.0};
        const double analytic_m0 = nPart;
        const double analytic_m3 = nPart * pow(mu, 3.0) * exp((9.0 / 2.0) * pow(log(sigma), 2.0));

        for (double vrat : vrats) {
            const UInt nbins = static_cast<UInt>(std::floor(1 + 3.0 * log(r_hig / r_low) / log(vrat)));
            const double rrat = cbrt(vrat);
            Vector_1D edges(nbins + 1);
            Vector_1D centers(nbins);
            edges[0] = r_low;
            for (UInt i = 1; i <= nbins; i++) {
                edges[i] = edges[i - 1] * rrat;
            }
            for (UInt i = 0; i < nbins; i++) {
                centers[i] = 0.5 * (edges[i] + edges[i + 1]);
            }

            Aerosol aer(centers, edges, nPart, mu, sigma, "lognormal");

            // Zeroth moment (total particle number) should be preserved to < 1%
            double m0 = aer.Moment(0);
            REQUIRE(m0 == Catch::Approx(analytic_m0).epsilon(0.01));

            // Third moment (volume / mass proportional) should be preserved to < 5%
            double m3 = aer.Moment(3);
            REQUIRE(m3 == Catch::Approx(analytic_m3).epsilon(0.05));
        }
    }

    SECTION("Size-Sensitive APC_Scheme Growth Continuity Across Bin Boundaries") {
        const double vrat = 1.80;
        const UInt nbins = static_cast<UInt>(std::floor(1 + 3.0 * log(r_hig / r_low) / log(vrat)));
        const double rrat = cbrt(vrat);
        Vector_1D edges(nbins + 1);
        Vector_1D centers(nbins);
        edges[0] = r_low;
        for (UInt i = 1; i <= nbins; i++) {
            edges[i] = edges[i - 1] * rrat;
        }
        for (UInt i = 0; i < nbins; i++) {
            centers[i] = 0.5 * (edges[i] + edges[i + 1]);
        }

        // Create a 2x1 grid where cell (0,0) has crystals in bin k just below V_edge[k+1]
        // and cell (0,1) has the same number of crystals in bin k+1 just above V_edge[k+1].
        Grid_Aerosol gridAer(2, 1, centers, edges, 0.0, 1.0e-6, 1.5, "lognormal");
        const UInt k = 25; // ~15 um bin
        const Vector_1D& vEdges = gridAer.getBinVEdges();
        const Vector_1D& logEdges = gridAer.getLogBinEdges();
        const double vBoundary = vEdges[k + 1];
        const double nConc = 0.5; // #/cm^3

        Vector_3D pdf = gridAer.getPDF();
        for (UInt b = 0; b < nbins; ++b) {
            pdf[b][0][0] = 0.0;
            pdf[b][0][1] = 0.0;
        }
        pdf[k][0][0] = nConc / (logEdges[k + 1] - logEdges[k]);
        pdf[k + 1][0][1] = nConc / (logEdges[k + 2] - logEdges[k + 1]);
        gridAer.updatePdf(pdf);

        Vector_3D iceV(nbins, Vector_2D(1, Vector_1D(2, 0.0)));
        iceV[k][0][0] = nConc * (0.9999 * vBoundary);
        iceV[k + 1][0][1] = nConc * (1.0001 * vBoundary);
        gridAer.UpdateCenters(iceV, gridAer.getPDF());

        const double T_val = 217.0;
        const double P_val = 26500.0;
        const double pSat = physFunc::pSat_H2Os(T_val);
        const double h2o_init = 1.40 * pSat / (physConst::kB * 1.0e6 * T_val); // RHi = 140%
        Vector_2D H2O(1, Vector_1D(2, h2o_init));
        Vector_2D T(1, Vector_1D(2, T_val));
        Vector_1D P(1, P_val);

        Vector_2D volBefore = gridAer.TotalVolume();
        gridAer.Grow(10.0, H2O, T, P, 2, 0);
        Vector_2D volAfter = gridAer.TotalVolume();

        double dV_left = volAfter[0][0] - volBefore[0][0];
        double dV_right = volAfter[0][1] - volBefore[0][1];
        REQUIRE(dV_left > 0.0);
        // With moving-center radius in APC_Scheme, growth rate across the bin edge is continuous (< 0.1% diff, rather than ~22% jump)
        REQUIRE(dV_right == Catch::Approx(dV_left).epsilon(1e-3));
    }
}

