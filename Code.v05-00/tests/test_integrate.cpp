#include <cmath>

#include <boost/numeric/odeint.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "EPM/Models/Original/Indexes.hpp"
#include "EPM/Models/Original/Integrate.hpp"
#include "EPM/Models/Original/StateObserver.hpp"
#include "EPM/Utils/Physics.hpp"
#include "Util/PhysFunction.hpp"

using namespace EPM::Utils::Physics;

TEST_CASE("EPM integrate", "[single-file]") {

    SECTION("Integrate") {


    }

    SECTION("RunMicrophysics") {


    }


    SECTION("Entrainment rate") {
        // Data from Karcher 1995 Figure 1
        // Act more like order of magnitude tests
        REQUIRE(entrainmentRate(0.0144) == Catch::Approx(40.82).epsilon(0.5));
        REQUIRE(entrainmentRate(0.311) == Catch::Approx(2.488).epsilon(0.5));
        REQUIRE(entrainmentRate(21.54) == Catch::Approx(0.01789).epsilon(0.5));
        REQUIRE(entrainmentRate(377.67) == Catch::Approx(0.00328).epsilon(0.5));
        REQUIRE(entrainmentRate(14762.3) == Catch::Approx(0.00006645468400721386).epsilon(0.5));
    }

    SECTION("Deposition rate") {


    }

    SECTION("odeRHS") {

    }

    SECTION("isFreezable") {

    }

    SECTION("Condensation rate") {

    }


}



namespace {
    using namespace EPM::Models::OriginalImpl;

    /* Same stepper as the EPM outer loop in Integrate.cpp */
    using epm_error_stepper = boost::numeric::odeint::runge_kutta_fehlberg78<Vector_1D>;

    /* Toy EPM state: T and P stay constant, and RH_w follows a parabola:
     *   RH_w(t) = rhEnd + (rhPeak - rhEnd) * 4 t (tEnd - t) / tEnd^2
     * RH_w equals rhEnd at t = 0 and t = tEnd, and rhPeak at t = tEnd / 2.
     */

    /* Water mixing ratio [-] at time t: satMixingRatio * RH_w(t) */
    double h2oProfile(double t, double satMixingRatio, double rhEnd, double rhPeak, double tEnd) {
        return satMixingRatio * (rhEnd + (rhPeak - rhEnd) * 4.0 * t * (tEnd - t) / (tEnd * tEnd));
    }

    /* Time derivative of h2oParabola [1/s] */
    double h2oProfileDerivative(double t, double satMixingRatio, double rhEnd, double rhPeak, double tEnd) {
        return satMixingRatio * (rhPeak - rhEnd) * 4.0 * (tEnd - 2.0 * t) / (tEnd * tEnd);
    }

    /* Toy EPM RHS: only the water mixing ratio changes, along h2oParabola */
    struct WaterParabolaRHS {
        double satMixingRatio, rhEnd, rhPeak, tEnd;

        void operator()(const Vector_1D &x, Vector_1D &dxdt, double t) const {
            std::fill(dxdt.begin(), dxdt.end(), 0.0);
            dxdt[EPM_ind_H2O] = h2oProfileDerivative(t, satMixingRatio, rhEnd, rhPeak, tEnd);
        }
    };

    double rhw(const Vector_1D &x) {
        return x[EPM_ind_H2O] * x[EPM_ind_P] / physFunc::pSat_H2Ol(x[EPM_ind_T]);
    }

    /* Integrate one outer interval [0, 10] s with the toy RHS, as the EPM
     * outer loop does (dt0 = interval / 100). Record only the endpoints and
     * return whether the observer saw saturation. */
    bool runWaterParabola(double rhEnd, double rhPeak, StateHistory &history) {
        const double T = 220.0;
        const double P = 2.65e4;
        const double tEnd = 10.0;
        const double satMixingRatio = physFunc::pSat_H2Ol(T) / P;
        const WaterParabolaRHS rhs{satMixingRatio, rhEnd, rhPeak, tEnd};

        Vector_1D x(13, 0.0);
        x[EPM_ind_T] = T;
        x[EPM_ind_P] = P;
        x[EPM_ind_H2O] = h2oProfile(0.0, satMixingRatio, rhEnd, rhPeak, tEnd);

        bool reached = false;
        const WaterSaturationObserver observer(reached);

        history.record(x, 0.0);
        // Take a single outer step from 0 to tEnd
        boost::numeric::odeint::integrate_adaptive(
            boost::numeric::odeint::make_controlled<epm_error_stepper>(1.0E-07, 1.0E-05),
            rhs, x, 0.0, tEnd, tEnd / 100.0, observer);
        history.record(x, tEnd);
        return reached;
    }
}

TEST_CASE("EPM water saturation observer", "[epm]") {
    SECTION("Flags a saturation peak between two recorded endpoints") {
        /* RH_w is 0.9 at both endpoints and >= 1 for t in about
         * [0.26, 9.74] s. At least one observed step is saturated. */
        StateHistory history({}, "");
        const bool reached = runWaterParabola(0.9, 1.9, history);

        /* The endpoints are not saturated */
        REQUIRE(history.states().size() == 2);
        for (const auto &x : history.states())
            REQUIRE(rhw(x) < 1.0);
        
        /* but at least one adaptive steps inside the interval is which
        we should have detected. */
        REQUIRE(reached);
    }

    SECTION("Stays false when the plume never saturates") {
        StateHistory history({}, "");
        REQUIRE_FALSE(runWaterParabola(0.5, 0.8, history));
    }
}

TEST_CASE("EPM outer time grid", "[epm]") {
    SECTION("EPM grid has an exact 180 s point") {
        const Vector_1D grid = makeEPMTimeGrid(301, 1.0E-04, 2.0E+03, 180.0);

        REQUIRE(grid.size() == 302);
        REQUIRE(grid.front() == 1.0E-04);
        REQUIRE(grid.back() == 2.0E+03);
        REQUIRE(std::count(grid.begin(), grid.end(), 180.0) == 1);
    }

}
