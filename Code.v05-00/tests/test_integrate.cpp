#include <cmath>

#include <boost/numeric/odeint.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "EPM/Models/Original/Indexes.hpp"
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

    /* Toy EPM state: T and P stay constant, and the water mixing ratio
     * follows RH_w(t) = rhEnd + (rhPeak - rhEnd) * sin^2(pi t / tEnd).
     * RH_w equals rhEnd at both ends of [0, tEnd] and peaks at mid-interval. */
    struct WaterHumpRHS {
        double satMixingRatio, rhEnd, rhPeak, tEnd;

        double h2o(double t) const {
            const double s = std::sin(M_PI * t / tEnd);
            return satMixingRatio * (rhEnd + (rhPeak - rhEnd) * s * s);
        }

        void operator()(const Vector_1D &x, Vector_1D &dxdt, double t) const {
            std::fill(dxdt.begin(), dxdt.end(), 0.0);
            dxdt[EPM_ind_H2O] = satMixingRatio * (rhPeak - rhEnd) * M_PI / tEnd * std::sin(2.0 * M_PI * t / tEnd);
        }
    };

    double rhw(const Vector_1D &x) {
        return x[EPM_ind_H2O] * x[EPM_ind_P] / physFunc::pSat_H2Ol(x[EPM_ind_T]);
    }

    /* Integrate one outer interval [0, 10] s with the toy RHS, as the EPM
     * outer loop does (dt0 = interval / 100). Record only the endpoints and
     * return whether the observer saw saturation. */
    bool runWaterHump(double rhEnd, double rhPeak, StateHistory &history) {
        const double T = 220.0;
        const double P = 2.65e4;
        const double tEnd = 10.0;
        const WaterHumpRHS rhs{physFunc::pSat_H2Ol(T) / P, rhEnd, rhPeak, tEnd};

        Vector_1D x(13, 0.0);
        x[EPM_ind_T] = T;
        x[EPM_ind_P] = P;
        x[EPM_ind_H2O] = rhs.h2o(0.0);

        bool reached = false;
        const WaterSaturationObserver observer(reached);

        history.record(x, 0.0);
        boost::numeric::odeint::integrate_adaptive(
            boost::numeric::odeint::make_controlled<epm_error_stepper>(1.0E-07, 1.0E-05),
            rhs, x, 0.0, tEnd, tEnd / 100.0, observer);
        history.record(x, tEnd);
        return reached;
    }
}

TEST_CASE("EPM water saturation observer", "[epm]") {
    SECTION("Flags a saturation peak between two recorded endpoints") {
        /* RH_w >= 1 for t in about [1, 9] s. The controller grows the step
         * by at most 5x from dt0 = 0.1 s, so it cannot jump from below 1 s
         * to beyond 9 s: at least one observed step is saturated. */
        StateHistory history({}, "");
        const bool reached = runWaterHump(0.9, 1.9, history);

        /* The endpoints never reach saturation ... */
        REQUIRE(history.states().size() == 2);
        for (const auto &x : history.states())
            REQUIRE(rhw(x) < 1.0);
        /* ... but the adaptive steps inside the interval do. odeint copies
         * the observer by value, so this also checks that the copies write
         * to the caller's flag. */
        REQUIRE(reached);
    }

    SECTION("Stays false when the plume never saturates") {
        StateHistory history({}, "");
        REQUIRE_FALSE(runWaterHump(0.5, 0.8, history));
    }
}
