#ifndef EPM_MODELS_ORIGINALIMPL_INTEGRATE_H_INCLUDED
#define EPM_MODELS_ORIGINALIMPL_INTEGRATE_H_INCLUDED

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

#include <boost/numeric/odeint/integrate/integrate_adaptive.hpp>

#include "Util/ForwardDecl.hpp"
#include "EPM/Models/Original/StateObserver.hpp"

namespace EPM::Models::OriginalImpl
{
    /* Outer EPM time grid: nLog log-spaced points from tStart to tEnd, plus
     * tInsert as an extra exact grid point (unless it already is one). The
     * first and last points are exactly tStart and tEnd; pow/log10 rounding
     * would otherwise leave them a few ULP off. tInsert may equal tStart or
     * tEnd, and then adds no point. The EPM uses tInsert = 180 s, so that the
     * 3-minute hand-off falls exactly on a grid point. */
    inline Vector_1D makeEPMTimeGrid(UInt nLog, double tStart, double tEnd, double tInsert)
    {
        if ( nLog < 2 || !( tStart > 0.0 && tStart < tEnd && tStart <= tInsert && tInsert <= tEnd ) )
            throw std::invalid_argument("makeEPMTimeGrid: need nLog >= 2 and 0 < tStart <= tInsert <= tEnd, tStart < tEnd");

        Vector_1D grid(nLog);
        const double log10_tStart = std::log10(tStart);
        const double log10_tEnd = std::log10(tEnd);
        const double inv_denominator = 1/double( nLog - 1.0 );
        for ( UInt i = 0; i < nLog; i++ )
            grid[i] = std::pow( 10.0, log10_tStart + i * ( log10_tEnd - log10_tStart ) * inv_denominator );
        grid.front() = tStart;
        grid.back() = tEnd;

        /* grid.back() == tEnd >= tInsert, so it is never end() */
        auto const it = std::lower_bound(grid.begin(), grid.end(), tInsert);
        if ( *it != tInsert )
            grid.insert(it, tInsert);
        return grid;
    }

    /* Advance the EPM state x over one outer interval [t0, t1], then record
     * the end state x(t1) at time t1. The observer sees every adaptive step
     * but does not touch the history. Returns the number of adaptive steps. */
    template <class Stepper, class System, class Observer>
    std::size_t advanceInterval(
        Stepper stepper, System system, Vector_1D &x,
        double t0, double t1, double dt0,
        Observer observer, StateHistory &history)
    {
        const std::size_t nSteps = boost::numeric::odeint::integrate_adaptive(
            stepper, system, x, t0, t1, dt0, observer);
        history.record(x, t1);
        return nSteps;
    }
}

#endif
