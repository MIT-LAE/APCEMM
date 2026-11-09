#ifndef EPM_MODELS_ORIGINALIMPL_INTEGRATE_H_INCLUDED
#define EPM_MODELS_ORIGINALIMPL_INTEGRATE_H_INCLUDED

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "Util/ForwardDecl.hpp"

namespace EPM::Models::OriginalImpl
{
    /* Outer EPM time grid: nLog log-spaced points from tStart to tEnd, plus
     * tInsert as an extra exact grid point (unless it already is one). */
    inline Vector_1D makeEPMTimeGrid(UInt nLog, double tStart, double tEnd, double tInsert)
    {
        if ( nLog < 2 || !( tStart > 0.0 && tStart < tEnd && tStart <= tInsert && tInsert <= tEnd ) )
            throw std::invalid_argument("makeEPMTimeGrid: need nLog >= 2 and 0 < tStart <= tInsert <= tEnd, tStart < tEnd");

        Vector_1D grid(nLog);
        const double log10_tStart = std::log10(tStart);
        const double log10_tEnd = std::log10(tEnd);
        const double inv_denominator = 1.0 / ( nLog - 1.0 );
        for ( UInt i = 0; i < nLog; i++ )
            grid[i] = std::pow( 10.0, log10_tStart + i * ( log10_tEnd - log10_tStart ) * inv_denominator );

        // Ensure that start and end times are exact (can be slight off due to rounding in log/exp)
        grid.front() = tStart;
        grid.back() = tEnd;

        // Add tInsert if it's not already in the grid
        auto const it = std::lower_bound(grid.begin(), grid.end(), tInsert);
        if ( *it != tInsert )
            grid.insert(it, tInsert);
        return grid;
    }
}

#endif
