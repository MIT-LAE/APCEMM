#ifndef PLUMEMODELUTILS_H
#define PLUMEMODELUTILS_H

#include <vector>
#include <cmath>
#include <Core/Parameters.hpp>

namespace PlumeModelUtils {
    double UpdateTime( double time, const double tStart, \
                    const double DYN_DT, double& nextTimeStep );

    std::vector<double> BuildTime( const double tStart, const double tEnd, \
                                const double DYN_DT );


    void AdvGlobal( const double time, const double T_UPDRAFT, \
                    const double V_UPDRAFT,                    \
                    double &v_x, double &v_y,                  \
                    double &dTrav_x, double &dTrav_y );
                    
    void DiffParam( const double time, double &d_x, double &d_y, \
                const double D_X, const double D_Y );
}



#endif