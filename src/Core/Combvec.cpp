/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Combvec Program File                                             */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : Combvec.cpp                               */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <iostream>
#include <vector>


std::vector<std::vector<double> > Copy_blocked( std::vector<std::vector<double> >& m, \
                                                int n );
std::vector<std::vector<double> > Copy_interleaved( std::vector<std::vector<double> >& m, \
                                                    int n );
std::vector<std::vector<double> > Reshape_Vector( std::vector<std::vector<double> >& vector_2D, int n_x, int n_y );

std::vector<std::vector<double> > CombVec( std::vector<double>& temperature_K, \
                                 std::vector<double>& pressure_Pa,   \
                                 std::vector<double>& relHumidity_w, \
                                 std::vector<double>& longitude_deg, \
                                 std::vector<double>& latitude_deg )
{
    std::vector<std::vector<double> > combinations;

    unsigned int counter = 1;
    unsigned int nCases = temperature_K.size();

    unsigned int i, j;

    std::vector<std::vector<double> > y, z;
    std::vector<std::vector<double> > u, v;

    y.push_back(std::vector<double>(temperature_K.size()));
    for ( i = 0; i < temperature_K.size(); i++ )
        y[0][i] = temperature_K[i];

    /* z = pressure_Pa */
    z.push_back(std::vector<double>(pressure_Pa.size()));
    for ( i = 0; i < pressure_Pa.size(); i++ )
        z[0][i] = pressure_Pa[i];

    u = Copy_blocked(y,z[0].size());
    v = Copy_interleaved(z,y[0].size());

    for ( i = 0; i < counter; i++ ) {
        y[i].clear();
    }
    z[0].clear();
    y.clear(); z.clear();

    nCases *= pressure_Pa.size();
    counter += 1;
    for ( i = 0; i < counter; i++ )
        y.push_back(std::vector<double>( nCases ));

    for ( i = 0; i < nCases; i++ ) {
        for ( j = 0; j < counter - 1; j++ ) {
            y[j][i] = u[j][i];
        }
        y[counter-1][i] = v[0][i];
    }

    for ( i = 0; i < counter - 1; i++ ) {
        u[i].clear();
    }
    u.clear();
    v[0].clear(); v.clear();

    /* z = relHumidity_w */
    z.push_back(std::vector<double>(relHumidity_w.size()));
    for ( i = 0; i < relHumidity_w.size(); i++ )
        z[0][i] = relHumidity_w[i];

    u = Copy_blocked(y,z[0].size());
    v = Copy_interleaved(z,y[0].size());

    for ( i = 0; i < counter; i++ ) {
        y[i].clear();
    }
    z[0].clear();
    y.clear(); z.clear();

    nCases *= relHumidity_w.size();
    counter += 1;
    for ( i = 0; i < counter; i++ )
        y.push_back(std::vector<double>( nCases ));

    for ( i = 0; i < nCases; i++ ) {
        for ( j = 0; j < counter - 1; j++ ) {
            y[j][i] = u[j][i];
        }
        y[counter-1][i] = v[0][i];
    }

    for ( i = 0; i < counter - 1; i++ ) {
        u[i].clear();
    }
    u.clear();
    v[0].clear(); v.clear();

    /* z = longitude_deg */
    z.push_back(std::vector<double>(longitude_deg.size()));
    for ( i = 0; i < longitude_deg.size(); i++ )
        z[0][i] = longitude_deg[i];

    u = Copy_blocked(y,z[0].size());
    v = Copy_interleaved(z,y[0].size());

    for ( i = 0; i < counter; i++ ) {
        y[i].clear();
    }
    z[0].clear();
    y.clear(); z.clear();

    nCases *= longitude_deg.size();
    counter += 1;
    for ( i = 0; i < counter; i++ )
        y.push_back(std::vector<double>( nCases ));

    for ( i = 0; i < nCases; i++ ) {
        for ( j = 0; j < counter - 1; j++ ) {
            y[j][i] = u[j][i];
        }
        y[counter-1][i] = v[0][i];
    }

    for ( i = 0; i < counter - 1; i++ ) {
        u[i].clear();
    }
    u.clear();
    v[0].clear(); v.clear();

    /* z = latitude_deg */
    z.push_back(std::vector<double>(latitude_deg.size()));
    for ( i = 0; i < latitude_deg.size(); i++ )
        z[0][i] = latitude_deg[i];

    u = Copy_blocked(y,z[0].size());
    v = Copy_interleaved(z,y[0].size());

    for ( i = 0; i < counter; i++ ) {
        y[i].clear();
    }
    z[0].clear();
    y.clear(); z.clear();

    nCases *= latitude_deg.size();
    counter += 1;
    for ( i = 0; i < counter; i++ )
        y.push_back(std::vector<double>( nCases ));

    for ( i = 0; i < nCases; i++ ) {
        for ( j = 0; j < counter - 1; j++ ) {
            y[j][i] = u[j][i];
        }
        y[counter-1][i] = v[0][i];
    }

    return y;

} /* End of CombVec */

std::vector<std::vector<double> > Copy_blocked( std::vector<std::vector<double> >& m, \
                                                int n )
{
    std::vector<std::vector<double> > b;
    const unsigned int mr = m.size();
    const unsigned int mc = m[0].size();
    unsigned int i, j, k;

    for ( i = 0; i < mr; i++ )
        b.push_back(std::vector<double>(mc*n));

//    unsigned int ind[mc];
    unsigned int *ind = new unsigned int[mc];

    for ( i = 0; i < mc; i++ )
        ind[i] = i;

    for ( j = 0; j < (n-1) * mc + 1; j+=mc ) {
        for ( i = 0; i < mr; i++ ) {
            for ( k = 0; k < mc; k++ ) {
                b[i][ind[k] + j] = m[i][k];
            }
        }
    }

    return b;
} /* End of Copy_blocked */

std::vector<std::vector<double> > Copy_interleaved( std::vector<std::vector<double> >& m, \
                                                    int n )
{
    std::vector<std::vector<double> > b;
    const unsigned int mr = m.size();
    const unsigned int mc = m[0].size();
    unsigned int i, j, k;

    for ( i = 0; i < mr*n; i++ )
        b.push_back(std::vector<double>(mc));

//    unsigned int ind[mr];
    unsigned int *ind = new unsigned int[mr];

    for ( i = 0; i < mr; i++ )
        ind[i] = i;

    for ( i = 0; i < (n-1) * mr + 1; i+=mr ) {
        for ( j = 0; j < mc; j++ ) {
            for ( k = 0; k < mr; k++ ) {
                b[ind[k] + i][j] = m[k][j]; //m[k][j]
            }
        }
    }

    return Reshape_Vector( b, mr, n*mc );

} /* End of Copy_interleaved */
