/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Reshape Program File                                             */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : Reshape.cpp                               */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <iostream>
#include <vector>

std::vector<std::vector<double> > Reshape_Vector( std::vector<std::vector<double> >& vector_2D, int n_x, int n_y )
{
    std::vector<std::vector<double> > output;
    int size_x = vector_2D.size();
    int size_y = vector_2D[0].size();

    if ( n_x * n_y != size_x * size_y ) {
        std::cout << "Invalid dimensions specified" << std::endl;
        return output;
    }
    else {
        int counter;
        for ( counter = 0; counter < n_x; counter++ )
            output.push_back(std::vector<double>(n_y));

        int counter_x = 0;
        int counter_y = 0;
        int orig_counter_x = 0;
        int orig_counter_y = 0;
        for ( counter = 0; counter < n_x * n_y; counter++ ) {
            counter_x = counter%n_x;
            //counter_y = counter%n_y;
            orig_counter_x = counter%size_x;
            //orig_counter_y = counter%size_y;
            output[counter_x][counter_y] = vector_2D[orig_counter_x][orig_counter_y];
            //std::cout << orig_counter_x << ", " << orig_counter_y << ": " << vector_2D[orig_counter_x][orig_counter_y] << " -> " << counter_x << ", " << counter_y << std::endl;
            if ( counter_x == n_x - 1 )
                counter_y += 1;
            if ( orig_counter_x == size_x - 1 )
                orig_counter_y += 1;
//            if ( counter_y == n_y - 1 )
//                counter_x += 1;
//            if ( orig_counter_y == size_y - 1 )
//                orig_counter_x += 1;
        }
//        int k, l;
//        for ( i = 0; i < n_x; i++ ) {
//            for ( j = 0; j < n_y; j++ ) {
//                k = i%size_x;
//                l = j%size_y;
//                output[i][j] = vector_2D[k][l];
//            }
//        }
    }

    return output;

} /* End of Reshape_Vector */
