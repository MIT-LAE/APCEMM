#include "Structure.h"

//using namespace std;

//double pSat_H2Ol( double T );
//double pSat_H2Os( double T );

Solution::Solution( const int nVar, const int n_x, const int n_y ) : \
        nVariables( nVar ), size_x( n_x ), size_y( n_y )
{
    /* Constructor */
}

Solution::~Solution()
{
    /* Destructor */
}

void Solution::Clear ( std::vector<std::vector<double> >& vector_2D )
{
    for ( unsigned int i = 0; i < vector_2D.size(); i++ ) {
        vector_2D[i].clear();
    }
    vector_2D.clear();

}

void Solution::SetShape ( std::vector<std::vector<double> >& vector_2D, unsigned int n_x, unsigned int n_y, double value )
{
    Clear( vector_2D );

    /* Dimensions are transposed! */
    for ( unsigned int i = 0; i < n_y; i++ ) {
        vector_2D.push_back( std::vector<double>( n_x, value ) );
    }

}

void Solution::SetToValue( std::vector<std::vector<double> >& vector_2D, double value )
{
    for ( unsigned int i = 0; i < vector_2D.size(); i++ ) {
        for ( unsigned int j = 0; j < vector_2D[0].size(); j++ ) {
            vector_2D[i][j] = (double) value;
        }
    }

}

void Solution::Print( std::vector<std::vector<double> >& vector_2D, unsigned int i_max, unsigned int j_max )
{
    for ( unsigned int i = 0; i < i_max; i++ ) {
        for ( unsigned int j = 0; j < j_max; j++ ) {
            std::cout << vector_2D[i][j];
        }
        std::cout << "" << std::endl;
    }
}
