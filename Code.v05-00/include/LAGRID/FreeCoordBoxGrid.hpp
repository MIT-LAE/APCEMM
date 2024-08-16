#ifndef LAGRID_FREECOORDBOXGRID_H
#define LAGRID_FREECOORDBOXGRID_H

#include "Util/ForwardDecl.hpp"

namespace LAGRID {
    using std::vector;
    struct MassBox {
        MassBox(double tlx, double tly, double brx, double bry, double mass):
            topLeftX(tlx), topLeftY(tly), botRightX(brx), botRightY(bry), mass(mass)
        {
    
        }
        double topLeftX;
        double topLeftY;
        double botRightX;
        double botRightY;
        double mass; // mass = phi * dy * dx where phi is a concentration
        inline const double area() const {
            return (botRightX - topLeftX) * (topLeftY - botRightY);
        }
    };
    struct FreeCoordBoxGrid {
        FreeCoordBoxGrid() = delete;
        //Can apply a mask to the cutoff, and then can easily find the boundary nodes in O(N) time using the 0/1 mask.
        FreeCoordBoxGrid(const Vector_1D& dx, const Vector_1D& dy, const Vector_2D& phi, const Vector_1D& x0, double y0, double cutoff);
        FreeCoordBoxGrid(const Vector_1D& dx, const Vector_1D& dy, const Vector_2D& phi, const Vector_1D& x0, double y0, const vector<vector<int>>& mask);
        
        void setMinMaxCoords();
        std::vector<MassBox> boxes;
        std::vector<int> boundaryBoxIndices; //To speed up the remapping algorithm
        double minX;
        double minY;
        double maxX;
        double maxY;

            

    };
}

#endif