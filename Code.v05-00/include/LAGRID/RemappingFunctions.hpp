#ifndef LAGRID_REMAPPINGFUNCTIONS_H
#define LAGRID_REMAPPINGFUNCTIONS_H

#include "LAGRID/FreeCoordBoxGrid.hpp"
#include <algorithm>
#include <type_traits>
#include <functional>

namespace LAGRID {
    
    struct Remapping {
        double x0;
        double y0;
        double dx;
        double dy;
        int nx;
        int ny;
        Remapping(double x0, double y0, double dx, double dy, int nx, int ny)
        : x0(x0),
          y0(y0),
          dx(dx),
          dy(dy),
          nx(nx),
          ny(ny)
        {

        }
    };

    struct twoDGridVariable {
        twoDGridVariable() = delete;

        //A little bit of template / universal reference magic to minimize the # of copy operations on these huge vectors.
        //Enforces that Vec2D is of type Vector_2D and that Vec1D is of type Vec1D. And then we need the decay on the types because 
        //with a forwarding reference, passing in a Vector_2D by value gets it deduced as a Vector_2D&
        //-Michael
        template <typename Vec2D, typename Vec1D,
                    typename = std::enable_if_t<std::is_same_v<std::decay_t<Vec2D>, Vector_2D>>,
                    typename = std::enable_if_t<std::is_same_v<std::decay_t<Vec1D>, Vector_1D>>
                 >
        twoDGridVariable(Vec2D&& phi, Vec1D&& xCoords, Vec1D&& yCoords) :
            phi(std::forward<Vector_2D>(phi)),
            xCoords(std::forward<Vector_1D>(xCoords)),
            yCoords(std::forward<Vector_1D>(yCoords)),
            dx(this->xCoords[1] - this->xCoords[0]),
            dy(this->yCoords[1] - this->yCoords[0])
        {
        }
        Vector_2D phi;
        Vector_1D xCoords;
        Vector_1D yCoords;
        double dx;
        double dy;
        void addBuffer(double bufLen_left, double bufLen_right, double bufLen_top, double bufLen_bot, double fillValue);
    };

    inline double coveredArea(const Remapping& remapping, const MassBox& b, int i, int j) {
        //Breaks if the box and the gridcell don't overlap at all.
        double gridCell_topLeftX = remapping.x0 + remapping.dx * i;
        double gridCell_topLeftY = remapping.y0 + remapping.dy * (j + 1);
        double gridCell_botRightX = remapping.x0 + remapping.dx * (i + 1);
        double gridCell_botRightY = remapping.y0 + remapping.dy * j;

        //X,Y coords for the portion of the gridcell the box covers
        double boxCover_topLeftX = std::max(b.topLeftX, gridCell_topLeftX);
        double boxCover_topLeftY = std::min(b.topLeftY, gridCell_topLeftY);
        double boxCover_botRightX = std::min(b.botRightX, gridCell_botRightX);
        double boxCover_botRightY = std::max(b.botRightY, gridCell_botRightY);

        return (boxCover_botRightX - boxCover_topLeftX) * (boxCover_topLeftY - boxCover_botRightY);
    }

    FreeCoordBoxGrid rectToBoxGrid(double dy_old, const Vector_1D& dy_new, double dx_old, double x0_old, double y0_new, const Vector_2D& phi_old, const vector<vector<int>>& mask); 

    double diffusionLossFunctionExact(const FreeCoordBoxGrid& boxGrid, const Remapping& remapping);
    double diffusionLossFunctionBoundaryEstimate(const FreeCoordBoxGrid& boxGrid, const Remapping& remapping);
    twoDGridVariable mapToStructuredGrid(const FreeCoordBoxGrid& boxGrid, const Remapping& remapping);
    twoDGridVariable getUnusedFraction(const FreeCoordBoxGrid& boxGrid, const Remapping& remapping);

    Vector_2D initVarToGrid(double mass, const Vector_1D& xEdges, const Vector_1D& yEdges,
                            std::function<double(double, double)> weightFunction, double logBinRatio = 1 );

    Vector_2D initVarToGridGaussian(double mass, const Vector_1D& xEdges, const Vector_1D& yEdges, double x0, double y0,
                            double sigmaX, double sigmaY, double logBinRatio = 1);

    Vector_2D initVarToGridBimodalY(double mass, const Vector_1D& xEdges, const Vector_1D& yEdges, double x0, double y0,
                                    double width, double depth, double logBinRatio = 1);
}

#endif
