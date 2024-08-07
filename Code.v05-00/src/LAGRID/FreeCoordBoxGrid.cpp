#include "LAGRID/FreeCoordBoxGrid.hpp"
namespace LAGRID {
    FreeCoordBoxGrid::FreeCoordBoxGrid(const Vector_1D& dx, const Vector_1D& dy, const Vector_2D& phi, const Vector_1D& x0, double y0, double cutoff) {
        auto isInteriorCell = [&](std::size_t i, std::size_t j) -> bool {
            bool edgeIndex = (i == 0 || j == 0 || i == dx.size() - 1 || j == dy.size() - 1);
            return !edgeIndex && (phi[j+1][i] > cutoff && phi[j-1][i] > cutoff && 
                                  phi[j][i+1] > cutoff && phi[j][i-1] > cutoff);
        };

        double curr_y = y0; // Y coordinate of bottom edge of the cell row
        for(std::size_t j = 0; j < dy.size(); j++) {
            double curr_x = x0[j]; // X coord of left edge of cell column
            for(std::size_t i = 0; i < dx.size(); i++) {
                if(phi[j][i] < cutoff) {
                    curr_x += dx[j];
                    continue;
                } 
                boxes.emplace_back(curr_x, curr_y + dy[j], curr_x + dx[j], curr_y, phi[j][i] * dx[j] * dy[j]);
                if(!isInteriorCell(i, j)) {
                    boundaryBoxIndices.push_back(boxes.size() - 1);
                }
                curr_x += dx[j];
            }
            curr_y += dy[j];
        }
        setMinMaxCoords();
    }

    FreeCoordBoxGrid::FreeCoordBoxGrid(const Vector_1D& dx, const Vector_1D& dy, const Vector_2D& phi, const Vector_1D& x0, double y0, const vector<vector<int>>& mask) {
        std::size_t nx = phi[0].size();
        auto isInteriorCell = [&](std::size_t i, std::size_t j) -> bool {
            bool edgeIndex = (i == 0 || j == 0 || i == nx - 1 || j == dy.size() - 1);
            return !edgeIndex && (mask[j+1][i] == 1 && mask[j-1][i] == 1 && 
                                  mask[j][i+1] == 1 && mask[j][i-1] == 1);
        };

        double curr_y = y0; // Y coordinate of bottom edge of the cell row

        for(std::size_t j = 0; j < dy.size(); j++) {
            double curr_x = x0[j]; // X coord of left edge of cell column
            for(std::size_t i = 0; i < nx; i++) {
                if(mask[j][i] == 0) {
                    curr_x += dx[j];
                    continue;
                } 
                boxes.emplace_back(curr_x, curr_y + dy[j], curr_x + dx[j], curr_y, phi[j][i] * dx[j] * dy[j]);
                if(!isInteriorCell(i, j)) {
                    boundaryBoxIndices.push_back(boxes.size() - 1);
                }
                curr_x += dx[j];
            }
            curr_y += dy[j];
        }

        setMinMaxCoords();
    }
    
    void FreeCoordBoxGrid::setMinMaxCoords() {
        minX = std::numeric_limits<double>::max();
        minY = std::numeric_limits<double>::max();
        maxX = std::numeric_limits<double>::min();
        maxY = std::numeric_limits<double>::min();
        for (auto& b: boxes) {
            if(b.topLeftX < minX) {
                minX = b.topLeftX;
            }
            if(b.botRightY < minY) {
                minY = b.botRightY;
            }
            if(b.topLeftY > maxY) {
                maxY = b.topLeftY;
            }
            if(b.botRightX > maxX) {
                maxX = b.botRightX;
            }
        }
    }
}