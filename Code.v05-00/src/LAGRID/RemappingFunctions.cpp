#include "LAGRID/RemappingFunctions.hpp"
namespace LAGRID {

    void twoDGridVariable::addBuffer(double bufLen_left, double bufLen_right, double bufLen_top, double bufLen_bot, double fillValue) {
        int nx_before = xCoords.size();
        // int ny_before = yCoords.size();

        //Generate coords of buffer areas
        int numRows_topBuffer = std::floor(bufLen_top / dy);
        Vector_1D topBufferCoords(numRows_topBuffer);
        std::generate(topBufferCoords.begin(), topBufferCoords.end(), [this, j = 1] () mutable { return yCoords[yCoords.size() - 1] + (j++) * dy;});
       
        int numRows_botBuffer = std::floor(bufLen_bot / dy);
        Vector_1D botBufferCoords(numRows_botBuffer);
        std::generate(botBufferCoords.rbegin(), botBufferCoords.rend(), [this, j = 1] () mutable { return yCoords[0] - (j++) * dy;});

        int numCols_leftBuffer = std::floor(bufLen_left / dx);
        Vector_1D leftBufferCoords(numCols_leftBuffer);
        std::generate(leftBufferCoords.rbegin(), leftBufferCoords.rend(), [this, i = 1] () mutable { return xCoords[0] - (i++) * dx;});
       
        int numCols_rightBuffer = std::floor(bufLen_right / dx);
        Vector_1D rightBufferCoords(numCols_rightBuffer);
        std::generate(rightBufferCoords.begin(), rightBufferCoords.end(), [this, i = 1] () mutable { return xCoords[xCoords.size() - 1] + (i++) * dx;});

        //Add generated coords to xCoords and yCoords
        yCoords.insert(yCoords.begin(), botBufferCoords.begin(), botBufferCoords.end());
        yCoords.insert(yCoords.end(), topBufferCoords.begin(), topBufferCoords.end());
        xCoords.insert(xCoords.begin(), leftBufferCoords.begin(), leftBufferCoords.end());
        xCoords.insert(xCoords.end(), rightBufferCoords.begin(), rightBufferCoords.end());

        //Expand phi 2d array
        for(int j = 0; j < numRows_topBuffer + numRows_botBuffer; j++) {
            phi.push_back(Vector_1D(nx_before, fillValue));
        }

        //rotate right by numRows_botBuffer to get those on the other side.
        std::rotate(phi.rbegin(), phi.rbegin() + numRows_botBuffer, phi.rend());
        for(std::size_t j = 0; j < yCoords.size(); j++) {
            phi[j].insert(phi[j].begin(), numCols_leftBuffer, fillValue);
            phi[j].insert(phi[j].end(), numCols_rightBuffer, fillValue);
        }

    }

    /*
    Calculating the exact value of the loss function (defined by the amount of rectangular grid space not filled)
    would take O(n_boxes) ~ O(nx*ny) time. This is horribly slow once the contrail gets large.

    Therefore, when nx and ny are small, we will perform the exact calculation
    because it will be more sensitive and still cheap, 
    and just perform an estimate after nx * ny exceeds some value.
        */
    double diffusionLossFunctionExact(const FreeCoordBoxGrid& boxGrid, const Remapping& remapping) {
        
        Vector_1D fillRatios(remapping.ny * remapping.nx, 0.0);
        double cellArea = remapping.dx * remapping.dy;
        for(auto& b: boxGrid.boxes) {
            int startGridIdx_x = std::floor((b.topLeftX - remapping.x0) / remapping.dx);
            int endGridIdx_x = std::floor((b.botRightX - remapping.x0) / remapping.dx);
            int startGridIdx_y = std::floor((b.botRightY - remapping.y0) / remapping.dy);
            int endGridIdx_y = std::floor((b.topLeftY - remapping.y0) / remapping.dy);

            for (int j = startGridIdx_y; j <= endGridIdx_y; j++) {
                for(int i = startGridIdx_x; i <= endGridIdx_x; i++) {
                    double area = coveredArea(remapping, b, i, j);
                    int flattened_idx = j * remapping.nx + i; 
                    fillRatios[flattened_idx] += area / cellArea;
                }
            }
        }
        //Count completely untouched cells as 1 (perfect) when evaluating loss function.
        //Assumes zero overlap (completely valid in APCEMM's case)
        return std::accumulate(fillRatios.begin(), fillRatios.end(), 0.0,
                                    [](double current_sum, double fillValue) -> double {
                                        return fillValue == 0.0
                                                ? current_sum
                                                : current_sum + (1 - fillValue);
                                    }
                                );
    }
    double diffusionLossFunctionBoundaryEstimate(const FreeCoordBoxGrid& boxGrid, const Remapping& remapping) {
        /*
            TODO: Write approximation algorithm that's more thought out, if any optimization is actually needed at all. 
            In its current form this is not a great solution. -MX
         */
        Vector_1D fillRatios(remapping.ny * remapping.nx, 0.0);
        double cellArea = remapping.dx * remapping.dy;
        for(int box_idx: boxGrid.boundaryBoxIndices) {
            const MassBox& b = boxGrid.boxes[box_idx];
            int startGridIdx_x = std::floor((b.topLeftX - remapping.x0) / remapping.dx);
            int endGridIdx_x = std::floor((b.botRightX - remapping.x0) / remapping.dx);
            // int startGridIdx_y = std::floor((b.topLeftY - remapping.y0) / remapping.dy);
            int endGridIdx_y = std::floor((b.botRightY - remapping.y0) / remapping.dy);

            for (int j = endGridIdx_y; j <= endGridIdx_y; j++) {
                for(int i = startGridIdx_x; i <= endGridIdx_x; i++) {
                    double area = coveredArea(remapping, b, i, j);
                    int flattened_idx = j * remapping.nx + i; 

                    //numerator takes fraction of mass in the box allocated cell at i, j 
                    //then, divide by grid cell area to convert to concentration.
                    fillRatios[flattened_idx] += area/ cellArea;
                }
            }
        }
        //Count completely untouched cells as 1 (perfect) when evaluating loss function.
        //Assumes zero overlap (completely valid in APCEMM's case)
        return std::accumulate(fillRatios.begin(), fillRatios.end(), 0.0,
                                    [](double current_sum, double fillValue) -> double {
                                            return fillValue == 0.0
                                                    ? current_sum
                                                    : current_sum + (1 - fillValue);
                                    }
                                );
    }


    FreeCoordBoxGrid rectToBoxGrid(double dy_old, const Vector_1D& dy_new, double dx_old, double x0_old, double y0_new, const Vector_2D& phi_old, const vector<vector<int>>& mask) {
        int ny = phi_old.size();
        int nx = phi_old[0].size();

        //Preserving the phi_old for now, can optimize this to directly overwrite phi_old later if performance is an issue here
        Vector_2D phi_new(ny, Vector_1D(nx));
        Vector_1D dx_new(ny);
        Vector_1D x0_new(ny);
        /*
        for(int j = 0; j < ny; j++) {
            dx_new[j] = dx_old * dy_new[j] / dy_old;
            x0_new[j] = x0_old + nx * (dx_old - dx_new[j]) / 2;
            double cellAreaRatio = dy_new[j] * dy_new[j] / (dy_old * dy_old);
            for(int i = 0; i < nx; i++) {
                phi_new[j][i] = phi_old[j][i] / cellAreaRatio;
            }
        }
        */
        // SDE 2024-08-12: Kludge - above is causing significant and cumulative noise
        // due to errors propagating from the met_.Update() calculation, even when
        // pressure velocity is zero
        for(int j = 0; j < ny; j++) {
            dx_new[j] = dx_old;
            x0_new[j] = x0_old;
            double cellAreaRatio = dy_new[j] * dy_new[j] / (dy_old * dy_old);
            for(int i = 0; i < nx; i++) {
                phi_new[j][i] = phi_old[j][i];
            }
        }
        return FreeCoordBoxGrid(dx_new, dy_new, phi_new, x0_new, y0_new, mask);
    }

    twoDGridVariable mapToStructuredGrid(const FreeCoordBoxGrid& boxGrid, const Remapping& remapping) {
        Vector_1D xCoords(remapping.nx);
        Vector_1D yCoords(remapping.ny);

        //Using some type-deduced local variables i, j that are default const (added with C++14 generic lambdas), so need the mutable to be able to modify them in the lambda.
        std::generate(xCoords.begin(), xCoords.end(), [&remapping, i = 0] () mutable { return remapping.x0 + remapping.dx * ( (i++) + 0.5);});
        std::generate(yCoords.begin(), yCoords.end(), [&remapping, j = 0] () mutable { return remapping.y0 + remapping.dy * ( (j++) + 0.5);});
        Vector_2D phi(remapping.ny, Vector_1D(remapping.nx, 0));
        double cellArea = remapping.dx * remapping.dy;

        for(auto& b: boxGrid.boxes) {
            //Need the std::min bounding to deal with edge cases of some boundary cell being included in the remapping.
            int startGridIdx_x = std::max(std::floor((b.topLeftX - remapping.x0) / remapping.dx), 0.0);
            int endGridIdx_x = std::min(std::floor((b.botRightX - remapping.x0) / remapping.dx), static_cast<double>(phi[0].size() - 1));
            int startGridIdx_y = std::max(std::floor((b.botRightY - remapping.y0) / remapping.dy), 0.0);
            int endGridIdx_y = std::min(std::floor((b.topLeftY - remapping.y0) / remapping.dy), static_cast<double>(phi.size() - 1));

            for (int j = startGridIdx_y; j <= endGridIdx_y; j++) {
                for(int i = startGridIdx_x; i <= endGridIdx_x; i++) {
                    double area = coveredArea(remapping, b, i, j);
                    phi[j][i] += (b.mass * area / b.area()) / cellArea;
                }
            }
        }
        return twoDGridVariable(std::move(phi), std::move(xCoords), std::move(yCoords));
    }

    // Just return the unused fraction
    // One day we will replace all this with linear algebra. It will be glorious
    twoDGridVariable getUnusedFraction(const FreeCoordBoxGrid& boxGrid, const Remapping& remapping) {
        Vector_1D xCoords(remapping.nx);
        Vector_1D yCoords(remapping.ny);

        //Using some type-deduced local variables i, j that are default const (added with C++14 generic lambdas), so need the mutable to be able to modify them in the lambda.
        std::generate(xCoords.begin(), xCoords.end(), [&remapping, i = 0] () mutable { return remapping.x0 + remapping.dx * ( (i++) + 0.5);});
        std::generate(yCoords.begin(), yCoords.end(), [&remapping, j = 0] () mutable { return remapping.y0 + remapping.dy * ( (j++) + 0.5);});
        double cellArea = remapping.dx * remapping.dy;

        // How much of each grid cell has not been written to?
        Vector_2D frac_unused(remapping.ny, Vector_1D(remapping.nx));
        for (int j = 0; j < remapping.ny; j++) {
            for (int i = 0; i < remapping.nx; i++) {
                frac_unused[j][i] = 1.0;
            }
        }

        for(auto& b: boxGrid.boxes) {
            //Need the std::min bounding to deal with edge cases of some boundary cell being included in the remapping.
            int startGridIdx_x = std::max(std::floor((b.topLeftX - remapping.x0) / remapping.dx), 0.0);
            int endGridIdx_x = std::min(std::floor((b.botRightX - remapping.x0) / remapping.dx), static_cast<double>(frac_unused[0].size() - 1));
            int startGridIdx_y = std::max(std::floor((b.botRightY - remapping.y0) / remapping.dy), 0.0);
            int endGridIdx_y = std::min(std::floor((b.topLeftY - remapping.y0) / remapping.dy), static_cast<double>(frac_unused.size() - 1));

            for (int j = startGridIdx_y; j <= endGridIdx_y; j++) {
                for(int i = startGridIdx_x; i <= endGridIdx_x; i++) {
                    double area = coveredArea(remapping, b, i, j);
                    frac_unused[j][i] -= area / cellArea;
                }
            }
        }

        return twoDGridVariable(std::move(frac_unused), std::move(xCoords), std::move(yCoords));
    }

    Vector_2D initVarToGrid( double mass, const Vector_1D& xEdges, const Vector_1D& yEdges,
                            std::function<double(double, double)> weightFunction, double logBinRatio ) 
    {
        /* 
        @param mass: Integral of quantity over 2D grid [original units * m2]
        @param xEdges: X-edges of the grid [m]
        @param yEdges: Y-edges of the grid [m]
        @param weightFunction: function that takes parameters (x, y) and outputs a weight. e.g. gaussian
        @param logBinRatio: optional parameter to specify the log of the bin size ratio if mapping aerosol quantities.
        @return 2D grid of PDF desired distribution that conserves # particles.
         */
        int nx = xEdges.size() - 1;
        int ny = yEdges.size() - 1;

        double newMass = 0;
        Vector_2D gridPDF(ny, Vector_1D(nx));
        for(int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                double yCenter = (yEdges[j+1] + yEdges[j]) / 2;
                double xCenter = (xEdges[i+1] + xEdges[i]) / 2;
                double cellArea = (yEdges[j+1] - yEdges[j]) * (xEdges[i+1] - xEdges[i]);
                double weight = weightFunction(xCenter, yCenter);
                gridPDF[j][i] = weight;
                newMass += weight * cellArea * logBinRatio;
            }
        }
        //We have no guarantees on the integral of the function, so need to scale to conserve mass
        double scalingFactor = mass / newMass;
        #pragma omp parallel for
        for(int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                gridPDF[j][i] *= scalingFactor;
            }
        }
        return gridPDF;
    }
    Vector_2D initVarToGridGaussian(double mass, const Vector_1D& xEdges, const Vector_1D& yEdges, double x0, double y0,
                            double sigmaX, double sigmaY, double logBinRatio ) 
    {
        auto gaussianFunc = [x0, y0, sigmaX, sigmaY] (double x, double y) -> double {
            return exp(- (pow(x - x0, 2.0) / (2.0 * sigmaX * sigmaX) + pow(y - y0, 2.0) / (2.0 * sigmaY * sigmaY)));
        };
        return initVarToGrid(mass, xEdges, yEdges, gaussianFunc, logBinRatio);
    }

    Vector_2D initVarToGridBimodalY(double mass, const Vector_1D& xEdges, const Vector_1D& yEdges, double x0, double y0,
                                    double width, double depth, double logBinRatio)
    {
        double sigmaX = width / 8;
        double omega = 2 * (physConst::PI / depth);
        auto func = [omega, x0, y0, sigmaX, depth] (double x, double y) -> double {
            return exp(- (pow(x - x0, 2.0) / (2.0 * sigmaX * sigmaX) )) * ( std::abs(y - y0) < depth/2 ) * std::abs( std::sin(omega * (y - y0)) );
        };
        return initVarToGrid(mass, xEdges, yEdges, func, logBinRatio);    
    }

}
