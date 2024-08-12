#include "Util/VectorUtils.hpp"
namespace VectorUtils {
    Vector_2D cellAreas (const Vector_1D& xEdges, const Vector_1D& yEdges) {
        int nx = xEdges.size() - 1;
        int ny = yEdges.size() - 1;

        Vector_2D areas(ny, Vector_1D(nx));
        #pragma omp parallel for default(shared)
        for(int j = 0; j < ny; j++) {
            for(int i = 0; i < nx; i++) {
                areas[j][i] = (yEdges[j+1] - yEdges[j]) * (xEdges[i+1] - xEdges[i]);
            }
        }
        return areas;
    }

    double VecMin2D (const Vector_2D& vec) {
        double min = std::numeric_limits<double>::max();
        for (std::size_t j = 0; j < vec.size(); j++) {
            for (std::size_t i = 0; i < vec[0].size(); i++) {
                if(std::isinf(vec[j][i]) || std::isnan(vec[j][i])) throw std::out_of_range("in VectorUtils::VecMin2D: inf or nan found in vector");
                if(vec[j][i] < min) min = vec[j][i];
            }
        }
        return min;
    }

    double VecMax2D (const Vector_2D& vec) {
        double max = std::numeric_limits<double>::min();
        for (std::size_t j = 0; j < vec.size(); j++) {
            for (std::size_t i = 0; i < vec[0].size(); i++) {
                if(std::isinf(vec[j][i]) || std::isnan(vec[j][i])) throw std::out_of_range("in VectorUtils::VecMax2D: inf or nan found in vector");
                if(vec[j][i] > max) max = vec[j][i];
            }
        }
        return max;
    }

    Vector_1D VecMax2D (const Vector_2D& vec, int axis) {
        Vector_1D axis_max;

        if(axis == 0) {
            axis_max = Vector_1D(vec.size());
            for (std::size_t j = 0; j < vec.size(); j++) {
                double rowMax = std::numeric_limits<double>::min();
                for (std::size_t i = 0; i < vec[0].size(); i++) {
                    if(std::isinf(vec[j][i]) || std::isnan(vec[j][i])) throw std::out_of_range("in VectorUtils::VecMax2D: inf or nan found in vector");
                    if(vec[j][i] > rowMax) rowMax = vec[j][i];
                }
                axis_max[j] = rowMax;
            }
        }
        else {
            axis_max = Vector_1D(vec[0].size());
            for (std::size_t i = 0; i < vec[0].size(); i++) {
                double colMax = std::numeric_limits<double>::min();
                for (std::size_t j = 0; j < vec.size(); j++) {
                    if(std::isinf(vec[j][i]) || std::isnan(vec[j][i])) throw std::out_of_range("in VectorUtils::VecMax2D: inf or nan found in vector");
                    if(vec[j][i] > colMax) colMax = vec[j][i];
                }
                axis_max[i] = colMax;
            }
        }
        return axis_max;

    }

    double Vec2DSum (const Vector_2D& vec) {
        double sum = 0;
        for (std::size_t j = 0; j < vec.size(); j++) {
            for (std::size_t i = 0; i < vec[0].size(); i++) {
                sum += vec[j][i];
            }
        }
        return sum;
    }

    std::pair<vector<vector<int>>, int> Vec2DMask (const Vector_2D& vec, std::function<bool (double)> maskFunc) {
        vector<vector<int>> mask(vec.size(), vector<int>(vec[0].size()));
        int nonMaskedElems = 0;
        for (std::size_t j = 0; j < vec.size(); j++) {
            for (std::size_t i = 0; i < vec[0].size(); i++) {
                mask[j][i] = maskFunc(vec[j][i]);
                if(mask[j][i]) {
                    nonMaskedElems++;
                }
            }
        }

        return std::make_pair(std::move(mask), nonMaskedElems);
    }

    std::pair<vector<vector<int>>, MaskInfo> Vec2DMask (const Vector_2D& vec, const Vector_1D& xEdges, const Vector_1D& yEdges, std::function<bool (double)> maskFunc) {
        vector<vector<int>> mask(vec.size(), vector<int>(vec[0].size()));
        MaskInfo info;
        int nonMaskedElems = 0;
        double minX = std::numeric_limits<double>::max();
        double minY = std::numeric_limits<double>::max();
        double maxX = std::numeric_limits<double>::lowest();
        double maxY = std::numeric_limits<double>::lowest();

        for (std::size_t j = 0; j < vec.size(); j++) {
            for (std::size_t i = 0; i < vec[0].size(); i++) {
                mask[j][i] = maskFunc(vec[j][i]);
                if(mask[j][i]){
                    nonMaskedElems++;
                    if(xEdges[i] < minX) 
                        minX = xEdges[i];
                    else if(xEdges[i + 1] > maxX) 
                        maxX = xEdges[i + 1];

                    if(yEdges[j] < minY) 
                        minY = yEdges[j];
                    else if(yEdges[j + 1] > maxY) 
                        maxY = yEdges[j + 1];
                }
            }
        }
        info.count = nonMaskedElems;
        info.maxX = maxX;
        info.minX = minX;
        info.maxY = maxY;
        info.minY = minY;
        return std::make_pair(std::move(mask), std::move(info));
    }

    Vector_2D vec2DOperation(const Vector_2D& vec1, const Vector_2D& vec2, std::function<double (double, double)> transformFunc) {
        Vector_2D result(vec1.size(), Vector_1D(vec1[0].size()));

        #pragma omp parallel for
        for(std::size_t j = 0; j < vec1.size(); j++) {
            for(std::size_t i = 0; i < vec1[0].size(); i++) {
                result[j][i] = transformFunc(vec1[j][i], vec2[j][i]);
            }
        }
        return result;
    }


}