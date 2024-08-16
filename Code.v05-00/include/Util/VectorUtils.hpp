#ifndef VECTORUTILS_H
#define VECTORUTILS_H

#include <functional>
#include "Util/ForwardDecl.hpp"

namespace VectorUtils {
    using std::vector;
    Vector_2D cellAreas (const Vector_1D& xEdges, const Vector_1D& yEdges);
    
    double VecMin2D (const Vector_2D& vec);

    double VecMax2D (const Vector_2D& vec);
    Vector_1D VecMax2D (const Vector_2D& vec, int axis);
    double Vec2DSum (const Vector_2D& vec);
    
    struct MaskInfo {
        int count;
        double maxX;
        double maxY;
        double minX;
        double minY;
    };

    std::pair<vector<vector<int>>, int> Vec2DMask (const Vector_2D& vec, std::function<bool (double)> maskFunc);
    std::pair<vector<vector<int>>, MaskInfo> Vec2DMask (const Vector_2D& vec, const Vector_1D& xEdges, const Vector_1D& yEdges, std::function<bool (double)> maskFunc);

    template<typename FillWithType>
    void fill2DVec(Vector_2D& toFill, const FillWithType& fillWith, std::function<bool (double)> fillCondFunction) {
        #pragma omp parallel for
        for(std::size_t j = 0; j < toFill.size(); j++) {
            for(std::size_t i = 0; i < toFill[0].size(); i++) {
                double fillWithValue;
                if constexpr(std::is_arithmetic_v<FillWithType>) {
                    fillWithValue = fillWith;
                }
                else {
                    fillWithValue = fillWith[j][i];
                }
                toFill[j][i] = fillCondFunction(toFill[j][i]) ? fillWithValue : toFill[j][i];
            }
        }
    }
    Vector_2D vec2DOperation(const Vector_2D& vec1, const Vector_2D& vec2, std::function<double (double, double)> transformFunc);
}

#endif