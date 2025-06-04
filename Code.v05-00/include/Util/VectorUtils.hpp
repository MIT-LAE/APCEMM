#ifndef VECTORUTILS_H
#define VECTORUTILS_H

#include <functional>
#include "Util/ForwardDecl.hpp"

namespace VectorUtils {
    using std::vector;

    Vector_2D cellAreas(const Vector_1D& xEdges, const Vector_1D& yEdges);

    double min(const Vector_2D& vec);

    double max(const Vector_2D& vec);
    Vector_1D max(const Vector_2D& vec, int axis);

    struct MaskInfo {
        int count;
        double maxX;
        double maxY;
        double minX;
        double minY;
    };

    std::pair<vector<vector<int>>, int> mask(const Vector_2D& vec, std::function<bool (double)> maskFunc);
    std::pair<vector<vector<int>>, MaskInfo> mask(const Vector_2D& vec, const Vector_1D& xEdges, const Vector_1D& yEdges, std::function<bool (double)> maskFunc);
}

#endif
