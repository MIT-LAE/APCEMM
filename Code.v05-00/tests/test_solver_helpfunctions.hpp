#include <tuple>
#include <limits>
#include <Util/ForwardDecl.hpp>
#include <Util/PhysConstant.hpp>

#ifndef TEST_SOLVER_HELPFUNCTIONS_H
#define TEST_SOLVER_HELPFUNCTIONS_H
double l2norm(const Vector_2D& vec){
    int ny = vec.size();
    int nx = vec[0].size();
    double norm = 0;
    for(int i = 0; i < nx; i++){
        for(int j = 0; j < ny; j++){
            norm += vec[j][i] * vec[j][i];
        }
    }
    return std::sqrt(norm);
}
double l2err(const Vector_2D& exact, const Vector_2D& solution){
    int ny = exact.size();
    int nx = exact[0].size();
    double err = 0;
    for(int i = 0; i < nx; i++){
        for(int j = 0; j < ny; j++){
            err += std::pow(solution[j][i] - exact[j][i],2);
        }
    }
    return std::sqrt(err)/l2norm(exact);
}
//Can't use this to test SANDS because nontrivial to add required source terms
void setVecToPrescribedSoln(Vector_2D& v, double t, double dx, double dy){
    int ny = v.size();
    int nx = v[0].size();    
    for(int i = 0; i < nx; i++){
        for(int j = 0; j < ny; j++){
            double x = dx*i;
            double y = dy*j;
            v[j][i] = std::exp(-t) * std::sin(physConst::PI*x) * std::sin(physConst::PI*y);
        }
    }
}

void initTestVecAdvection(Vector_2D& v, double dx, double dy){
    int ny = v.size();
    int nx = v[0].size();    
    for(int i = 0; i < nx; i++){
        for(int j = 0; j < ny; j++){
            double x = dx*i;
            double y = dy*j;
            if(x < 0.25 && y < 0.25)
                v[j][i] = std::sin(physConst::PI*4*x) * std::sin(physConst::PI*4*y);
            else
                v[j][i] = 0;
        }
    }
}

void initTestVecShear(Vector_2D& v, double dx, double dy){
    int ny = v.size();
    int nx = v[0].size();    
    for(int i = 0; i < nx; i++){
        for(int j = 0; j < ny; j++){
            double x = dx*i;
            double y = dy*j;
            if(x > 0.25){
                v[j][i] = 0;
                continue;
            }
            v[j][i] = std::sin(physConst::PI * 4 * x);
        }
    }
}

std::tuple<double, double, double> vecMax2D(const Vector_2D& v, double dx, double dy){
    double max = std::numeric_limits<double>::min();
    double max_x = 0;
    double max_y = 0;
    int ny = v.size();
    int nx = v[0].size();    
    for(int i = 0; i < nx; i++){
        for(int j = 0; j < ny; j++){
            if(v[j][i] > max){
                max = v[j][i];
                max_x = dx*i;
                max_y = dy*i;
            }
        }
    }
    //std::cout << "max:" << max << ", max x: " << max_j << ", maxi: " << max_i << std::endl;
    return std::make_tuple(max, max_x, max_y);
}

std::tuple<double, double, double> vecMin2D(const Vector_2D& v, double dx, double dy){
    double min = std::numeric_limits<double>::max();
    double min_x = 0;
    double min_y = 0;
    int ny = v.size();
    int nx = v[0].size();    
    for(int i = 0; i < nx; i++){
        for(int j = 0; j < ny; j++){
            if(v[j][i] < min){
                min = v[j][i];
                min_x = dx*i;
                min_y = dy*i;
            }
        }
    }
    //std::cout << "max:" << max << ", max x: " << max_j << ", maxi: " << max_i << std::endl;
    return std::make_tuple(min, min_x, min_y);
}
#endif