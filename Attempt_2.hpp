#ifndef ATTEMPT_2_H
#define ATTEMPT_2_H

#include <vector>

double ffunc(double t, double &p, double &q, double &r);
void derivative_function(
    std::vector<double> x, std::vector<std::vector<double>> dcm, 
    std::vector<double> v_body, double time, 
    std::vector<double> &xdot, double &Cdot
);
double dotProduct(std::vector<std::vector<double>> a, std::vector<std::vector<double>> b);
std::vector<double> addVectors(std::vector<double> a, std::vector<double> b);
std::vector<double> scaleVector(std::vector<double> vec, double scalar);
std::vector<std::vector<double>> scaleMatrix(std::vector<std::vector<double>> vec, double scalar);
std::vector<std::vector<double>> addMatrix(
    std::vector<std::vector<double>> a, std::vector<std::vector<double>> b
);
std::vector<double> add4Vectors(
    std::vector<double> a, std::vector<double> b, 
    std::vector<double> c, std::vector<double> d
);
std::vector<std::vector<double>> add4Matrix(
    std::vector<std::vector<double>> a, std::vector<std::vector<double>> b, 
    std::vector<std::vector<double>> c, std::vector<std::vector<double>> d
);

#endif // ATTEMPT_2_H
