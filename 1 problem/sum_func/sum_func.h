#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <vector>
#include <chrono>
#include <iostream>
#include <cmath>
#include <iomanip>

void Product(const std::vector<std::vector<double>>& A, const std::vector<double>& x, std::vector<double>& Ax, const int n);
void Difference(const std::vector<double>& a, const std::vector<double>& b, std::vector<double>& c, const int n);
double Norm(const std::vector<double>& a, const int n);
double calculate_residual_norm(const std::vector<std::vector<double>>& A, const std::vector<double>& x, const std::vector<double>& b, const int n);
double calculate_norm_error(const std::vector<double>& x, const int n);
void print_matrix(const std::vector<std::vector<double>>& A, int m);
void print_vector(const std::vector<double>& vec, int m);
void print_all(const std::vector<std::vector<double>>& A, const std::vector<double>& b, const std::vector<double>& x_slow, const std::vector<double>& x_fast, const int m, const int n, const std::chrono::duration<double> slow_method_elapsed, const std::chrono::duration<double> fast_method_elapsed);

#endif
