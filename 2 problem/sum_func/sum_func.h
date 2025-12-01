#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <vector>
#include <chrono>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <pthread.h>

void* calculateAx(void* arg);
void* calculateNorms(void* arg);
double calculateResidualNorm(const std::vector<std::vector<double>>& A, const std::vector<double>& x, const std::vector<double>& b, const int n, const int num_threads);
void* calculatePartialNormError(void* arg);
double calculateNormError(const std::vector<double>& x, const int n, const int num_threads);

void print_matrix(const std::vector<std::vector<double>>& A, int m);
void print_vector(const std::vector<double>& vec, int m);
void print_all(const std::vector<std::vector<double>>& A, const std::vector<double>& b, const std::vector<double>& x_slow, const int m, const int n, const std::chrono::duration<double> slow_method_elapsed, int p);

#endif
