#ifndef GAUSSIAN_ELIMINATION_H
#define GAUSSIAN_ELIMINATION_H

#include <vector>
#include <pthread.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#define EPS 1e-12

void* gaussianStep(void* arg);
bool gaussian_elimination(const std::vector<std::vector<double>>& A, const std::vector<double>& b, std::vector<double>& x, const int n, const int num_threads);

#endif
