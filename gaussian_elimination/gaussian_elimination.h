#ifndef GAUSSIAN_ELIMINATION_H
#define GAUSSIAN_ELIMINATION_H

#include <vector>

bool gaussian_elimination_slow(const std::vector<std::vector<double>>& A, const std::vector<double>& b, std::vector<double>& x, const int n);
bool gaussian_elimination_fast(const std::vector<std::vector<double>>& A, const std::vector<double>& b, std::vector<double>& x, const int n);

#endif
