#ifndef MATRIX_IO_H
#define MATRIX_IO_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <pthread.h>

bool read_matrix_from_file(const std::string& filename, std::vector<std::vector<double>>& A, int n);
void initialize_b(std::vector<std::vector<double>>& A, std::vector<double>& b, int n);
double f(int k, int n, int i, int j);
void* initializeMatrixThread(void* arg);
void initialize_matrix(std::vector<std::vector<double>>& A, int k, int n, int num_threads);

#endif
