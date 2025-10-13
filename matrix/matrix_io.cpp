#include "matrix_io.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>

bool read_matrix_from_file(const std::string& filename, std::vector<std::vector<double>>& A, int n) {
    std::ifstream file(filename);

    if (!file.is_open()) 
    {
        std::cerr << "error: can't open the file" << filename << std::endl;
        return false;
    }

    for (int i = 0; i < n; ++i) 
    {
        for (int j = 0; j < n; ++j) 
        {
            if (!(file >> A[i][j])) 
            {
                file.close();
                std::cout << "error: wrong input\n";
                return false; 
            }
        }
    }
    file.close();
    return true;
}

void initialize_matrix(std::vector<std::vector<double>>& A, int k, int n) {
    for (int i = 1; i <= n; ++i) 
    {
        for (int j = 1; j <= n; ++j) 
        {    
            A[i - 1][j - 1] = f(k, n, i, j);
        }
    }
}

void initialize_b(std::vector<std::vector<double>>& A, std::vector<double>& b, int n)
{
    for (int i = 0; i < n; i++) 
    {
        double sum_value = 0.0;
        for (int k = 0; (2 * k + 1) < n; k++) 
        {
            sum_value += A[i][2 * k + 1];
        }
        b[i] = sum_value;
    }
}

double f(int k, int n, int i, int j) 
{
    switch (k) 
    {
        case 1:
            return n - std::max(i, j) + 1;
        case 2:
            return std::max(i, j);
        case 3:
            return std::abs(i - j);
        case 4:
            return 1.0 / (i + j - 1);
        default:
            throw std::invalid_argument("wrong formula number");
    }
}
