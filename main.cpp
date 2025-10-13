#include <iostream>
#include <chrono>
#include "matrix_io.h"
#include "gaussian_elimination.h"
#include "sum_func.h"

using namespace std;

int main(int argc, char* argv[]) 
{
    if (argc < 4 || argc > 5) 
    {
        cerr << "usage: " << argv[0] << " n m k filename\n";
        return -1;
    }

    int n = stoi(argv[1]); 
    int m = stoi(argv[2]);       
    int k = stoi(argv[3]);       

    if (m > n) 
    {
        cerr << "m must be less or equal than n\n";
        return -1;
    }

    bool ok = true;
    vector<vector<double>> A;
    vector<double> b, x_slow, x_fast;
    A.resize(n, vector<double>(n));
    b.resize(n);

    // Инициализация матрицы
    if (k == 0) 
    {
        string filename = argv[4];
        ok = read_matrix_from_file(filename, A, n);
        if (!ok) return -1;
    }
    else initialize_matrix(A, k, n);

    // Построение вектора b
    initialize_b(A, b, n);
    
    auto start = chrono::high_resolution_clock::now();
    ok = gaussian_elimination_slow(A, b, x_slow, n);
    auto end = chrono::high_resolution_clock::now();
    
    chrono::duration<double> slow_method_elapsed = end - start; 
 
    if (!ok)
    {
        cerr << "error: the matrix is singular!!!!!!\n";
        return -1;
    }

    start = chrono::high_resolution_clock::now();
    ok = gaussian_elimination_fast(A, b, x_fast, n);
    end = chrono::high_resolution_clock::now();

    chrono::duration<double> fast_method_elapsed = end - start;
   
    if (!ok) 
    {
        cerr << "error: the matrix is singular\n";
        return -1;
    }
    
    print_all(A, b, x_slow, x_fast, m, n, slow_method_elapsed, fast_method_elapsed);
    
    return 0;
}
