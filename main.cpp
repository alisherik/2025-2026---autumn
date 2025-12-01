#include <iostream>
#include <chrono>
#include "matrix_io.h"
#include "gaussian_elimination.h"
#include "sum_func.h"

using namespace std;

int main(int argc, char* argv[]) 
{
    if (argc < 5 || argc > 6)
    {
        cerr << "usage: " << argv[0] << " n p m k filename\n";
        return -1;
    }

    int n = std::stoi(argv[1]);       // Размерность матрицы
    int p = std::stoi(argv[2]);       // Количество потоков
    int m = std::stoi(argv[3]);       // Количество выводимых значений
    int k = std::stoi(argv[4]);       // Номер формулы или 0 для чтения из файла

    if (m > n) 
    {
        cerr << "m must be less or equal than n\n";
        return -1;
    }

    if (p <= 0) {
        std::cerr << "p must be natural number" << std::endl;
        return -1;
    }

    bool ok = true;
    vector<vector<double>> A;
    vector<double> b, x;
    A.resize(n, vector<double>(n));
    b.resize(n);

    // Инициализация матрицы
    if (k == 0) 
    {
        string filename = argv[5];
        ok = read_matrix_from_file(filename, A, n);
        if (!ok) return -1;
    }
    else initialize_matrix(A, k, n, p);

    // Построение вектора b
    initialize_b(A, b, n);
    
    auto start = chrono::high_resolution_clock::now();
    ok = gaussian_elimination(A, b, x, n, p);
    auto end = chrono::high_resolution_clock::now();
    
    chrono::duration<double> elapsed = end - start; 
 
    if (!ok) return -1;
    
    print_all(A, b, x, m, n, elapsed, p);
    
    return 0;
}
