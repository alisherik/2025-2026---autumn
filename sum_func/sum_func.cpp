#include "sum_func.h"

void Product(const std::vector<std::vector<double>>& A, const std::vector<double>& x, std::vector<double>& Ax, const int n) 
{
    Ax.resize(n, 0.0);
    for (int i = 0; i < n; ++i) 
    {
        for (int j = 0; j < n; ++j) 
        {
            Ax[i] += A[i][j] * x[j];
        }
    }
}

void Difference(const std::vector<double>& a, const std::vector<double>& b, std::vector<double>& c, const int n) 
{
    c.resize(n);
    for (int i = 0; i < n; ++i) 
    {
        c[i] = a[i] - b[i];
    }
}

double Norm(const std::vector<double>& a, const int n) 
{
    double norm = 0;
    for (int i = 0; i < n; ++i) 
    {
        norm += a[i] * a[i];
    }
    return norm;
}

double calculate_residual_norm(const std::vector<std::vector<double>>& A, const std::vector<double>& x, const std::vector<double>& b, const int n) 
{
    std::vector<double> Ax;
    std::vector<double> c;
    double norm_b, norm_diff;
    Product(A, x, Ax, n);
    Difference(Ax, b, c, n);
    norm_diff = Norm(c, n);
    norm_b = Norm(b, n);
    return std::sqrt(norm_diff) / std::sqrt(norm_b);
}

double calculate_norm_error(const std::vector<double>& x, const int n) 
{
    double norm = 0.0;
    std::vector<double> c(n, 0.0);
    std::vector<double> d(n, 0.0);
    for (int k = 0; k < n; k++) 
    {
        c[k] = k % 2;
    }
    Difference(x, c, d, n);
    norm = Norm(d, n);
    return std::sqrt(norm);
}

void print_matrix(const std::vector<std::vector<double>>& A, int m) 
{
    for (int i = 0; i < m; ++i) 
    {
        for (int j = 0; j < m; ++j) 
        {
            std::cout << std::setw(10) << std::setprecision(3) << std::scientific << A[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

void print_vector(const std::vector<double>& vec, int m) 
{
    for (int i = 0; i < m; ++i) 
    {
        std::cout << std::setw(10) << std::setprecision(3) << std::scientific << vec[i] << " ";
    }
    std::cout << std::endl;
}

void print_all(const std::vector<std::vector<double>>& A, const std::vector<double>& b, const std::vector<double>& x_slow, const std::vector<double>& x_fast, const int m, const int n, const std::chrono::duration<double> slow_method_elapsed, const std::chrono::duration<double> fast_method_elapsed)
{
    std::cout << "matrix A:\n";
    print_matrix(A, m);
    std::cout << "vector b:\n";
    print_vector(b, m);
    std::cout << "solution x (slow version):\n";
    print_vector(x_slow, m);
    std::cout << "solution x (fast version):\n";
    print_vector(x_fast, m);   
    std::vector<double> Ax_slow(n, 0.0);
    Product(A, x_slow, Ax_slow, n);
    std::cout << "matrix Ax (slow version):\n";
    print_vector(Ax_slow, m);
    std::vector<double> Ax_fast(n, 0.0);
    Product(A, x_fast, Ax_fast, n);
    std::cout << "matrix Ax (fast version):\n";
    print_vector(Ax_fast, m);

    // Вычисление и вывод нормы невязки
    std::cout << "residual norm (slow version): " << std::scientific << calculate_residual_norm(A, x_slow, b, n) << std::endl;

    // Вычисление и вывод нормы погрешности
    std::cout << "error norm (slow version): " << std::scientific << calculate_norm_error(x_slow, n) << std::endl;

    // Вывод времени выполнения
    std::cout << "decision time (slow version): " << slow_method_elapsed.count() << " sec" << std::endl;
    std::cout << "decision time (fast version): " << fast_method_elapsed.count() << " sec" << std::endl;
    std::cout << "(decision time slow)/(decision time fast) = " << slow_method_elapsed.count()/fast_method_elapsed.count() << std::endl;
    
}
