#include "gaussian_elimination.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#define EPS 1e-12

bool gaussian_elimination_slow(const std::vector<std::vector<double>>& A, const std::vector<double>& b, std::vector<double>& x, const int n) 
{
    x.resize(n, 0.0);

    // Копирование данных
    std::vector<std::vector<double>> A_copy = A;
    std::vector<double> b_copy = b;

    // Вектор перестановок для столбцов (чтобы отслеживать swap)
    std::vector<int> perm(n);
    for (int j = 0; j < n; ++j) perm[j] = j;  //оригинальный порядок

    for (int i = 0; i < n; ++i) 
    {
        // Поиск главного элемента в строке
        int maxCol = i;
        double maxVal = std::fabs(A_copy[i][i]);
        for (int j = i + 1; j < n; ++j) 
        {
            double val = std::fabs(A_copy[i][j]);
            if (val > maxVal) 
            {
                maxVal = val;
                maxCol = j;
            }
        }
        
        if (std::fabs(A_copy[i][maxCol]) < EPS)
        {
            std::cerr << "bruh\n";
            return false;
        }


        // Перестановка столбцов
        for (int k = 0; k < n; ++k) 
        {  
            std::swap(A_copy[k][i], A_copy[k][maxCol]);
        }

        std::swap(perm[i], perm[maxCol]);  //отслеживаем перестановку столбцов

        double factor = 1 / A_copy[i][i];
        b_copy[i] *= factor;  
        for (int j = i; j < n; ++j) 
        {
            A_copy[i][j] *= factor;
        }  

        // Прямой ход метода Гаусса
        for (int k = i + 1; k < n; ++k) 
        {
            double factor = A_copy[k][i] / A_copy[i][i];
            b_copy[k] -= factor * b_copy[i];
            for (int j = i; j < n; ++j) 
            {
                A_copy[k][j] -= factor * A_copy[i][j];
            }
        }
    }

    // Обратный ход метода Гаусса
    for (int i = n - 1; i >= 0; --i) 
    {
        x[i] = b_copy[i];
        for (int k = i - 1; k >= 0; --k) 
        {
            b_copy[k] -= A_copy[k][i] * x[i];
        }
    }

    std::vector<double> x_original(n);
    for (int j = 0; j < n; ++j) 
    {
        x_original[perm[j]] = x[j];
    }
    x = x_original;

    return true;
}

bool gaussian_elimination_fast(const std::vector<std::vector<double>>& A, const std::vector<double>& b, std::vector<double>& x, const int n) 
{
    x.resize(n, 0.0);

    // Копирование данных
    std::vector<std::vector<double>> A_copy = A;
    std::vector<double> b_copy = b;

    for (int i = 0; i < n; ++i) 
    {
        // Поиск главного элемента в строке
        int maxCol = i;
        double maxVal = std::fabs(A_copy[i][1]);
        for (int j = 0; j < n; ++j)
        {
            double val = std::fabs(A_copy[i][j]);
            if (val > maxVal) 
            {
                maxVal = val;
                maxCol = j;
            }
        }
        
        if (std::fabs(A_copy[i][maxCol]) < EPS)
            return false;

        double factor = 1 / A_copy[i][maxCol];
        b_copy[i] *= factor;  
        for (int j = 0; j < n; ++j)
        {
            A_copy[i][j] *= factor;
        }  

        // Прямой ход метода Гаусса
        for (int k = 0; k < n; ++k) 
        {
	    if (k == i) continue;
            double factor = A_copy[k][maxCol] / A_copy[i][maxCol];
            b_copy[k] -= factor * b_copy[i];
            for (int j = 0; j < n; ++j)
            {
                A_copy[k][j] -= factor * A_copy[i][j];
            }
        }
    }

    for (int i = 0; i < n; ++i)
    {
	int count = 0;
	for (int j = 0; j < n; ++j) 
        {
            if (std::fabs(A_copy[i][j] - 1.0) < EPS)
	    {
		x[j] = b_copy[i];
		count += 1;
	    }	
        }
        if (count != 1)
	{
            std::cerr << "error: the process gone wrong\n";
	    return false;        
        }
    }

    return true;
}

