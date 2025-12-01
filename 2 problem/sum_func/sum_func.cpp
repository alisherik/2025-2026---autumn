#include "sum_func.h"

struct ThreadDataResidual {
    const std::vector<std::vector<double>>& A;
    const std::vector<double>& x;
    const std::vector<double>& b;
    std::vector<double>& Ax;
    double norm_b = 0.0;
    double norm_residual = 0.0;
    int n;
    int tid;
    int num_threads;

    ThreadDataResidual(const std::vector<std::vector<double>>& A_, const std::vector<double>& x_, const std::vector<double>& b_, std::vector<double>& Ax_, int n_, int tid_, int num_threads_)
        : A(A_), x(x_), b(b_), Ax(Ax_), n(n_), tid(tid_), num_threads(num_threads_) {}
};

void* calculateAx(void* arg) {
    ThreadDataResidual* data = static_cast<ThreadDataResidual*>(arg);
    int tid = data->tid;
    int num_threads = data->num_threads;
    for (int i = tid; i < data->n; i += num_threads) {
        for (int j = 0; j < data->n; ++j) {
            data->Ax[i] += data->A[i][j] * data->x[j];
        }
    }
    return nullptr;
}

void* calculateNorms(void* arg) {
    ThreadDataResidual* data = static_cast<ThreadDataResidual*>(arg);
    int tid = data->tid;
    int num_threads = data->num_threads;
    for (int i = tid; i < data->n; i += num_threads) {
        data->norm_b += data->b[i] * data->b[i];
        data->norm_residual += (data->Ax[i] - data->b[i]) * (data->Ax[i] - data->b[i]);
    }
    return nullptr;
}

double calculateResidualNorm(const std::vector<std::vector<double>>& A, const std::vector<double>& x, const std::vector<double>& b, const int n, const int num_threads) {
    std::vector<double> Ax(n, 0.0);
    std::vector<pthread_t> threads(num_threads);
    std::vector<ThreadDataResidual> thread_data;

    thread_data.reserve(num_threads);
    for (int tid = 0; tid < num_threads; ++tid) {
        thread_data.emplace_back(A, x, b, Ax, n, tid, num_threads);
        pthread_create(&threads[tid], nullptr, calculateAx, &thread_data[tid]);
    }
    for (int tid = 0; tid < num_threads; ++tid) {
        pthread_join(threads[tid], nullptr);
    }

    for (int tid = 0; tid < num_threads; ++tid) {
        pthread_create(&threads[tid], nullptr, calculateNorms, &thread_data[tid]);
    }
    for (int tid = 0; tid < num_threads; ++tid) {
        pthread_join(threads[tid], nullptr);
    }

    double norm_b = 0.0, norm_residual = 0.0;
    for (int tid = 0; tid < num_threads; ++tid) {
        norm_b += thread_data[tid].norm_b;
        norm_residual += thread_data[tid].norm_residual;
    }

    return std::sqrt(norm_residual) / std::sqrt(norm_b);
}

/////////////////////////////////////

struct ThreadDataError {
    const std::vector<double>& x;
    double norm = 0.0;
    int n;
    int tid;
    int num_threads;

    ThreadDataError(const std::vector<double>& x_, int n_, int tid_, int num_threads_)
        : x(x_), n(n_), tid(tid_), num_threads(num_threads_) {}
};

void* calculatePartialNormError(void* arg) {
    ThreadDataError* data = static_cast<ThreadDataError*>(arg);
    int tid = data->tid;
    int num_threads = data->num_threads;

    for (int k = tid; k < data->n; k += num_threads) {
        data->norm += (data->x[k] - (k % 2)) * (data->x[k] - (k % 2));
    }
    return nullptr;
}

double calculateNormError(const std::vector<double>& x, const int n, const int num_threads) {
    std::vector<pthread_t> threads(num_threads);
    std::vector<ThreadDataError> thread_data;
    thread_data.reserve(num_threads);

    for (int tid = 0; tid < num_threads; ++tid) {
        thread_data.emplace_back(x, n, tid, num_threads);
        pthread_create(&threads[tid], nullptr, calculatePartialNormError, &thread_data[tid]);
    }
    for (int tid = 0; tid < num_threads; ++tid) {
        pthread_join(threads[tid], nullptr);
    }

    double norm = 0.0;
    for (int tid = 0; tid < num_threads; ++tid) {
        norm += thread_data[tid].norm;
    }

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

void print_all(const std::vector<std::vector<double>>& A, const std::vector<double>& b, const std::vector<double>& x,const int m, const int n, const std::chrono::duration<double> elapsed, int p)
{
    std::cout << "number of threads: " << p << std::endl;
    std::cout << "matrix A:\n";
    print_matrix(A, m);
    std::cout << "vector b:\n";
    print_vector(b, m);
    std::cout << "solution x:\n";
    print_vector(x, m);
    
    // Вычисление и вывод нормы невязки
    std::cout << "residual norm: " << std::scientific << calculateResidualNorm(A, x, b, n, p) << std::endl;

    // Вычисление и вывод нормы погрешности
    std::cout << "error norm: " << std::scientific << calculateNormError(x, n, p) << std::endl;

    // Вывод времени выполнения
    std::cout << "decision time: " << elapsed.count() << " sec" << std::endl;
    
}
