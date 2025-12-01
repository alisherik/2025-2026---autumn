#include "matrix_io.h"

bool read_matrix_from_file(const std::string& filename, std::vector<std::vector<double>>& A, int n) 
{
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


struct ThreadData 
{
    std::vector<std::vector<double>>* A;
    int k, n;
    int start_row, end_row; 
};



void* initializeMatrixThread(void* arg) 
{
    ThreadData* data = static_cast<ThreadData*>(arg);

    for (int i = data->start_row; i < data->end_row; ++i) 
    {
        for (int j = 0; j < data->n; ++j) 
	{
            (*data->A)[i][j] = f(data->k, data->n, i + 1, j + 1);
        }
    }
    return nullptr;
}

void initialize_matrix(std::vector<std::vector<double>>& A, int k, int n, int num_threads) 
{
    std::vector<pthread_t> threads(num_threads);
    std::vector<ThreadData> thread_data(num_threads);

    int rows_per_thread = n / num_threads;
    int remaining_rows = n % num_threads;
    int start_row = 0;

    
    for (int tid = 0; tid < num_threads; ++tid) 
    {
        int end_row = start_row + rows_per_thread + (tid < remaining_rows ? 1 : 0);
        thread_data[tid] = {&A, k, n, start_row, end_row};  //                         A
        pthread_create(&threads[tid], nullptr, initializeMatrixThread, &thread_data[tid]);
        start_row = end_row;
    }

    
    for (int tid = 0; tid < num_threads; ++tid) 
    {
        pthread_join(threads[tid], nullptr);
    }
}
