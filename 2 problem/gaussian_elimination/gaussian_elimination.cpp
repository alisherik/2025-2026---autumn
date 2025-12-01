#include "gaussian_elimination.h"

pthread_barrier_t barrier; // Барьер для синхронизации
pthread_mutex_t mutex; // Мьютекс для синхронизации доступа к глобальному максимуму

struct ThreadData 
{
    std::vector<std::vector<double>>& A;
    std::vector<double>& b;
    std::vector<double>& x;
    std::vector<double>& x_temp;
    std::vector<int>& perm;
    int n;
    int tid; // ID потока
    int num_threads;
    int* globalMaxCol; // Указатель на переменную для строки с максимальным элементом

    ThreadData(std::vector<std::vector<double>>& A_, std::vector<double>& b_, std::vector<double>& x_, std::vector<double>& x_temp_, std::vector<int>& perm_, int n_, int tid_, int num_threads_, int* globalMaxCol_)
        : A(A_), b(b_), x(x_), x_temp(x_temp_), perm(perm_), n(n_), tid(tid_), num_threads(num_threads_), globalMaxCol(globalMaxCol_) {}
};

void* gaussianStep(void* arg) 
{
    struct ThreadData *data = (struct ThreadData*)arg;
    int n = data->n;
    int tid = data->tid;
    int num_threads = data->num_threads;
    
    //std::cout << "Thread " << tid << " started" << std::endl;
    for (int i = 0; i < n; ++i) 
    {
        //std::cout << "Thread " << tid << " at iteration " << i << std::endl;
        if (tid == 0) *(data->globalMaxCol) = i;
        
        // 1. Поиск строки с максимальным элементом в строке i
        int localMaxCol = i;
        for (int k = i + tid; k < n; k += num_threads) 
	{
            if (std::fabs(data->A[i][k]) > std::fabs(data->A[i][localMaxCol])) 
	    {
                localMaxCol = k;
            }
        }
        // Синхронизация всех потоков
        pthread_barrier_wait(&barrier);

        // Обновляем глобальный максимум, если локальный больше
        pthread_mutex_lock(&mutex);
        if (std::fabs(data->A[i][localMaxCol]) > std::fabs(data->A[i][*data->globalMaxCol])) 
        {
            *(data->globalMaxCol) = localMaxCol;
        }
        pthread_mutex_unlock(&mutex);

        // Синхронизация всех потоков после поиска максимума
        pthread_barrier_wait(&barrier);

	// Проверка на невырожденность матрицы
	if (tid == 0) 
	{
		if (std::fabs(data->A[i][*data->globalMaxCol]) < EPS) 
	        {
			std::cout << "The matrix is singular\n";
			return (void*)1;
            	}
	}

	// Перестановка столбцов
        if (tid == 0) 
	{
            for (int k = 0; k < n; ++k) 
            {  
		std::swap(data->A[k][i], data->A[k][*data->globalMaxCol]);
            }
            std::swap(data->perm[i], data->perm[*data->globalMaxCol]);  //отслеживаем перестановку столбцов
            double factor = 1 / data->A[i][i];
            data->b[i] *= factor;
            for (int j = i; j < n; ++j) 
	    {
                data->A[i][j] *= factor;
	    }            
        }

        // Синхронизация перед началом прямого хода
        pthread_barrier_wait(&barrier);

        // 2. Прямой ход
        for (int k = i + 1 + tid; k < n; k += num_threads) 
	{
            double factor = data->A[k][i];
            data->b[k] -= factor * data->b[i];
            for (int j = i; j < n; ++j) 
	    {
                data->A[k][j] -= factor * data->A[i][j];
            }
        }

        // Синхронизация после завершения прямого хода
        pthread_barrier_wait(&barrier);
        
        //std::cout << "Thread " << tid << " passed barrier after iteration " << i << std::endl;
    }

    // 3. Обратный ход
    for (int i = n - 1; i >= 0; --i)
    {
        if (tid == 0) data->x_temp[i] = data->b[i];

        // Синхронизация после вычисления x[i]
        pthread_barrier_wait(&barrier);

        for (int k = i - 1 - tid; k >= 0; k -= num_threads) 
	{
            data->b[k] -= data->A[k][i] * data->x_temp[i];
        }

        // Синхронизация перед переходом к следующей строке
        pthread_barrier_wait(&barrier);
    }
	
    if (tid == 0) 
    {
	    for (int j = 0; j < n; ++j) 
	    {
	        data->x[data->perm[j]] = data->x_temp[j];
	    }
    }
    
    //std::cout << "Thread " << tid << " finished" << std::endl;
    return (void*)0;
}

bool gaussian_elimination(const std::vector<std::vector<double>>& A, const std::vector<double>& b, std::vector<double>& x, const int n, const int num_threads)
{
    x.resize(n, 0.0);

    // Копирование данных
    std::vector<std::vector<double>> A_copy = A;
    std::vector<double> b_copy = b;

    // Инициализация барьера для синхронизации потоков
    pthread_barrier_init(&barrier, nullptr, num_threads);
    pthread_mutex_init(&mutex, nullptr);

    int globalMaxCol = 0;

    std::vector<double> x_temp(n, 0.0);

    // Вектор перестановок
    std::vector<int> perm(n);
    for (int j = 0; j < n; ++j) perm[j] = j;  //оригинальный порядок

    std::vector<pthread_t> threads(num_threads);
    std::vector<ThreadData*> thread_data(num_threads);

    for (int tid = 0; tid < num_threads; ++tid) 
    {
        thread_data[tid] = new ThreadData(A_copy, b_copy, x, x_temp, perm, n, tid, num_threads, &globalMaxCol);
        pthread_create(&threads[tid], nullptr, gaussianStep, thread_data[tid]);
    }

    bool success = true;
    // Ожидание завершения работы всех потоков
    for (int tid = 0; tid < num_threads; ++tid) 
    {
        void* retval;
        pthread_join(threads[tid], &retval);
	if (retval != (void*)0)
	{
	    success = false;
	}
        delete thread_data[tid];
    }
    
    // Уничтожение барьера и мьютекса
    pthread_barrier_destroy(&barrier);
    pthread_mutex_destroy(&mutex);

    return success;
}