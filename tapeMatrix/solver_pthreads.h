#define HAVE_STRUCT_TIMESPEC

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include "solver_serial.h"
#include "matrix.h"

namespace band_matrix_pthreads
{
	//typedef struct {
	//    int32_t tid;
	//    int32_t num_threads;
	//    int32_t n;
	//    int32_t b;
	//    struct DecomposeMatrix* result;
	//    pthread_barrier_t* barrier;
	//} LUContext;

	//typedef struct {
	//    int32_t tid;                      // ID потока
	//    int32_t num_threads;              // Количество потоков
	//    int32_t n;                        // Размер матрицы
	//    struct DecomposeMatrix* decomp;   // Разложение матрицы (L и U)
	//    struct Matrix* matrix;            // Исходная матрица (C и X)
	//    double* y;                        // Промежуточный вектор y
	//    double* local_sums;               // Локальные суммы для каждого потока
	//    pthread_barrier_t* barrier;       // Барьер для синхронизации
	//} SolveLUContext;

	//void* forward_substitution_worker(void* arg) {
	//    SolveLUContext* ctx = (SolveLUContext*)arg;
	//    int32_t n = ctx->n;
	//    int32_t num_threads = ctx->num_threads;
	//    int32_t tid = ctx->tid;

	//    for (int32_t i = 0; i < n; i++) {
	//        // Разделяем вычисление суммы s между потоками
	//        double local_s = 0.0;
	//        int32_t chunk = i / num_threads;
	//        int32_t remainder = i % num_threads;
	//        int32_t start = tid * chunk + (tid < remainder ? tid : remainder);
	//        int32_t count = chunk + (tid < remainder ? 1 : 0);
	//        int32_t end = start + count;

	//        // Вычисляем локальную сумму
	//        for (int32_t j = start; j < end; j++) {
	//            local_s += ctx->decomp->l[i][j] * ctx->y[j];
	//        }

	//        // Сохраняем локальную сумму
	//        ctx->local_sums[tid] = local_s;

	//        // Синхронизируем потоки перед редукцией
	//        pthread_barrier_wait(ctx->barrier);

	//        // Главный поток (tid = 0) выполняет редукцию
	//        if (tid == 0) {
	//            double global_s = 0.0;
	//            for (int32_t t = 0; t < num_threads; t++) {
	//                global_s += ctx->local_sums[t];
	//            }
	//            ctx->y[i] = ctx->matrix->C[i] - global_s;
	//        }

	//        // Синхронизируем потоки после редукции
	//        pthread_barrier_wait(ctx->barrier);
	//    }
	//    return NULL;
	//}

	//void* backward_substitution_worker(void* arg) {
	//    SolveLUContext* ctx = (SolveLUContext*)arg;
	//    int32_t n = ctx->n;
	//    int32_t num_threads = ctx->num_threads;
	//    int32_t tid = ctx->tid;

	//    for (int32_t i = n - 1; i >= 0; i--) {
	//        // Разделяем вычисление суммы s между потоками
	//        double local_s = 0.0;
	//        int32_t chunk = (n - i - 1) / num_threads;
	//        int32_t remainder = (n - i - 1) % num_threads;
	//        int32_t start = i + 1 + tid * chunk + (tid < remainder ? tid : remainder);
	//        int32_t count = chunk + (tid < remainder ? 1 : 0);
	//        int32_t end = start + count;

	//        // Вычисляем локальную сумму
	//        for (int32_t j = start; j < end; j++) {
	//            local_s += ctx->decomp->u[i][j] * ctx->matrix->X[j];
	//        }

	//        // Сохраняем локальную сумму
	//        ctx->local_sums[tid] = local_s;

	//        // Синхронизируем потоки перед редукцией
	//        pthread_barrier_wait(ctx->barrier);

	//        // Главный поток (tid = 0) выполняет редукцию
	//        if (tid == 0) {
	//            double global_s = 0.0;
	//            for (int32_t t = 0; t < num_threads; t++) {
	//                global_s += ctx->local_sums[t];
	//            }
	//            ctx->matrix->X[i] = (ctx->y[i] - global_s) / ctx->decomp->u[i][i];
	//        }

	//        // Синхронизируем потоки после редукции
	//        pthread_barrier_wait(ctx->barrier);
	//    }
	//    return NULL;
	//}

	//void* lu_decomposition_worker(void* arg) {
	//    LUContext* ctx = (LUContext*)arg;
	//    int32_t n = ctx->n;
	//    int32_t b = ctx->b;
	//    int32_t num_threads = ctx->num_threads;
	//    int32_t tid = ctx->tid;

	//    for (int32_t k = 0; k < n - 1; ++k) {
	//        int32_t upper_bound = (k + b + 1 < n) ? (k + b + 1) : n;
	//        int32_t total_rows = (upper_bound > k + 1) ? (upper_bound - (k + 1)) : 0;

	//        if (total_rows > 0) {
	//            int32_t base = k + 1;
	//            int32_t chunk = total_rows / num_threads;
	//            int32_t remainder = total_rows % num_threads;
	//            int32_t start = base + tid * chunk + (tid < remainder ? tid : remainder);
	//            int32_t count = chunk + (tid < remainder ? 1 : 0);
	//            int32_t end = start + count;

	//            for (int32_t i = start; i < end; ++i) {
	//                ctx->result->l[i][k] = ctx->result->u[i][k] / ctx->result->u[k][k];
	//                for (int32_t j = k; j < upper_bound; ++j) {
	//                    ctx->result->u[i][j] -= ctx->result->l[i][k] * ctx->result->u[k][j];
	//                }
	//            }
	//        }
	//        pthread_barrier_wait(ctx->barrier);
	//    }
	//    return NULL;
	//}

	//struct DecomposeMatrix lu_decomposition(struct Matrix matrix, int32_t num_threads) {
	//    struct DecomposeMatrix result;
	//    int32_t n = matrix.n;
	//    int32_t b = matrix.b;

	//    result.l = (double**)malloc(n * sizeof(double*));
	//    result.u = (double**)malloc(n * sizeof(double*));

	//    for (int32_t i = 0; i < n; i++) {
	//        result.l[i] = (double*)calloc(n, sizeof(double));
	//        result.u[i] = (double*)malloc(n * sizeof(double));
	//    }

	//    for (int32_t i = 0; i < n; i++) {
	//        for (int32_t j = 0; j < n; j++) {
	//            result.u[i][j] = matrix.A[i][j];
	//        }
	//        result.l[i][i] = 1.0;
	//    }

	//    if (n < 1000 || num_threads <= 1) {
	//        return band_matrix_serial::lu_decomposition(matrix);
	//    }

	//    pthread_t* threads = (pthread_t*)malloc(num_threads * sizeof(pthread_t));
	//    LUContext* contexts = (LUContext*)malloc(num_threads * sizeof(LUContext));
	//    pthread_barrier_t barrier;
	//    pthread_barrier_init(&barrier, NULL, num_threads);

	//    for (int32_t t = 0; t < num_threads; t++) {
	//        contexts[t].tid = t;
	//        contexts[t].num_threads = num_threads;
	//        contexts[t].n = n;
	//        contexts[t].b = b;
	//        contexts[t].result = &result;
	//        contexts[t].barrier = &barrier;
	//        pthread_create(&threads[t], NULL, lu_decomposition_worker, &contexts[t]);
	//    }

	//    for (int32_t t = 0; t < num_threads; t++) {
	//        pthread_join(threads[t], NULL);
	//    }
	//    pthread_barrier_destroy(&barrier);
	//    free(threads);
	//    free(contexts);

	//    return result;
	//}

	//void solve_lu(struct DecomposeMatrix decompose_matrix, struct Matrix* matrix, int32_t num_threads) 
	//{
	//    int32_t n = matrix->n;
	//    double* y = (double*)malloc(n * sizeof(double));
	//    matrix->X = (double*)malloc(n * sizeof(double));

	//    if (n < 1000 || num_threads <= 1) {
	//        // Для малых матриц или одного потока используем последовательную версию
	//        band_matrix_serial::solve_lu(decompose_matrix, matrix);
	//        free(y);
	//        return;
	//    }

	//    pthread_t* threads = (pthread_t*)malloc(num_threads * sizeof(pthread_t));
	//    SolveLUContext* contexts = (SolveLUContext*)malloc(num_threads * sizeof(SolveLUContext));
	//    pthread_barrier_t barrier;
	//    pthread_barrier_init(&barrier, NULL, num_threads);

	//    // Выделяем память для локальных сумм
	//    double* local_sums = (double*)malloc(num_threads * sizeof(double));

	//    // Инициализация контекстов
	//    for (int32_t t = 0; t < num_threads; t++) {
	//        contexts[t].tid = t;
	//        contexts[t].num_threads = num_threads;
	//        contexts[t].n = n;
	//        contexts[t].decomp = &decompose_matrix;
	//        contexts[t].matrix = matrix;
	//        contexts[t].y = y;
	//        contexts[t].local_sums = local_sums;
	//        contexts[t].barrier = &barrier;
	//    }

	//    // Прямой ход
	//    for (int32_t t = 0; t < num_threads; t++) {
	//        pthread_create(&threads[t], NULL, forward_substitution_worker, &contexts[t]);
	//    }
	//    for (int32_t t = 0; t < num_threads; t++) {
	//        pthread_join(threads[t], NULL);
	//    }

	//    // Обратный ход
	//    for (int32_t t = 0; t < num_threads; t++) {
	//        pthread_create(&threads[t], NULL, backward_substitution_worker, &contexts[t]);
	//    }
	//    for (int32_t t = 0; t < num_threads; t++) {
	//        pthread_join(threads[t], NULL);
	//    }

	//    // Освобождаем ресурсы
	//    pthread_barrier_destroy(&barrier);
	//    free(threads);
	//    free(contexts);
	//    free(local_sums);
	//    free(y);
	//}

	//void free_decompose_matrix(struct DecomposeMatrix* dm, int32_t n) {
	//    for (int32_t i = 0; i < n; i++) {
	//        free(dm->l[i]);
	//        free(dm->u[i]);
	//    }
	//    free(dm->l);
	//    free(dm->u);
	//}

	typedef struct {
		int32_t tid;                      // ID потока
		int32_t num_threads;              // Количество потоков
		int32_t n;                        // Размер матрицы
		int32_t b;                        // Полоса матрицы (для LU-разложения)
		struct DecomposeMatrix* decomp;   // Разложение матрицы (L и U)
		struct Matrix* matrix;            // Исходная матрица (A, C, X)
		double* y;                        // Промежуточный вектор y (для solve_lu)
		double* local_sums;               // Локальные суммы для каждого потока (для solve_lu)
		pthread_barrier_t* barrier;       // Барьер для синхронизации
	} ThreadContext;

	static void run_threads(pthread_t* threads, ThreadContext* contexts, int32_t num_threads,
		void* (*worker)(void*)) {
		for (int32_t t = 0; t < num_threads; t++) {
			pthread_create(&threads[t], NULL, worker, &contexts[t]);
		}
		for (int32_t t = 0; t < num_threads; t++) {
			pthread_join(threads[t], NULL);
		}
	}

	void* lu_decomposition_worker(void* arg) {
		ThreadContext* ctx = (ThreadContext*)arg;
		int32_t n = ctx->n;
		int32_t b = ctx->b;
		int32_t num_threads = ctx->num_threads;
		int32_t tid = ctx->tid;

		for (int32_t k = 0; k < n - 1; k++) {
			int32_t upper_bound = (k + b + 1 < n) ? (k + b + 1) : n;
			int32_t total_rows = upper_bound - (k + 1);
			if (total_rows <= 0) {
				pthread_barrier_wait(ctx->barrier);
				continue;
			}

			int32_t base = k + 1;
			int32_t chunk = total_rows / num_threads;
			int32_t remainder = total_rows % num_threads;
			int32_t start = base + tid * chunk + (tid < remainder ? tid : remainder);
			int32_t count = chunk + (tid < remainder ? 1 : 0);
			int32_t end = start + count;

			for (int32_t i = start; i < end; i++) {
				ctx->decomp->l[i][k] = ctx->decomp->u[i][k] / ctx->decomp->u[k][k];
				for (int32_t j = k; j < upper_bound; j++) {
					ctx->decomp->u[i][j] -= ctx->decomp->l[i][k] * ctx->decomp->u[k][j];
				}
			}
			pthread_barrier_wait(ctx->barrier);
		}
		return NULL;
	}

	void* forward_substitution_worker(void* arg) {
		ThreadContext* ctx = (ThreadContext*)arg;
		int32_t n = ctx->n;
		int32_t num_threads = ctx->num_threads;
		int32_t tid = ctx->tid;

		for (int32_t i = 0; i < n; i++) {
			double local_s = 0.0;
			int32_t chunk = i / num_threads;
			int32_t remainder = i % num_threads;
			int32_t start = tid * chunk + (tid < remainder ? tid : remainder);
			int32_t count = chunk + (tid < remainder ? 1 : 0);
			int32_t end = start + count;

			for (int32_t j = start; j < end; j++) {
				local_s += ctx->decomp->l[i][j] * ctx->y[j];
			}
			ctx->local_sums[tid] = local_s;
			pthread_barrier_wait(ctx->barrier);

			if (tid == 0) {
				double global_s = 0.0;
				for (int32_t t = 0; t < num_threads; t++) {
					global_s += ctx->local_sums[t];
				}
				ctx->y[i] = ctx->matrix->C[i] - global_s;
			}
			pthread_barrier_wait(ctx->barrier);
		}
		return NULL;
	}

	void* backward_substitution_worker(void* arg) {
		ThreadContext* ctx = (ThreadContext*)arg;
		int32_t n = ctx->n;
		int32_t num_threads = ctx->num_threads;
		int32_t tid = ctx->tid;

		for (int32_t i = n - 1; i >= 0; i--) {
			double local_s = 0.0;
			int32_t chunk = (n - i - 1) / num_threads;
			int32_t remainder = (n - i - 1) % num_threads;
			int32_t start = i + 1 + tid * chunk + (tid < remainder ? tid : remainder);
			int32_t count = chunk + (tid < remainder ? 1 : 0);
			int32_t end = start + count;

			for (int32_t j = start; j < end; j++) {
				local_s += ctx->decomp->u[i][j] * ctx->matrix->X[j];
			}
			ctx->local_sums[tid] = local_s;
			pthread_barrier_wait(ctx->barrier);

			if (tid == 0) {
				double global_s = 0.0;
				for (int32_t t = 0; t < num_threads; t++) {
					global_s += ctx->local_sums[t];
				}
				ctx->matrix->X[i] = (ctx->y[i] - global_s) / ctx->decomp->u[i][i];
			}
			pthread_barrier_wait(ctx->barrier);
		}
		return NULL;
	}

	struct DecomposeMatrix lu_decomposition(struct Matrix matrix, int32_t num_threads) {
		int32_t n = matrix.n;
		int32_t b = matrix.b;
		struct DecomposeMatrix result;

		if (n < 1000 || num_threads <= 1) {
			return band_matrix_serial::lu_decomposition(matrix);
		}

		// Выделение памяти для L и U
		result.l = (double**)malloc(n * sizeof(double*));
		result.u = (double**)malloc(n * sizeof(double*));
		for (int32_t i = 0; i < n; i++) {
			result.l[i] = (double*)calloc(n, sizeof(double));
			result.u[i] = (double*)malloc(n * sizeof(double));
			result.l[i][i] = 1.0;
			for (int32_t j = 0; j < n; j++) {
				result.u[i][j] = matrix.A[i][j];
			}
		}

		// Создание потоков
		pthread_t* threads = (pthread_t*)malloc(num_threads * sizeof(pthread_t));
		ThreadContext* contexts = (ThreadContext*)malloc(num_threads * sizeof(ThreadContext));
		pthread_barrier_t barrier;
		pthread_barrier_init(&barrier, NULL, num_threads);

		for (int32_t t = 0; t < num_threads; t++) {
			contexts[t].tid = t;
			contexts[t].num_threads = num_threads;
			contexts[t].n = n;
			contexts[t].b = b;
			contexts[t].decomp = &result;
			contexts[t].barrier = &barrier;
			contexts[t].matrix = NULL;  // Не используется в LU-разложении
			contexts[t].y = NULL;       // Не используется в LU-разложении
			contexts[t].local_sums = NULL;  // Не используется в LU-разложении
		}

		run_threads(threads, contexts, num_threads, lu_decomposition_worker);

		// Освобождение ресурсов
		pthread_barrier_destroy(&barrier);
		free(threads);
		free(contexts);

		return result;
	}

	void solve_lu(struct DecomposeMatrix decompose_matrix, struct Matrix* matrix, int32_t num_threads) {
		int32_t n = matrix->n;
		double* y = (double*)malloc(n * sizeof(double));
		matrix->X = (double*)malloc(n * sizeof(double));

		if (n < 1000 || num_threads <= 1) {
			band_matrix_serial::solve_lu(decompose_matrix, matrix);
			free(y);
			return;
		}

		// Создание потоков
		pthread_t* threads = (pthread_t*)malloc(num_threads * sizeof(pthread_t));
		ThreadContext* contexts = (ThreadContext*)malloc(num_threads * sizeof(ThreadContext));
		pthread_barrier_t barrier;
		pthread_barrier_init(&barrier, NULL, num_threads);

		// Выделение памяти для локальных сумм
		double* local_sums = (double*)malloc(num_threads * sizeof(double));

		// Инициализация контекстов
		for (int32_t t = 0; t < num_threads; t++) {
			contexts[t].tid = t;
			contexts[t].num_threads = num_threads;
			contexts[t].n = n;
			contexts[t].b = 0;  // Не используется в solve_lu
			contexts[t].decomp = &decompose_matrix;
			contexts[t].matrix = matrix;
			contexts[t].y = y;
			contexts[t].local_sums = local_sums;
			contexts[t].barrier = &barrier;
		}

		// Прямой ход
		run_threads(threads, contexts, num_threads, forward_substitution_worker);

		// Обратный ход
		run_threads(threads, contexts, num_threads, backward_substitution_worker);

		// Освобождение ресурсов
		pthread_barrier_destroy(&barrier);
		free(threads);
		free(contexts);
		free(local_sums);
		free(y);
	}

	void free_decompose_matrix(struct DecomposeMatrix* dm, int32_t n) {
		for (int32_t i = 0; i < n; i++) {
			free(dm->l[i]);
			free(dm->u[i]);
		}
		free(dm->l);
		free(dm->u);
	}
}