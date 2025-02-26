#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "matrix.h"

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif

namespace band_matrix_omp
{
    void reverse_array(double* array, int32_t n)
    {
        double temp;
#pragma omp parallel for private(temp) schedule(static)
        for (int32_t i = 0; i < n / 2; ++i)
        {
            temp = array[i];
            array[i] = array[n - 1 - i];
            array[n - 1 - i] = temp;
        }
    }

    struct DecomposeMatrix lu_decomposition(struct Matrix matrix)
    {
        struct DecomposeMatrix result;
        int32_t n = matrix.n;
        int32_t b = matrix.b;

        result.l = (double**)malloc(n * sizeof(double*));
        result.u = (double**)malloc(n * sizeof(double*));

        // Параллельное выделение памяти
#pragma omp parallel for schedule(static)
        for (int32_t i = 0; i < n; i++)
        {
            result.l[i] = (double*)calloc(n, sizeof(double));
            result.u[i] = (double*)malloc(n * sizeof(double));
        }

        // Копирование матрицы в U
#pragma omp parallel for collapse(2) schedule(static)
        for (int32_t i = 0; i < n; i++)
        {
            for (int32_t j = 0; j < n; j++)
            {
                result.u[i][j] = matrix.A[i][j];
            }
        }

        // Инициализация диагонали L
#pragma omp parallel for schedule(static)
        for (int32_t i = 0; i < n; i++) {
            result.l[i][i] = 1.0;
        }

        // Основной цикл разложения
        for (int32_t k = 0; k < n - 1; ++k)
        {
            int32_t upper_bound = (k + b + 1 < n) ? (k + b + 1) : n;
            // Параллелизация внутреннего цикла
#pragma omp parallel for schedule(dynamic)
            for (int32_t i = k + 1; i < upper_bound; ++i)
            {
                result.l[i][k] = result.u[i][k] / result.u[k][k];
                for (int32_t j = k; j < upper_bound; ++j)
                {
                    result.u[i][j] -= result.l[i][k] * result.u[k][j];
                }
            }
        }
        return result;
    }

  /*  struct DecomposeMatrix lu_decomposition(struct Matrix matrix)
    {
        struct DecomposeMatrix result;

        int32_t n = matrix.n;
        int32_t b = matrix.b;

        result.l = (double**)malloc(n * sizeof(double*));
        result.u = (double**)malloc(n * sizeof(double*));

        #pragma omp parallel for
        for (int32_t i = 0; i < n; i++)
        {
            result.l[i] = (double*)calloc(n, sizeof(double));
            result.u[i] = (double*)malloc(n * sizeof(double));
        }

        #pragma omp parallel for collapse(2)
        for (int32_t i = 0; i < n; i++)
        {
            for (int32_t j = 0; j < n; j++)
            {
                result.u[i][j] = matrix.A[i][j];
            }
        }

        #pragma omp parallel for
        for (int32_t i = 0; i < n; i++) {
            result.l[i][i] = 1.0;
        }

        for (int32_t k = 0; k < n - 1; ++k)
        {
            int32_t upper_bound = (k + b + 1 < n) ? (k + b + 1) : n;
            #pragma omp parallel for shared(result, upper_bound, k)
            for (int32_t i = k + 1; i < upper_bound; ++i)
            {
                result.l[i][k] = result.u[i][k] / result.u[k][k];
                for (int32_t j = k; j < upper_bound; ++j)
                {
                    result.u[i][j] -= result.l[i][k] * result.u[k][j];
                }
            }
        }
        return result;
    }*/

    void solve_lu(DecomposeMatrix decomp, struct Matrix* matrix)
    {
        int32_t n = matrix->n;
        double* y = (double*)malloc(n * sizeof(double));

        // Решение Ly = b
        for (int32_t i = 0; i < n; i++)
        {
            double s = 0.0;
#pragma omp parallel for reduction(+:s) schedule(static)
            for (int32_t j = 0; j < i; j++)
            {
                s += decomp.l[i][j] * y[j];
            }
            y[i] = matrix->C[i] - s;
        }

        matrix->X = (double*)malloc(n * sizeof(double));
        // Решение Ux = y
        for (int32_t i = n - 1, k = 0; i >= 0; --i, k++)
        {
            double s = 0.0;
#pragma omp parallel for reduction(+:s) schedule(static)
            for (int32_t j = n - 1; j > i; --j)
            {
                s += decomp.u[i][j] * matrix->X[n - j - 1];
            }
            matrix->X[k] = (y[i] - s) / decomp.u[i][i];
        }

        reverse_array(matrix->X, n);
        free(y);
    }

   /* void solve_lu(DecomposeMatrix decomp, struct Matrix* matrix)
    {
        int32_t n = matrix->n;
        double* y = (double*)malloc(n * sizeof(double));

        for (int32_t i = 0; i < n; i++)
        {
            double s = 0.0;
            #pragma omp parallel for reduction(+:s)
            for (int32_t j = 0; j < i; j++)
            {
                s += decomp.l[i][j] * y[j];
            }
            y[i] = matrix->C[i] - s;
        }

        matrix->X = (double*)malloc(n * sizeof(double));
        for (int32_t i = n - 1, k = 0; i >= 0; --i, k++)
        {
            double s = 0.0;
            #pragma omp parallel for reduction(+:s)
            for (int32_t j = n - 1; j > i; --j)
            {
                s += decomp.u[i][j] * matrix->X[n - j - 1];
            }
            matrix->X[k] = (y[i] - s) / decomp.u[i][i];
        }

        reverse_array(matrix->X, n);
        free(y);
    }*/
}
