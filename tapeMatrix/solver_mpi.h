#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <utility>
#include <mpi.h>
#include "matrix.h"

namespace band_matrix_mpi
{
    void reverse_array(double* array, int32_t n)
    {
        double temp;
        for (int32_t i = 0; i < n / 2; ++i)
        {
            temp = array[i];
            array[i] = array[n - 1 - i];
            array[n - 1 - i] = temp;
        }
    }

    struct DecomposeMatrix lu_decomposition(struct Matrix matrix, int32_t rank, int32_t size)
    {
        struct DecomposeMatrix result;
        int32_t n = matrix.n;
        int32_t b = matrix.b;

        // Выделение памяти
        result.l = (double**)malloc(n * sizeof(double*));
        result.u = (double**)malloc(n * sizeof(double*));
        for (int32_t i = 0; i < n; i++) {
            result.l[i] = (double*)calloc(n, sizeof(double));
            result.u[i] = (double*)malloc(n * sizeof(double));
            result.l[i][i] = 1.0;  // Диагональ L = 1
            for (int32_t j = 0; j < n; j++) {
                result.u[i][j] = matrix.A[i][j];
            }
        }

        // Определение блока работы для каждого процесса
        int32_t rows_per_process = (n + size - 1) / size;  // Округление вверх
        int32_t start_row = rank * rows_per_process;
        int32_t end_row = (rank + 1) * rows_per_process;
        if (end_row > n) end_row = n;

        // Буферы для коммуникации
        double* pivot_row = (double*)malloc(n * sizeof(double));
        double* l_column = (double*)malloc(n * sizeof(double));

        // Инициализация буфера l_column нулями
        for (int32_t i = 0; i < n; i++) {
            l_column[i] = 0.0;
        }

        for (int32_t k = 0; k < n - 1; k++) {
            int32_t pivot_owner = k / rows_per_process;
            int32_t upper_bound = (k + b + 1 < n) ? (k + b + 1) : n;

            // Процесс-владелец строки k рассылает её всем
            if (rank == pivot_owner) {
                for (int32_t j = k; j < n; j++) {
                    pivot_row[j] = result.u[k][j];
                }
            }
            MPI_Bcast(&pivot_row[k], n - k, MPI_DOUBLE, pivot_owner, MPI_COMM_WORLD);

            // Каждый процесс обрабатывает свои строки
            for (int32_t i = start_row; i < end_row; i++) {
                if (i > k && i < upper_bound) {
                    result.l[i][k] = result.u[i][k] / pivot_row[k];
                    for (int32_t j = k; j < upper_bound; j++) {
                        result.u[i][j] -= result.l[i][k] * pivot_row[j];
                    }
                }
            }

            // Сбрасываем буфер перед использованием
            for (int32_t i = 0; i < n; i++) {
                l_column[i] = 0.0;
            }

            // Заполняем буфер локальными значениями L
            for (int32_t i = start_row; i < end_row; i++) {
                if (i > k && i < upper_bound) {
                    l_column[i] = result.l[i][k];
                }
            }

            // Объединяем значения L со всех процессов
            double* global_l = (double*)malloc(n * sizeof(double));
            MPI_Allreduce(l_column, global_l, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

            // Обновляем локальные копии L
            for (int32_t i = k + 1; i < upper_bound; i++) {
                result.l[i][k] = global_l[i];
            }

            free(global_l);
        }

        free(pivot_row);
        free(l_column);
        return result;
    }

    void solve_lu(struct DecomposeMatrix decompose_matrix, struct Matrix* matrix, int32_t rank, int32_t size)
    {
        int32_t n = matrix->n;
        double* y = (double*)calloc(n, sizeof(double));
        double* global_y = (double*)calloc(n, sizeof(double));

        int32_t rows_per_process = (n + size - 1) / size;
        int32_t start_row = rank * rows_per_process;
        int32_t end_row = (rank + 1) * rows_per_process;
        if (end_row > n) end_row = n;

        // Forward substitution (Ly = b)
        for (int32_t i = 0; i < n; i++) {
            // Синхронизация после каждой итерации
            MPI_Bcast(global_y, i, MPI_DOUBLE, 0, MPI_COMM_WORLD);

            if (i >= start_row && i < end_row) {
                double s = 0.0;
                for (int32_t j = 0; j < i; j++) {
                    s += decompose_matrix.l[i][j] * global_y[j];
                }
                y[i] = matrix->C[i] - s;
            }

            // Собираем текущий элемент y
            MPI_Allreduce(&y[i], &global_y[i], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        }

        // Backward substitution (Ux = y)
        if (rank == 0) {
            matrix->X = (double*)malloc(n * sizeof(double));
        }
        else {
            matrix->X = (double*)calloc(n, sizeof(double));
        }

        for (int32_t i = n - 1; i >= 0; i--) {
            double local_x = 0.0;
            if (i >= start_row && i < end_row) {
                double s = 0.0;
                for (int32_t j = i + 1; j < n; j++) {
                    s += decompose_matrix.u[i][j] * matrix->X[j];
                }
                local_x = (global_y[i] - s) / decompose_matrix.u[i][i];
            }

            double global_x = 0.0;
            MPI_Allreduce(&local_x, &global_x, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            matrix->X[i] = global_x;
        }

        free(y);
        free(global_y);
    }
}