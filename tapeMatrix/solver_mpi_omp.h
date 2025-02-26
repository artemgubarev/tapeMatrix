#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "matrix.h"

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif

namespace band_matrix_mpi_omp
{
    DecomposeMatrix lu_decomposition(const Matrix matrix)
    {
        int n = matrix.n;
        int b = matrix.b;

        // ��������� ������ ��� L � U.
        DecomposeMatrix result;
        result.l = (double**)malloc(n * sizeof(double*));
        result.u = (double**)malloc(n * sizeof(double*));
        for (int i = 0; i < n; i++) {
            result.l[i] = (double*)calloc(n, sizeof(double));
            result.u[i] = (double*)malloc(n * sizeof(double));
        }

        // ������������ ������������� ������ L � U.
#pragma omp parallel for schedule(static)
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                result.u[i][j] = matrix.A[i][j];
            }
            result.l[i][i] = 1.0;
        }

        // ������� ���� LU-����������.
        for (int k = 0; k < n - 1; k++) {
            double pivotVal = result.u[k][k];

            // ���������� ����� �� k+1 �� k+b (�������� ������� n).
#pragma omp parallel for schedule(static)
            for (int i = k + 1; i < MIN(k + b + 1, n); i++)
            {
                result.l[i][k] = result.u[i][k] / pivotVal;
                for (int j = k; j < MIN(k + b + 1, n); j++)
                {
                    result.u[i][j] -= result.l[i][k] * result.u[k][j];
                }
            }
        }

        return result;
    }

    // ������� ������� ������� ������� LU-����������.
    void solve_lu(const DecomposeMatrix decomp, Matrix* matrix)
    {
        int n = matrix->n;
        double* y = (double*)malloc(n * sizeof(double));

        // ������ ���: ������� ������� L*y = C.
        // ������ ���� �������� ���� ����������������, �.�. y[i] ������� �� y[0..i-1].
        for (int i = 0; i < n; i++)
        {
            double s = 0.0;
            // ���������� ���� ������������ ����� ��������������, ���� ����� ��������� ���������� ������.
#pragma omp parallel for reduction(+:s) schedule(static)
            for (int j = 0; j < i; j++)
            {
                s += decomp.l[i][j] * y[j];
            }
            y[i] = matrix->C[i] - s;
        }

        // �������� ���: ������� ������� U*x = y.
        matrix->X = (double*)malloc(n * sizeof(double));
        for (int i = (int)n - 1; i >= 0; i--)
        {
            double s = 0.0;
#pragma omp parallel for reduction(+:s) schedule(static)
            for (int j = i + 1; j < (int)n; j++)
            {
                s += decomp.u[i][j] * matrix->X[j];
            }
            matrix->X[i] = (y[i] - s) / decomp.u[i][i];
        }

        free(y);
    }
}
