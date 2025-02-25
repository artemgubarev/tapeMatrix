//#include <iostream>
//#include <cstdlib>
//#include <stdio.h>
//#include <stdlib.h>
//#include <cstring>
//#include <time.h>
//#include <mpi.h>
//#include <omp.h>
//#include "../tapeMatrix/matrix.h"
//#include "../tapeMatrix/printer.h"
//
//double** allocate_matrix_mpi(int n)
//{
//	double** matrix = (double**)malloc(n * sizeof(double*));
//	if (!matrix)
//	{
//		perror("Error allocating matrix pointers");
//		exit(EXIT_FAILURE);
//	}
//	matrix[0] = (double*)malloc(n * n * sizeof(double));
//	if (!matrix[0])
//	{
//		perror("Error allocating matrix memory");
//		exit(EXIT_FAILURE);
//	}
//	for (int i = 1; i < n; i++)
//	{
//		matrix[i] = matrix[0] + i * n;
//	}
//	return matrix;
//}
//
//Matrix read_matrix_mpi(const char* filename) {
//	Matrix mat;
//	FILE* file = fopen(filename, "r");
//	if (!file) {
//		perror("Error opening file");
//		exit(EXIT_FAILURE);
//	}
//
//	// Read matrix size first!
//	if (fscanf(file, "%ld", &mat.n) != 1) {
//		fprintf(stderr, "Error reading matrix size\n");
//		exit(EXIT_FAILURE);
//	}
//
//	mat.A = allocate_matrix_mpi(mat.n);
//	mat.C = (double*)malloc(mat.n * sizeof(double));
//	if (!mat.C) {
//		perror("Error allocating vector C");
//		exit(EXIT_FAILURE);
//	}
//	mat.X = NULL;
//
//	long total_matrix_elements = mat.n * mat.n;
//	for (long idx = 0; idx < total_matrix_elements; idx++)
//	{
//		if (fscanf(file, "%lf", &mat.A[0][idx]) != 1)
//		{
//			fprintf(stderr,
//				"Error reading matrix data at element %ld (expected %ld total)\n",
//				idx, total_matrix_elements);
//			exit(EXIT_FAILURE);
//		}
//	}
//
//	for (long i = 0; i < mat.n; i++)
//	{
//		if (fscanf(file, "%lf", &mat.C[i]) != 1)
//		{
//			fprintf(stderr,
//				"Error reading vector C at index %ld (expected %ld elements)\n",
//				i, mat.n);
//			exit(EXIT_FAILURE);
//		}
//	}
//
//	fclose(file);
//	return mat;
//}
//
//DecomposeMatrix lu_decomposition(Matrix matrix) 
//{
//	int rank, size;
//	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//	MPI_Comm_size(MPI_COMM_WORLD, &size);
//
//	int n = matrix.n;
//	int b = matrix.b;
//
//	DecomposeMatrix result;
//	result.l = (double**)malloc(n * sizeof(double*));
//	result.u = (double**)malloc(n * sizeof(double*));
//	for (int i = 0; i < n; i++) 
//	{
//		result.l[i] = (double*)calloc(n, sizeof(double));
//		result.u[i] = (double*)malloc(n * sizeof(double));
//	}
//
//	for (int i = 0; i < n; i++) 
//	{
//		for (int j = 0; j < n; j++) 
//		{
//			result.u[i][j] = matrix.A[i][j];
//		}
//	}
//
//	for (int i = 0; i < n; i++) 
//	{
//		result.l[i][i] = 1.0;
//	}
//
//	int start_row = (n * rank) / size;
//	int end_row = (n * (rank + 1)) / size;
//
//	for (int k = 0; k < n - 1; k++) 
//	{
//		int upper_bound = (k + b + 1 < n) ? (k + b + 1) : n;
//		int owner = (k * size) / n;
//		int segment_length = upper_bound - k;
//		double* u_row_segment = (double*)malloc(segment_length * sizeof(double));
//		if (rank == owner) {
//			for (int j = k; j < upper_bound; j++) {
//				u_row_segment[j - k] = result.u[k][j];
//			}
//		}
//		MPI_Request req;
//		MPI_Ibcast(u_row_segment, segment_length, MPI_DOUBLE, owner, MPI_COMM_WORLD, &req);
//		MPI_Wait(&req, MPI_STATUS_IGNORE);
//
//		int local_start = (k + 1 > start_row) ? k + 1 : start_row;
//		for (int i = local_start; i < end_row && i < upper_bound; i++) {
//			result.l[i][k] = result.u[i][k] / u_row_segment[0];
//			for (int j = k; j < upper_bound; j++) {
//				result.u[i][j] -= result.l[i][k] * u_row_segment[j - k];
//			}
//		}
//		free(u_row_segment);
//	}
//	return result;
//}
//
//void solve_lu(DecomposeMatrix decompose_matrix, Matrix* matrix) 
//{
//	int rank, size;
//	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//	MPI_Comm_size(MPI_COMM_WORLD, &size);
//	int n = matrix->n;
//
//	int start_row = (n * rank) / size;
//	int end_row = (n * (rank + 1)) / size;
//
//	double* y = (double*)malloc(n * sizeof(double));
//
//	for (int i = 0; i < n; i++) 
//	{
//		if (i >= start_row && i < end_row) 
//		{
//			double s = 0.0;
//			for (int j = 0; j < i; j++) 
//			{
//				s += decompose_matrix.l[i][j] * y[j];
//			}
//			y[i] = matrix->C[i] - s;
//		}
//		int owner = (i * size) / n;
//		MPI_Request req;
//		MPI_Ibcast(&y[i], 1, MPI_DOUBLE, owner, MPI_COMM_WORLD, &req);
//		MPI_Wait(&req, MPI_STATUS_IGNORE);
//	}
//
//	if (matrix->X == NULL) 
//	{
//		matrix->X = (double*)malloc(n * sizeof(double));
//	}
//	for (int i = n - 1; i >= 0; i--) 
//	{
//		if (i >= start_row && i < end_row) 
//		{
//			double s = 0.0;
//			for (int j = i + 1; j < n; j++) 
//			{
//				s += decompose_matrix.u[i][j] * matrix->X[j];
//			}
//			matrix->X[i] = (y[i] - s) / decompose_matrix.u[i][i];
//		}
//		int owner = (i * size) / n;
//		MPI_Request req;
//		MPI_Ibcast(&matrix->X[i], 1, MPI_DOUBLE, owner, MPI_COMM_WORLD, &req);
//		MPI_Wait(&req, MPI_STATUS_IGNORE);
//	}
//
//	free(y);
//}
//
//int main(int argc, char* argv[])
//{
//	/*const char* filename = getenv("INPUT_MATRIX_FILE");
//	if (!filename)
//	{
//		fprintf(stderr, "Error: environment variable INPUT_MATRIX_FILE not set.\n");
//		return 1;
//	}*/
//
//	MPI_Init(&argc, &argv);
//
//	int rank, size;
//	MPI_Comm comm = MPI_COMM_WORLD;
//	MPI_Comm_rank(comm, &rank);
//	MPI_Comm_size(comm, &size);
//
//	Matrix matrix;
//	matrix.X = NULL;
//	long n = 0, b = 0;
//
//	if (rank == 0)
//	{
//		matrix = read_matrix_mpi("testData\\matrix100.txt");
//		//matrix = read_matrix_mpi(filename);
//		n = matrix.n;
//		b = matrix.b;
//	}
//
//	MPI_Bcast(&n, 1, MPI_LONG, 0, comm);
//	MPI_Bcast(&b, 1, MPI_LONG, 0, comm);
//
//	if (rank != 0) {
//		matrix.n = n;
//		matrix.b = b;
//		matrix.A = allocate_matrix_mpi(n);
//		matrix.C = (double*)malloc(n * sizeof(double));
//		if (!matrix.C)
//		{
//			perror("Ошибка выделения вектора C на не-root процессе");
//			MPI_Abort(comm, EXIT_FAILURE);
//		}
//		matrix.X = NULL;
//	}
//
//	MPI_Bcast(matrix.A[0], n * n, MPI_DOUBLE, 0, comm);
//	MPI_Bcast(matrix.C, n, MPI_DOUBLE, 0, comm);
//
//	double start_time_mpi = MPI_Wtime();
//
//	DecomposeMatrix decomp = lu_decomposition(matrix);
//	solve_lu(decomp, &matrix);
//
//	double end_time_mpi = MPI_Wtime();
//	if (rank == 0)
//	{
//		printf("Time: %f sec.\n", end_time_mpi - start_time_mpi);
//	}
//
//	free(matrix.A[0]);
//	free(matrix.A);
//	free(matrix.C);
//	free(matrix.X);
//
//	for (int i = 0; i < n; i++)
//	{
//		free(decomp.l[i]);
//		free(decomp.u[i]);
//	}
//	free(decomp.l);
//	free(decomp.u);
//
//	MPI_Finalize();
//	return 0;
//}


#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <mpi.h>
#include <omp.h>
#include "../tapeMatrix/matrix.h"
#include "../tapeMatrix/printer.h"

double** allocate_matrix_mpi(uint32_t n)
{
    double** matrix = (double**)malloc(n * sizeof(double*));
    if (!matrix)
    {
        perror("Error allocating matrix pointers");
        exit(EXIT_FAILURE);
    }
    matrix[0] = (double*)malloc(n * n * sizeof(double));
    if (!matrix[0])
    {
        perror("Error allocating matrix memory");
        exit(EXIT_FAILURE);
    }
    for (size_t i = 1; i < n; i++)
    {
        matrix[i] = matrix[0] + i * n;
    }
    return matrix;
}

Matrix read_matrix_mpi(const char* filename)
{
    Matrix mat;
    FILE* file = fopen(filename, "r");
    if (!file) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    // Читаем размер матрицы и ширину ленты
    if (fscanf(file, "%u %u", &mat.n, &mat.b) != 2) {
        fprintf(stderr, "Error reading matrix size and band width\n");
        exit(EXIT_FAILURE);
    }

    mat.A = allocate_matrix_mpi(mat.n);
    mat.C = (double*)malloc(mat.n * sizeof(double));
    if (!mat.C) {
        perror("Error allocating vector C");
        exit(EXIT_FAILURE);
    }
    mat.X = NULL;

    // Читаем матрицу
    uint64_t total_matrix_elements = mat.n * mat.n;
    for (uint64_t idx = 0; idx < total_matrix_elements; idx++)
    {
        if (fscanf(file, "%lf", &mat.A[0][idx]) != 1)
        {
            fprintf(stderr, "Error reading matrix data at element %llu\n", idx);
            exit(EXIT_FAILURE);
        }
    }

    // Читаем вектор C
    for (size_t i = 0; i < mat.n; i++)
    {
        if (fscanf(file, "%lf", &mat.C[i]) != 1)
        {
            fprintf(stderr, "Error reading vector C at index %zu\n", i);
            exit(EXIT_FAILURE);
        }
    }

    fclose(file);
    return mat;
}

DecomposeMatrix lu_decomposition(Matrix matrix)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    uint32_t n = matrix.n;
    uint32_t b = matrix.b;

    DecomposeMatrix result;
    result.l = (double**)malloc(n * sizeof(double*));
    result.u = (double**)malloc(n * sizeof(double*));
    for (size_t i = 0; i < n; i++)
    {
        result.l[i] = (double*)calloc(n, sizeof(double));
        result.u[i] = (double*)malloc(n * sizeof(double));
    }

    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < n; j++)
        {
            result.u[i][j] = matrix.A[i][j];
        }
    }

    for (size_t i = 0; i < n; i++)
    {
        result.l[i][i] = 1.0;
    }

    uint32_t start_row = (n * rank) / size;
    uint32_t end_row = (n * (rank + 1)) / size;

    for (size_t k = 0; k < n - 1; k++)
    {
        uint32_t upper_bound = (k + b + 1 < n) ? (k + b + 1) : n;
        uint32_t owner = (k * size) / n;
        uint32_t segment_length = upper_bound - k;

        if (segment_length > 0) {
            double* u_row_segment = (double*)malloc(segment_length * sizeof(double));
            if (rank == owner) {
                for (size_t j = k; j < upper_bound; j++) {
                    u_row_segment[j - k] = result.u[k][j];
                }
            }
            MPI_Request req;
            MPI_Ibcast(u_row_segment, segment_length, MPI_DOUBLE, owner, MPI_COMM_WORLD, &req);
            MPI_Wait(&req, MPI_STATUS_IGNORE);
            free(u_row_segment);
        }
    }
    return result;
}

void solve_lu(DecomposeMatrix decompose_matrix, Matrix* matrix)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    uint32_t n = matrix->n;

    uint32_t start_row = (n * rank) / size;
    uint32_t end_row = (n * (rank + 1)) / size;

    double* y = (double*)malloc(n * sizeof(double));
    if (!y) {
        perror("Error allocating vector y");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    for (uint32_t i = 0; i < n; i++)
    {
        if (i >= start_row && i < end_row)
        {
            double s = 0.0;
            for (uint32_t j = 0; j < i; j++)
            {
                s += decompose_matrix.l[i][j] * y[j];
            }
            y[i] = matrix->C[i] - s;
        }

        uint32_t owner = (i * size) / n;
        MPI_Request req;
        MPI_Ibcast(&y[i], 1, MPI_DOUBLE, owner, MPI_COMM_WORLD, &req);
        MPI_Wait(&req, MPI_STATUS_IGNORE);
    }

    if (matrix->X == NULL)
    {
        matrix->X = (double*)malloc(n * sizeof(double));
        if (!matrix->X) {
            perror("Error allocating vector X");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
    }

    for (int32_t i = n - 1; i >= 0; i--)
    {
        if (i >= static_cast<int32_t>(start_row) && i < static_cast<int32_t>(end_row))
        {
            double s = 0.0;
            for (uint32_t j = i + 1; j < n; j++)
            {
                s += decompose_matrix.u[i][j] * matrix->X[j];
            }
            matrix->X[i] = (y[i] - s) / decompose_matrix.u[i][i];
        }

        uint32_t owner = (i * size) / n;
        MPI_Request req;
        MPI_Ibcast(&matrix->X[i], 1, MPI_DOUBLE, owner, MPI_COMM_WORLD, &req);
        MPI_Wait(&req, MPI_STATUS_IGNORE);
    }

    free(y);
}

int main(int argc, char* argv[])
{
    const char* filename = getenv("INPUT_MATRIX_FILE");
    if (!filename)
    {
        fprintf(stderr, "Error: environment variable INPUT_MATRIX_FILE not set.\n");
        return 1;
    }

    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    Matrix matrix;
    matrix.X = NULL;
    uint32_t n = 0, b = 0;

    if (rank == 0)
    {
        //matrix = read_matrix_mpi("testData/matrix100.txt");
        matrix = read_matrix_mpi(filename);
        n = matrix.n;
        b = matrix.b;
    }

    int32_t n_int = static_cast<int32_t>(n);
    int32_t b_int = static_cast<int32_t>(b);

    MPI_Bcast(&n_int, 1, MPI_INT, 0, comm);
    MPI_Bcast(&b_int, 1, MPI_INT, 0, comm);

    n = static_cast<uint32_t>(n_int);
    b = static_cast<uint32_t>(b_int);

    if (b == 0 || b > n) {
        fprintf(stderr, "Error: Invalid block size b = %u (n = %u)\n", b, n);
        MPI_Abort(comm, EXIT_FAILURE);
    }

    if (rank != 0) {
        matrix.n = n;
        matrix.b = b;
        matrix.A = allocate_matrix_mpi(n);
        matrix.C = (double*)malloc(n * sizeof(double));
        if (!matrix.C)
        {
            perror("Ошибка выделения вектора C на не-root процессе");
            MPI_Abort(comm, EXIT_FAILURE);
        }
        matrix.X = NULL;
    }

    MPI_Bcast(matrix.A[0], n * n, MPI_DOUBLE, 0, comm);
    MPI_Bcast(matrix.C, n, MPI_DOUBLE, 0, comm);

    double start_time_mpi = MPI_Wtime();

    DecomposeMatrix decomp = lu_decomposition(matrix); 
    solve_lu(decomp, &matrix);

    double end_time_mpi = MPI_Wtime();
    if (rank == 0)
    {
        printf("Time: %f sec.\n", end_time_mpi - start_time_mpi);
    }

    free(matrix.A[0]);
    free(matrix.A);
    free(matrix.C);
    free(matrix.X);

    MPI_Finalize();
    return 0;
}