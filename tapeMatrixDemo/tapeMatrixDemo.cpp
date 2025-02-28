#define HAVE_STRUCT_TIMESPEC
#define _NO_DEBUG_HEAP 1

#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <mpi.h>
#include <omp.h>
#include <time.h>
#include <pthread.h>
#include "comparator.h"
#include "../tapeMatrix/matrix.h"
#include "../tapeMatrix/printer.h"
#include "../tapeMatrix/writer.h"
#include "../tapeMatrix/solver_serial.h"
#include "../tapeMatrix/solver_omp.h"
#include "../tapeMatrix/solver_mpi.h"
#include "../tapeMatrix/solver_pthreads.h"

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

void get_output_filename(const char* input_file, char* output_filename, long size)
{
	const char* filename = strrchr(input_file, '/');
	filename = (filename) ? filename + 1 : input_file;

	const char* name_part = strstr(filename, "matrix");
	if (name_part)
	{
		name_part += 6;
	}
	else
	{
		name_part = filename;
	}

	char name_only[256];
	strncpy(name_only, name_part, sizeof(name_only) - 1);
	name_only[sizeof(name_only) - 1] = '\0';

	char* dot = strrchr(name_only, '.');
	if (dot)
	{
		*dot = '\0';
	}

	snprintf(output_filename, size, "mSolutions/msolution%s.txt", name_only);
}

double get_time()
{
#ifdef _WIN32
	LARGE_INTEGER frequency, start;
	QueryPerformanceFrequency(&frequency);
	QueryPerformanceCounter(&start);
	return (double)start.QuadPart / frequency.QuadPart;
#else
	struct timespec t;
	clock_gettime(CLOCK_MONOTONIC, &t);
	return t.tv_sec + t.tv_nsec * 1e-9;
#endif
}

int main(int argc, char* argv[])
{
#pragma region get slurm variables

	const char* filename = getenv("INPUT_MATRIX_FILE");
	const char* num_threads_str = getenv("SLURM_CPUS_PER_TASK");
	const char* mode_str = getenv("MODE");

	if (!filename || !num_threads_str || !mode_str) {
		std::cerr << "Error: Required environment variables are not set." << std::endl;
		return 1;
	}

	int32_t num_threads = atoi(num_threads_str);
	int32_t mode = atoi(mode_str);

	if (num_threads <= 0) {
		std::cerr << "Error: Invalid number of threads." << std::endl;
		return 1;
	}

#pragma endregion

#pragma region run solver

    Matrix matrix;
    double start_time, end_time;

    if (mode == 3) 
    {
        MPI_Init(&argc, &argv);
    }

    if (mode != 3) 
    {
        matrix = read_matrix(filename);
    }

    switch (mode) 
    {
    case 0:  // serial
        start_time = (double)clock() / CLOCKS_PER_SEC;
        DecomposeMatrix decomp_serial = tape_matrix_serial::lu_decomposition(matrix);
        tape_matrix_serial::solve_lu(decomp_serial, &matrix);
        end_time = (double)clock() / CLOCKS_PER_SEC;
        break;

    case 1:  // pthreads
        start_time = get_time();
        DecomposeMatrix decomp_pthreads = tape_matrix_pthreads::lu_decomposition(matrix, num_threads);
        tape_matrix_pthreads::solve_lu(decomp_pthreads, &matrix, num_threads);
        end_time = get_time();
        break;

    case 2:  // openmp
        start_time = omp_get_wtime();
        DecomposeMatrix decomp_omp = tape_matrix_omp::lu_decomposition(matrix);
        tape_matrix_omp::solve_lu(decomp_omp, &matrix);
        end_time = omp_get_wtime();
        break;

    case 3:  // MPI
        int32_t rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        if (rank == 0) {
            matrix = tape_matrix_mpi::read_matrix(filename, rank);
        }

        MPI_Bcast(&matrix.n, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&matrix.b, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (rank != 0) 
        {
            matrix.A = (double**)malloc(matrix.n * sizeof(double*));
            matrix.A[0] = (double*)malloc(matrix.n * matrix.n * sizeof(double));
            for (int i = 1; i < matrix.n; i++) 
            {
                matrix.A[i] = matrix.A[0] + i * matrix.n;
            }
            matrix.C = (double*)malloc(matrix.n * sizeof(double));
        }

        MPI_Bcast(matrix.A[0], matrix.n * matrix.n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(matrix.C, matrix.n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        start_time = MPI_Wtime();
        DecomposeMatrix decomp_mpi = tape_matrix_mpi::lu_decomposition(matrix, rank, size);
        tape_matrix_mpi::solve_lu(decomp_mpi, &matrix, rank, size);
        end_time = MPI_Wtime();

        if (rank == 0) 
        {
            printf("MPI time: %.6f sec\n", end_time - start_time);
            write_1d("solution.txt", matrix.X, matrix.n);
        }
        break;

    default:
        printf("Incorrect mode value\n");
        if (mode == 3) MPI_Finalize();
        return 1;
    }

    if (mode != 3)
    {
        write_1d("solution.txt", matrix.X, matrix.n);
    }

    if (mode == 3) 
    {
        MPI_Finalize();
    }

    if (mode != 3) 
    {
        printf("Time spent: %.6f seconds\n", end_time - start_time);
    }


#pragma endregion

#pragma region compare
	double epsilon = 0.00001;
	double* numbers1 = new double[MAX_NUMBERS];
	double* numbers2 = new double[MAX_NUMBERS];
	size_t count1, count2;
	char output_filename[512];
	get_output_filename(filename, output_filename, sizeof(output_filename));
	if (!load_numbers("solution.txt", numbers1, &count1) ||
		!load_numbers(output_filename, numbers2, &count2)) {
		return 1;
	}
	if (compare_numbers(numbers1, numbers2, count1, count2, epsilon)) {
		printf("\033[32mTest Correct\033[0m\n");
	}
	else {
		printf("\033[31mTest Failed\033[0m\n");
	}
	delete[] numbers1;
	delete[] numbers2;
#pragma endregion

	return 0;
}