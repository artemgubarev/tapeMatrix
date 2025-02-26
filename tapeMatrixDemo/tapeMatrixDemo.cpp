
#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <mpi.h>
#include <omp.h>
#include "comparator.h"
#include "../tapeMatrix/matrix.h"
#include "../tapeMatrix/printer.h"
#include "../tapeMatrix/writer.h"
#include "../tapeMatrix/solver_serial.h"
#include "../tapeMatrix/solver_omp.h"

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

void reverse_array(double* array, size_t n)
{
	double temp;
	for (size_t i = 0; i < n / 2; ++i)
	{
		temp = array[i];
		array[i] = array[n - 1 - i];
		array[n - 1 - i] = temp;
	}
}


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

void random_fill(double** matrix, int32_t size)
{
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			matrix[i][j] = ((rand() % 10) + 1);
		}
	}

	//Ensure the matrix is diagonal dominant to guarantee invertible-ness
	//diagCount well help keep track of which column the diagonal is in
	int diagCount = 0;
	double sum = 0;
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			//Sum all column vaalues
			sum += abs(matrix[i][j]);
		}
		//Remove the diagonal  value from the sum
		sum -= abs(matrix[i][diagCount]);
		//Add a random value to the sum and place in diagonal position
		matrix[i][diagCount] = sum + ((rand() % 5) + 1);
		++diagCount;
		sum = 0;
	}
}

int main(int argc, char* argv[])
{
	const char* filename = getenv("INPUT_MATRIX_FILE");
	if (!filename)
	{
		fprintf(stderr, "Error: environment variable INPUT_MATRIX_FILE not set.\n");
		return 1;
	}

	Matrix matrix;
	matrix = read_matrix_mpi(filename);
	//matrix = read_matrix_mpi("testData/matrix2000.txt");

	clock_t start, end;
	double cpu_time_used;
	start = clock();
	DecomposeMatrix decomp = band_matrix_omp::lu_decomposition(matrix);
	band_matrix_omp::solve_lu(decomp, &matrix);
	end = clock();
	cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
	printf("time: %f sec\n", cpu_time_used);

	//print_1d(matrix.X, matrix.n);

		//MPI_Init(&argc, &argv);
	
		//int rank, size;
		//MPI_Comm comm = MPI_COMM_WORLD;
		//MPI_Comm_rank(comm, &rank);
		//MPI_Comm_size(comm, &size);
	
		//Matrix matrix;
		//matrix.X = NULL;
		//uint32_t n = 0, b = 0;
	
		//if (rank == 0)
		//{
		//	matrix = read_matrix_mpi("testData/matrix5.txt");
		//	//matrix = read_matrix_mpi(filename);
		//	n = matrix.n;
		//	b = matrix.b;
		//}
	
		////print_1d(matrix.C, matrix.n);
	
		//int32_t n_int = static_cast<int32_t>(n);
		//int32_t b_int = static_cast<int32_t>(b);
	
		//MPI_Bcast(&n_int, 1, MPI_INT, 0, comm);
		//MPI_Bcast(&b_int, 1, MPI_INT, 0, comm);
	
		//n = static_cast<uint32_t>(n_int);
		//b = static_cast<uint32_t>(b_int);
	
		//if (b == 0 || b > n) 
		//{
		//	fprintf(stderr, "Error: Invalid block size b = %u (n = %u)\n", b, n);
		//	MPI_Abort(comm, EXIT_FAILURE);
		//}
	
		//if (rank != 0) {
		//	matrix.n = n;
		//	matrix.b = b;
		//	matrix.A = allocate_matrix_mpi(n);
		//	matrix.C = (double*)malloc(n * sizeof(double));
		//	if (!matrix.C)
		//	{
		//		perror("Ошибка выделения вектора C на не-root процессе");
		//		MPI_Abort(comm, EXIT_FAILURE);
		//	}
		//	matrix.X = NULL;
		//}
	
		//MPI_Bcast(matrix.A[0], n * n, MPI_DOUBLE, 0, comm);
		//MPI_Bcast(matrix.C, n, MPI_DOUBLE, 0, comm);
	
		//double start_time_mpi = MPI_Wtime();
	
		///*DecomposeMatrix decomp = band_matrix_serial::lu_decomposition(matrix);
		//band_matrix_serial::solve_lu(decomp, &matrix);*/
	
		//double end_time_mpi = MPI_Wtime();
		//if (rank == 0)
		//{
		//	printf("Time: %f sec.\n", end_time_mpi - start_time_mpi);
		//	//write_1d("solution.txt", matrix.X, matrix.n);
		//	//print_1d(matrix.X, matrix.n);
		//}
	
	#pragma region compare
		/*double epsilon = 0.00001;
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
		delete[] numbers2;*/
	#pragma endregion
	
		/*free(matrix.A[0]);
		free(matrix.A);
		free(matrix.C);
		free(matrix.X);
	
		MPI_Finalize();*/
		return 0;

}