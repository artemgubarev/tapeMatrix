#pragma once
#include <cstdint>

struct DecomposeMatrix
{
    double** l;
    double** u;
};

struct Matrix
{
    int32_t n;
    int32_t b;
    double* C;
    double** A;
    double* X;
};

double** allocate_matrix(uint32_t n)
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

Matrix read_matrix(const char* filename)
{
	Matrix mat;
	FILE* file = fopen(filename, "r");
	if (!file)
	{
		perror("Error opening file");
		exit(EXIT_FAILURE);
	}

	if (fscanf(file, "%u %u", &mat.n, &mat.b) != 2)
	{
		fprintf(stderr, "Error reading matrix size and band width\n");
		exit(EXIT_FAILURE);
	}

	mat.A = allocate_matrix(mat.n);
	mat.C = (double*)malloc(mat.n * sizeof(double));
	if (!mat.C)
	{
		perror("Error allocating vector C");
		exit(EXIT_FAILURE);
	}
	mat.X = NULL;

	uint64_t total_matrix_elements = mat.n * mat.n;
	for (uint64_t idx = 0; idx < total_matrix_elements; idx++)
	{
		if (fscanf(file, "%lf", &mat.A[0][idx]) != 1)
		{
			fprintf(stderr, "Error reading matrix data at element %llu\n", idx);
			exit(EXIT_FAILURE);
		}
	}

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