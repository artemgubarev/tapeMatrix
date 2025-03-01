from mpi4py import MPI
import time
import sys
import os
sys.path.append("../tapeMatrix")
import matrix as mtrx
import solver_mpi as solver
import comparator as comp
import numpy as np

def get_output_filename(input_filename, suffix="_output"):
    base, ext = os.path.splitext(input_filename)
    return f"{base}{suffix}{ext}"

def test_compare_files(filename, solution_filename="solution.txt", epsilon=1e-5):
    numbers1, count1 = comp.load_numbers(solution_filename)
    if numbers1 is None:
        return 1
    
    output_filename = get_output_filename(filename)
    numbers2, count2 = comp.load_numbers(output_filename)
    if numbers2 is None:
        return 1

    if comp.compare_numbers(numbers1, numbers2, count1, count2, epsilon):
        print("\033[32mTest Correct\033[0m")
        return 0
    else:
        print("\033[31mTest Failed\033[0m")
        return 1

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

filename = os.environ['INPUT_MATRIX_FILE']
#filename = 'testData/matrix2000.txt'
matrix = mtrx.read_matrix(filename, rank)

if rank == 0:
    start_time = time.time()

decompose_matrix = solver.lu_decomposition(matrix, rank, size)
solver.solve_lu(decompose_matrix, matrix, rank, size)

if rank == 0:
    end_time = time.time()
    elapsed_time = end_time - start_time
    np.savetxt("output.txt", matrix.X, fmt='%.6f', delimiter=' ')
    #print(f"Solution X: {matrix.X}")
    print(f"Elapsed time: {elapsed_time:.6f} seconds")
    test_compare_files(filename)