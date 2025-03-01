from mpi4py import MPI
import time
import sys
import os
sys.path.append("../tapeMatrix")
import matrix as mtrx
import solver_mpi as solver
import comparator as comp

def get_output_filename(input_filename):
    base, ext = os.path.splitext(input_filename)
    return f"{base}{ext}"

def write_sol(filename, numbers, precision=5):
    with open(filename, 'w') as file:
        formatted_numbers = [f"{num:.{precision}f}" for num in numbers]
        file.write(" ".join(formatted_numbers))

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
    print(elapsed_time)
    write_sol("solution.txt", matrix.X)
    comp.compare_numbers("solution.txt", get_output_filename(filename))