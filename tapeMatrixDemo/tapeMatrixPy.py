from mpi4py import MPI
import time
import sys
import os
sys.path.append("../tapeMatrix")
import matrix as mtrx
import solver_mpi as solver

def write_sol(sol):
     with open('solution.txt', 'w') as f:
        formatted_values = [f"{x:.6f}" for x in sol]
        f.write(' '.join(formatted_values))

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
    write_sol(matrix.X)
    #print(f"Solution X: {matrix.X}")
    print(f"Elapsed time: {elapsed_time:.6f} seconds")