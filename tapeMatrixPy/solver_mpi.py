import numpy as np
from mpi4py import MPI
import matrix as mtrx

def reverse_array(array):
    array[:] = array[::-1]
    
def lu_decomposition(matrix, rank, size):
    comm = MPI.COMM_WORLD
    n = matrix.n
    b = matrix.b
    
    l = np.eye(n)
    u = np.copy(matrix.A)
    
    rows_per_process = (n + size - 1) // size
    start_row = rank * rows_per_process
    end_row = min((rank + 1) * rows_per_process, n)

    pivot_row = np.zeros(n)
    l_column = np.zeros(n)

    for k in range(n - 1):
        pivot_owner = k // rows_per_process
        upper_bound = min(k + b + 1, n)

        if rank == pivot_owner:
            pivot_row[k:] = u[k, k:]
        comm.Bcast(pivot_row[k:], root=pivot_owner)

        for i in range(start_row, end_row):
            if i > k and i < upper_bound:
                l[i, k] = u[i, k] / pivot_row[k]
                u[i, k:upper_bound] -= l[i, k] * pivot_row[k:upper_bound]

        l_column[:] = 0.0
        for i in range(start_row, end_row):
            if i > k and i < upper_bound:
                l_column[i] = l[i, k]

        global_l = np.zeros(n)
        comm.Allreduce(l_column, global_l, op=MPI.SUM)

        for i in range(k + 1, upper_bound):
            l[i, k] = global_l[i]
            
    return mtrx.DecomposeMatrix(l=l, u=u)


def solve_lu(decompose_matrix, matrix, rank, size):
    comm = MPI.COMM_WORLD
    n = matrix.n
    y = np.zeros(n)
    global_y = np.zeros(n)

    rows_per_process = (n + size - 1) // size
    start_row = rank * rows_per_process
    end_row = min((rank + 1) * rows_per_process, n)
    
    for i in range(n):
        comm.Bcast(global_y[:i], root=0)

        if start_row <= i < end_row:
            s = np.dot(decompose_matrix.l[i, :i], global_y[:i])
            y[i] = matrix.C[i] - s

        comm.Allreduce(y[i:i+1], global_y[i:i+1], op=MPI.SUM)
        
    matrix.X = np.zeros(n)

    for i in range(n - 1, -1, -1):
        local_x = np.zeros(1)
        if start_row <= i < end_row:
            s = np.dot(decompose_matrix.u[i, i+1:], matrix.X[i+1:])
            local_x[0] = (global_y[i] - s) / decompose_matrix.u[i, i]

        global_x = np.zeros(1)
        comm.Allreduce(local_x, global_x, op=MPI.SUM)
        matrix.X[i] = global_x[0]