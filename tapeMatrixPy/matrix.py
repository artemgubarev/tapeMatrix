import numpy as np
from mpi4py import MPI

class Matrix:
    def __init__(self, n=0, b=0, C=None, A=None, X=None):
        self.n = n
        self.b = b
        self.C = C
        self.A = A
        self.X = X

class DecomposeMatrix:
    def __init__(self, l=None, u=None):
        self.l = l
        self.u = u

def read_matrix(filename, rank):
    comm = MPI.COMM_WORLD
    try:
        with open(filename, "r") as file:
            # Read n and b from separate lines
            try:
                n = int(file.readline().strip())
                b = int(file.readline().strip())
            except ValueError:
                if rank == 0:
                    print(f"Error: 'n' and 'b' must be integers in '{filename}'.")
                comm.Abort(1)
            except Exception as e:
                if rank == 0:
                    print(f"Error reading n and b from '{filename}': {e}")
                comm.Abort(1)

            A = np.zeros((n, n))
            C = np.zeros(n)

            # Read matrix A
            for i in range(n):
                try:
                    line = file.readline().strip()
                    row = list(map(float, line.split()))
                    if len(row) != n:
                        if rank == 0:
                            print(f"Error: Row {i+1} in '{filename}' has {len(row)} elements (expected {n}).")
                        comm.Abort(1)
                    A[i, :] = row
                except ValueError:
                    if rank == 0:
                        print(f"Error: Invalid numeric value in row {i+1} of '{filename}'.")
                    comm.Abort(1)
                except Exception as e:
                    if rank == 0:
                        print(f"Error reading row {i+1} from '{filename}': {e}")
                    comm.Abort(1)

            # Read vector C
            try:
                line = file.readline().strip()
                C[:] = list(map(float, line.split()))
                if len(C) != n:
                    if rank == 0:
                        print(f"Error: Vector C in '{filename}' has incorrect number of elements (expected {n}).")
                    comm.Abort(1)
            except ValueError:
                if rank == 0:
                    print(f"Error: Invalid numeric value (or not enough values) for vector C in '{filename}'.")
                comm.Abort(1)
            except Exception as e:
                if rank == 0:
                    print(f"Error reading vector C from '{filename}': {e}")
                comm.Abort(1)

        return Matrix(n=n, b=b, C=C, A=A)

    except FileNotFoundError:
        if rank == 0:
            print(f"Error: File '{filename}' not found.")
        comm.Abort(1)
    except Exception as e:
        if rank == 0:
            print(f"An unexpected error occurred: {e}")
        comm.Abort(1)

    except FileNotFoundError:
        if rank == 0:
            print(f"Error: File '{filename}' not found.")
        comm.Abort(1)
    except Exception as e:
        if rank == 0:
           print(f"An unexpected error occurred while reading the file: {e}")
        comm.Abort(1)