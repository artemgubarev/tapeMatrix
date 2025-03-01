import numpy as np
from mpi4py import MPI

def compare_numbers(file1, file2, tolerance=1e-5):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    data1_local = None
    data2_local = None
    local_mismatched_indices = []
    global_mismatched_indices = []
    all_match = True

    try:
        if rank == 0:
            try:
                data1 = np.loadtxt(file1)
                data2 = np.loadtxt(file2)
            except FileNotFoundError:
                print("One of the files was not found.")
                raise
            except ValueError:
                print("Error reading data. Make sure the file format is correct.")
                raise
            except Exception as e:
                print(f"An unexpected error occurred: {e}")
                raise

            if data1.shape != data2.shape:
                print("The files contain a different number of elements.")
                raise ValueError("Array sizes do not match")

            counts = [len(data1) // size + (1 if i < len(data1) % size else 0) for i in range(size)]
            displacements = [sum(counts[:i]) for i in range(size)]

        else:
            counts = None
            displacements = None

        counts = comm.bcast(counts, root=0)
        displacements = comm.bcast(displacements, root=0)

        if rank != 0:
            data1_local = np.empty(counts[rank], dtype=np.float64)
            data2_local = np.empty(counts[rank], dtype=np.float64)
        else:
            data1_local = np.empty(counts[rank], dtype=np.float64)
            data2_local = np.empty(counts[rank], dtype=np.float64)

        comm.Scatterv([data1, counts, displacements, MPI.DOUBLE], data1_local, root=0)
        comm.Scatterv([data2, counts, displacements, MPI.DOUBLE], data2_local, root=0)

        difference = np.abs(data1_local - data2_local)
        local_mismatched_indices = np.where(difference > tolerance)[0] + displacements[rank]

        all_local_mismatched_indices = comm.gather(local_mismatched_indices, root=0)

        if rank == 0:
            for sublist in all_local_mismatched_indices:
                global_mismatched_indices.extend(sublist)
            all_match = len(global_mismatched_indices) == 0

    except Exception as e:
        if rank == 0:
            print(f"Error in rank 0: {e}")
        all_match = False
        global_mismatched_indices = []

    all_match = comm.bcast(all_match, root=0)

    if rank == 0:
        if all_match:
            print("Test correct")
        else:
            print("Test failed")
