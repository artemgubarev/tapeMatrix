import numpy as np
from mpi4py import MPI

def compare_numbers(file1, file2, tolerance=1e-5):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    all_match = True

    try:
        if rank == 0:
            try:
                with open(file1, 'r') as f1, open(file2, 'r') as f2:
                    data1 = [np.fromstring(line, dtype=float, sep=' ') for line in f1]
                    data2 = [np.fromstring(line, dtype=float, sep=' ') for line in f2]
                    
                max_len = max(max(len(row) for row in data1), max(len(row) for row in data2))
                
                data1 = np.array([np.pad(row, (0, max_len - len(row)), 'constant', constant_values=np.nan) for row in data1])
                data2 = np.array([np.pad(row, (0, max_len - len(row)), 'constant', constant_values=np.nan) for row in data2])
                
                if data1.shape[0] != data2.shape[0]:
                    print("The files contain a different number of rows.")
                    raise ValueError("Row counts do not match")


            except FileNotFoundError:
                print("One of the files was not found.")
                raise
            except ValueError as e:
                print(f"Error reading data: {e}")
                raise
            except Exception as e:
                print(f"An unexpected error occurred: {e}")
                raise

            counts = [len(data1) // size + (1 if i < len(data1) % size else 0) for i in range(size)]
            displacements = [sum(counts[:i]) for i in range(size)]

        else:
            counts = None
            displacements = None
            max_len = None

        counts = comm.bcast(counts, root=0)
        displacements = comm.bcast(displacements, root=0)
        max_len = comm.bcast(max_len if rank == 0 else None, root=0)

        if rank != 0:
            data1_local = np.empty((counts[rank], max_len), dtype=np.float64)
            data2_local = np.empty((counts[rank], max_len), dtype=np.float64)
        else:
            data1_local = np.empty((counts[rank], max_len), dtype=np.float64)
            data2_local = np.empty((counts[rank], max_len), dtype=np.float64)

        comm.Scatterv([data1, counts * max_len, [d * max_len for d in displacements], MPI.DOUBLE], data1_local, root=0)
        comm.Scatterv([data2, counts * max_len, [d * max_len for d in displacements], MPI.DOUBLE], data2_local, root=0)
        
        difference = np.abs(data1_local - data2_local)
        local_mismatched = (difference > tolerance) & ~ (np.isnan(data1_local) & np.isnan(data2_local))
        local_mismatched_indices = np.where(local_mismatched)[0] + displacements[rank]
        

        all_local_mismatched_indices = comm.gather(local_mismatched_indices, root=0)


        if rank == 0:
             global_mismatched_indices = []
             for sublist in all_local_mismatched_indices:
                 global_mismatched_indices.extend(sublist)
             all_match = len(global_mismatched_indices) == 0


    except Exception as e:
        if rank == 0:
            print(f"Error in rank 0: {e}")
        all_match = False

    all_match = comm.bcast(all_match, root=0)
    if rank == 0:
      if all_match:
        print("Test correct")
      else:
        print("Test failed")
