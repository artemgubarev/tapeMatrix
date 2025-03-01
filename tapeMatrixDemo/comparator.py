from mpi4py import MPI
import numpy as np

def compare_numbers(file1, file2, tolerance=1e-5):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    if rank == 0:
        numbers1 = np.loadtxt(file1)
        numbers2 = np.loadtxt(file2)
        
        if numbers1.shape != numbers2.shape:
            print("test failed: количество чисел в файлах различается")
            return
        
        if np.allclose(numbers1, numbers2, atol=tolerance):
            print("test correct")
        else:
            print("test failed: числа в файлах различаются")
    else:
        pass
    
    comm.Barrier()
