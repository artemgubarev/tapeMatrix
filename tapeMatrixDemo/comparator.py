from mpi4py import MPI

def compare_numbers(file1, file2, tolerance=1e-5):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    if rank == 0:
        with open(file1, 'r') as f1:
            numbers1 = list(map(float, f1.read().split()))
        with open(file2, 'r') as f2:
            numbers2 = list(map(float, f2.read().split()))
        if len(numbers1) != len(numbers2):
            print("test failed")
            return
        for i, (num1, num2) in enumerate(zip(numbers1, numbers2)):
            if abs(num1 - num2) > tolerance:
                print(f"test failed {i}: {num1} != {num2}")
                return
        print("test correct")
    else:
        pass
    comm.Barrier()
