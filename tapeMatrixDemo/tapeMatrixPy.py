from mpi4py import MPI
import time
import sys
import os
sys.path.append("../tapeMatrix")
import matrix as mtrx
import solver_mpi as solver
import comparator as comp
import numpy as np
import colorama

def get_output_filename(input_filename):
    base, ext = os.path.splitext(input_filename)
    return f"{base}{ext}"

def compare_numbers(file1_path, file2_path, tolerance=1e-5):
    colorama.init()
    try:
        with open(file1_path, 'r') as file1, open(file2_path, 'r') as file2:
            numbers1 = [float(num) for num in file1.read().split()]
            numbers2 = [float(num) for num in file2.read().split()]
            if len(numbers1) != len(numbers2):
                print(colorama.Fore.RED + "TestFailed: Files have different number of numbers." + colorama.Style.RESET_ALL)
                return
            for i, (num1, num2) in enumerate(zip(numbers1, numbers2)):
                if abs(num1 - num2) <= tolerance:
                    print(colorama.Fore.GREEN + f"TestCorrect: Number {i+1} - {num1:.5f} ≈ {num2:.5f}" + colorama.Style.RESET_ALL)
                else:
                    print(colorama.Fore.RED + f"TestFailed: Number {i+1} - {num1:.5f} != {num2:.5f} (Difference: {abs(num1 - num2):.5f})" + colorama.Style.RESET_ALL)
    except FileNotFoundError:
        print(colorama.Fore.RED + "TestFailed: One or both files not found." + colorama.Style.RESET_ALL)
    except ValueError:
        print(colorama.Fore.RED + "TestFailed: Invalid number format in one of the files." + colorama.Style.RESET_ALL)

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
    np.savetxt("solution.txt", [matrix.X], fmt='%.6f', delimiter=' ')
    #print(f"Solution X: {matrix.X}")
    print(f"Elapsed time: {elapsed_time:.6f} seconds")
    compare_numbers("solution.txt", get_output_filename(filename))