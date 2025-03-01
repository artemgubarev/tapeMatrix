import numpy as np

MAX_NUMBERS = 1000000

def load_numbers(filename):
    try:
        with open(filename, "r") as file:
            numbers = []
            count = 0
            for line in file:
                try:
                    number = float(line.strip())
                    numbers.append(number)
                    count += 1
                    if count >= MAX_NUMBERS:
                        break
                except ValueError:
                    print(f"Warning: Invalid number format in {filename} - skipping line: {line.strip()}")
                    continue
            return np.array(numbers), count
    except FileNotFoundError:
        print(f"Error: Cannot open file {filename}")
        return None, 0
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return None, 0

def compare_numbers(vec1, vec2, count1, count2, epsilon=1e-5):
    if count1 != count2:
        print("Error: Files contain different number of values!")
        return False
    diff = np.abs(vec1 - vec2)
    mismatches = diff > epsilon
    if np.any(mismatches):
        indices = np.where(mismatches)[0]
        for i in indices:
            print(f"Mismatch at index {i}: {vec1[i]} vs {vec2[i]} (diff = {diff[i]})")
        return False
    return True
