def compare_files(file1_path, file2_path, tolerance=1e-5):
    try:
        with open(file1_path, 'r') as f1:
            numbers1 = [float(num) for num in f1.read().split()]
        
        with open(file2_path, 'r') as f2:
            numbers2 = [float(num) for num in f2.read().split()]
        
        if len(numbers1) != len(numbers2):
            print(f"Failed diff number: {len(numbers1)} vs {len(numbers2)}")
            return False
        
        for i, (num1, num2) in enumerate(zip(numbers1, numbers2)):
            if abs(num1 - num2) > tolerance:
                print(f"Dismatch {i}: {num1} vs {num2}")
                return False
        
        return True
    
    except Exception as e:
        print(f"Failed: {e}")
        return False
    
file1_path = "solution.txt"
file2_path = "mSolutions/msolution100x60.txt"

if compare_files(file1_path, file2_path):
    print("Test correct")
else:
    print("Test failed")
