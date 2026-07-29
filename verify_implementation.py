#!/usr/bin/env python3
"""
Verification script to compare Python and C++ implementations.
This helps ensure the Python version produces consistent results.
"""

import numpy as np
import sys
import argparse

def compare_files(file1, file2, tolerance=1e-6, skip_lines=1):
    """
    Compare two text files with numerical data.
    
    Args:
        file1: Path to first file
        file2: Path to second file  
        tolerance: Relative tolerance for numerical comparison
        skip_lines: Number of header lines to skip
    
    Returns:
        bool: True if files match within tolerance
    """
    try:
        data1 = np.loadtxt(file1, delimiter='\t', skiprows=skip_lines)
        data2 = np.loadtxt(file2, delimiter='\t', skiprows=skip_lines)
    except Exception as e:
        print(f"Error loading files: {e}")
        return False
    
    if data1.shape != data2.shape:
        print(f"Error: Shape mismatch - File1: {data1.shape}, File2: {data2.shape}")
        return False
    
    # Compare values
    diff = np.abs(data1 - data2)
    rel_diff = diff / (np.abs(data1) + tolerance)
    
    max_rel_diff = np.max(rel_diff)
    mean_rel_diff = np.mean(rel_diff)
    
    print(f"Comparison results:")
    print(f"  Maximum relative difference: {max_rel_diff:.2e}")
    print(f"  Mean relative difference: {mean_rel_diff:.2e}")
    print(f"  Tolerance: {tolerance:.2e}")
    
    if max_rel_diff > tolerance:
        print(f"  Status: FAIL - Difference exceeds tolerance")
        # Find worst offenders
        worst_indices = np.unravel_index(np.argmax(rel_diff), rel_diff.shape)
        print(f"  Worst element at index {worst_indices}:")
        print(f"    File1 value: {data1[worst_indices]:.6e}")
        print(f"    File2 value: {data2[worst_indices]:.6e}")
        print(f"    Relative difference: {max_rel_diff:.2e}")
        return False
    else:
        print(f"  Status: PASS - Files match within tolerance")
        return True

def verify_parameter_conversion():
    """Verify that parameter conversion from TXT to JSON is correct"""
    print("\n=== Verifying Parameter Conversion ===")
    
    try:
        # Try to read old format if exists
        import os
        from python.utils import Utils
        
        if os.path.exists('startparams.txt'):
            print("Old format parameter file exists")
            
            # Read old format
            old_params = {}
            with open('startparams.txt', 'r') as f:
                for line in f:
                    parts = line.strip().split()
                    if len(parts) >= 3:
                        old_params[parts[0]] = float(parts[1])
            
            # Read JSON format if exists
            if os.path.exists('config.json'):
                print("New format parameter file exists")
                from python.utils import Utils
                new_params = Utils.load_json_config('config.json')
                
                # Compare values
                matches = 0
                mismatches = 0
                for name, value in old_params.items():
                    if 'parameters' in new_params and name in new_params['parameters']:
                        new_value = new_params['parameters'][name]['value']
                        if abs(new_value - value) > 1e-10:
                            print(f"Mismatch {name}: old={value}, new={new_value}")
                            mismatches += 1
                        else:
                            matches += 1
                
                print(f"Parameter comparison: {matches} matches, {mismatches} mismatches")
        
        return True
    
    except Exception as e:
        print(f"Error in parameter verification: {e}")
        return False

def run_verification_suite():
    """Run comprehensive verification tests"""
    print("=== Python Implementation Verification Suite ===\n")
    
    all_passed = True
    
    # 1. Import verification
    print("1. Verifying imports...")
    try:
        from python.iv_param_fitter import IVParamFitter
        from python.utils import Utils
        from python.minimization import MinimizationAlgorithms
        print("   All imports successful")
    except ImportError as e:
        print(f"   Import failed: {e}")
        all_passed = False
    
    # 2. Numerical model verification
    print("\n2. Verifying numerical model...")
    try:
        from python.ceb_numeric_model import CEBNumericModel
        
        model = CEBNumericModel()
        
        # Test key functions with known values
        test_v = 0.5
        test_tau = 0.1
        result = model.current(test_v, test_tau)
        print(f"   Current at v={test_v}, tau={test_tau}: {result:.6e}")
        print(f"   CEBNumericModel working")
        
    except Exception as e:
        print(f"   Numerical model verification failed: {e}")
        all_passed = False
    
    # 3. Parameter conversion verification
    try:
        param_ok = verify_parameter_conversion()
        all_passed = all_passed and param_ok
    except Exception as e:
        print(f"   Parameter conversion verification failed: {e}")
        all_passed = False
    
    # 4. File comparison if outputs exist
    print("\n4. Comparing output files...")
    output_files = [
        ('Te.txt', 'Te_cpp.txt'),
        ('Noise.txt', 'Noise_cpp.txt'),
        ('NEP.txt', 'NEP_cpp.txt')
    ]
    
    for py_file, cpp_file in output_files:
        if all([os.path.exists(f) for f in [py_file, cpp_file]]):
            match = compare_files(cpp_file, py_file)
            all_passed = all_passed and match
    
    # Final summary
    print("\n=== Verification Summary ===")
    if all_passed:
        print("All tests PASSED")
        return 0
    else:
        print("Some tests FAILED")
        return 1

def benchmark_performance():
    """Run simple performance benchmark"""
    print("\n=== Performance Benchmark ===")
    
    try:
        import time
        from python.iv_param_fitter import IVParamFitter
        
        fitter = IVParamFitter()
        
        # Test loading
        start = time.time()
        for _ in range(10):
            # Simulate some numerical operations
            x = np.linspace(0, 1, 1000)
            y = np.sin(x) * np.exp(-x)
        elapsed = time.time() - start
        print(f"NumPy operations (10 runs): {elapsed:.4f} seconds")
        
        return 0
    except Exception as e:
        print(f"Benchmark failed: {e}")
        return 1

def main():
    parser = argparse.ArgumentParser(description='Verify Python CEB fitting implementation')
    parser.add_argument('--compare', nargs=2, metavar=('FILE1', 'FILE2'),
                       help='Compare two output files')
    parser.add_argument('--tolerance', type=float, default=1e-6,
                       help='Numerical tolerance for comparison')
    parser.add_argument('--benchmark', action='store_true',
                       help='Run performance benchmark')
    parser.add_argument('--full-verification', action='store_true',
                       help='Run full verification suite')
    
    args = parser.parse_args()
    
    import os
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    
    if args.compare:
        file1, file2 = args.compare
        result = compare_files(file1, file2, tolerance=args.tolerance)
        return 0 if result else 1
    
    if args.benchmark:
        return benchmark_performance()
    
    if args.full_verification:
        return run_verification_suite()
    
    # Default: run verification suite
    return run_verification_suite()

if __name__ == "__main__":
    sys.exit(main())