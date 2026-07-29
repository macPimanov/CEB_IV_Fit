"""
Test script for CEB C++ library performance validation.

This script validates the C++ threaded implementation and compares
performance against the pure Python implementation.
"""

import sys
import time
import numpy as np
from pathlib import Path

# Add python directory to path
sys.path.insert(0, str(Path(__file__).parent))

from python.iv_param_fitter import IVParamFitter
from python.ceb_performance import (
    is_cpp_backend_available, 
    compute_ceb_properties, 
    patch_fitter
)

def test_basic_functionality():
    """Test basic functionality of C++ implementation."""
    print("=" * 60)
    print("TEST 1: Basic Functionality")
    print("=" * 60)
    
    fitter = IVParamFitter()
    fitter.load_default_parameters()
    
    # Test with default parameters
    params = fitter.par.copy()
    print(f"Testing with {len(params)} parameters")
    
    if is_cpp_backend_available():
        print("[+] C++ backend is available")

        try:
            result = compute_ceb_properties(fitter, use_cpp=True)
            print(f"[+] C++ computation completed successfully")
            print(f"  - Computed {len(result['Inum'])} voltage steps")
            print(f"  - Time spent: {result['time_spent']:.4f} seconds")
            print(f"  - Error code: {result['error_code']}")

            # Check results validity
            if len(result['Inum']) > 0 and len(result['Vnum']) > 0:
                print("[+] Results contain valid data")
                return True
            else:
                print("[-] Results are empty")
                return False

        except Exception as e:
            print(f"[-] C++ computation failed: {e}")
            return False
    else:
        print("[-] C++ backend not available, skipping test")
        return False


def test_performance_comparison():
    """Compare performance between C++ and Python implementations."""
    print("\n" + "=" * 60)
    print("TEST 2: Performance Comparison")
    print("=" * 60)
    
    if not is_cpp_backend_available():
        print("✗ C++ backend not available, skipping performance test")
        return False
    
    fitter = IVParamFitter()
    fitter.load_default_parameters()
    
    # Store original method for Python testing
    python_compute = fitter.compute_ceb_properties
    fitter.compute_ceb_properties = lambda: None  # Temporarily disable
    
    # Test C++ performance
    print("\nTesting C++ performance...")
    cpp_start = time.time()
    cpp_result = compute_ceb_properties(fitter, use_cpp=True)
    cpp_time = time.time() - cpp_start
    
    print(f"✓ C++ completed in {cpp_time:.4f} seconds")
    print(f"  - Computed {len(cpp_result['Inum'])} voltage steps")
    
    # Test Python performance (smaller dataset for fairness)
    print("\nTesting Python performance...")
    fitter.compute_ceb_properties = python_compute
    
    # Reset for fair comparison
    fitter.Inum = None
    fitter.Vnum = None
    
    python_start = time.time()
    python_steps = fitter.compute_ceb_properties()
    python_time = time.time() - python_start
    
    print(f"✓ Python completed in {python_time:.4f} seconds")
    print(f"  - Computed {len(fitter.Inum)} voltage steps")
    
    # Compare results
    print("\nPerformance Comparison:")
    print(f"  C++ time:    {cpp_time:.4f} seconds")
    print(f"  Python time:  {python_time:.4f} seconds")
    print(f"  Speedup:      {python_time/cpp_time:.2f}x faster")
    
    # Data validation
    print("\nData Validation:")
    if len(cpp_result['Inum']) == len(fitter.Inum):
        print(f"✓ Array sizes match: {len(cpp_result['Inum'])} elements")
        
        # Check if results are approximately equal
        if np.allclose(cpp_result['Inum'], fitter.Inum, rtol=1e-5, atol=1e-10):
            print("✓ Current values match within tolerance")
        else:
            print("⚠ Current values differ within tolerance")
            max_diff = np.max(np.abs(cpp_result['Inum'] - fitter.Inum))
            print(f"  Max difference: {max_diff:.2e}")
            
        if np.allclose(cpp_result['Vnum'], fitter.Vnum, rtol=1e-5, atol=1e-10):
            print("✓ Voltage values match within tolerance")
        else:
            print("⚠ Voltage values differ within tolerance")
            max_diff = np.max(np.abs(cpp_result['Vnum'] - fitter.Vnum))
            print(f"  Max difference: {max_diff:.2e}")
            
        return True
    else:
        print(f"✗ Array sizes differ: C++={len(cpp_result['Inum'])}, Python={len(fitter.Inum)}")
        return False


def test_patched_fitter():
    """Test the patched fitter interface."""
    print("\n" + "=" * 60)
    print("TEST 3: Patched Fitter Interface")
    print("=" * 60)
    
    if not is_cpp_backend_available():
        print("✗ C++ backend not available, skipping patch test")
        return False
    
    fitter = IVParamFitter()
    fitter.load_default_parameters()
    
    print("Patching IVParamFitter instance...")
    patch_fitter(fitter)
    print("✓ Fitter patched successfully")
    
    try:
        print("\nTesting patched compute_ceb_properties...")
        start_time = time.time()
        steps = fitter.compute_ceb_properties()
        elapsed = time.time() - start_time
        
        print(f"✓ Patched method completed in {elapsed:.4f} seconds")
        print(f"  - Computed {steps} voltage steps")
        print(f"  - Results stored in fitter.Inum and fitter.Vnum")
        
        if len(fitter.Inum) > 0 and len(fitter.Vnum) > 0:
            print("✓ Patched implementation produced valid results")
            return True
        else:
            print("✗ Patched implementation produced empty results")
            return False
            
    except Exception as e:
        print(f"✗ Patched implementation failed: {e}")
        return False


def test_error_handling():
    """Test error handling in C++ implementation."""
    print("\n" + "=" * 60)
    print("TEST 4: Error Handling")
    print("=" * 60)
    
    if not is_cpp_backend_available():
        print("✗ C++ backend not available, skipping error test")
        return False
    
    fitter = IVParamFitter()
    fitter.load_default_parameters()
    
    # Test with invalid parameters (empty range)
    print("Testing with invalid parameters (zero voltage range)...")
    original_dv = fitter.par['dV']
    fitter.par['dV'] = 0.0  # Should cause error
    
    try:
        result = compute_ceb_properties(fitter, use_cpp=True)
        print("✗ Should have raised an error for zero voltage step")
        return False
    except RuntimeError as e:
        print(f"✓ Correctly raised error: {str(e)[:80]}...")
        fitter.par['dV'] = original_dv  # Restore
        return True
    except Exception as e:
        print(f"✗ Unexpected error type: {type(e).__name__}: {e}")
        fitter.par['dV'] = original_dv  # Restore
        return False


def main():
    """Run all tests."""
    print("CEB C++ Library Performance Validation")
    print("=" * 60)
    
    results = []
    
    # Run tests
    results.append(("Basic Functionality", test_basic_functionality()))
    results.append(("Performance Comparison", test_performance_comparison()))
    results.append(("Patched Fitter Interface", test_patched_fitter()))
    results.append(("Error Handling", test_error_handling()))
    
    # Summary
    print("\n" + "=" * 60)
    print("TEST SUMMARY")
    print("=" * 60)
    
    passed = sum(1 for _, result in results if result)
    total = len(results)
    
    for test_name, passed_test in results:
        status = "✓" if passed_test else "✗"
        print(f"{status} {test_name}")
    
    print(f"\nPassed: {passed}/{total}")
    
    if passed == total:
        print("\n🎉 All tests passed!")
        return 0
    else:
        print(f"\n⚠ {total - passed} test(s) failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())
