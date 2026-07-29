#!/usr/bin/env python3
"""
Simple test script for the Python CEB fitting implementation.
This demonstrates the basic functionality without requiring full experimental data.
"""

import numpy as np
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def test_utils():
    """Test utility functions"""
    print("Testing Utils class...")
    from python.utils import Utils
    
    # Test data loading
    test_data_file = "test_data.txt"
    
    # Create simple test data
    test_V = np.linspace(0, 1e-3, 100)
    test_I = test_V * 1000  # Simple resistance of 1 ohm
    test_data = np.column_stack((test_V, test_I))
    np.savetxt(test_data_file, test_data, delimiter='\t')
    
    try:
        Iexp, Vexp = Utils.load_experimental_data(test_data_file)
        print(f"  Loaded {len(Iexp)} data points")
        
        # Test resampling
        Inum = Iexp * 0.95  # Simulated numerical model
        Vnum = Vexp
        
        Irex, Vrex = Utils.resample(Iexp, Vexp, Inum, Vnum)
        print(f"  Resampled {len(Irex)} points")
        
        # Test chi_sq calculation
        chi2 = Utils.chi_sq(Inum, Irex)
        print(f"  Chi-squared value: {chi2:.6e}")
        
        # Test golden minimization
        def test_func(x):
            return (x - 0.5) ** 2
        
        x_min = Utils.golden_minimize(test_func, 0.0, 1.0, 0.5)
        print(f"  Golden minimization test: found minimum at {x_min:.6f} (expected 0.5)")
        
        return True
    finally:
        # Clean up
        if os.path.exists(test_data_file):
            os.remove(test_data_file)

def test_numerical_model():
    """Test CEB numerical model"""
    print("\nTesting CEB Mathematical Model...")
    from python.ceb_numeric_model import CEBNumericModel
    
    model = CEBNumericModel()
    
    # Test current calculation
    v_test = 0.5
    tau_test = 0.1
    current_val = model.current(v_test, tau_test)
    print(f"  Current at v={v_test}, tau={tau_test}: {current_val:.6e}")
    
    # Test current integral
    DT_test = 0.5
    tauE_test = 0.2
    current_int_val = model.current_integral(DT_test, v_test, tau_test, tauE_test)
    print(f"  Current integral: {current_int_val:.6e}")
    
    # Test cooling power
    power_cool_val = model.power_cool(v_test, tau_test, tauE_test)
    print(f"  Cooling power: {power_cool_val:.6e}")
    
    # Test Andreev current
    Wt_test = 0.0001
    tm_test = 1.0
    and_current_val = model.and_current(DT_test, v_test, tauE_test, Wt_test, tm_test)
    print(f"  Andreev current: {and_current_val:.6e}")
    
    # Test power cooling integral
    po, ps = model.power_cool_integral(DT_test, v_test, tau_test, tauE_test)
    print(f"  Power cooling integral: po={po:.6e}, ps={ps:.6e}")
    
    return True

def test_minimization():
    """Test minimization algorithms"""
    print("\nTesting Minimization Algorithms...")
    from python.minimization import MinimizationAlgorithms
    
    # Test simple quadratic function
    def quadratic(x):
        return (x - 0.5) * (x - 0.5)
    
    # Golden section minimization
    x_golden, f_golden = MinimizationAlgorithms.golden_section_minimization(quadratic, 0.0, 1.0)
    print(f"  Golden section: minimum at {x_golden:.6f}, value {f_golden:.6e}")
    
    # Brent minimization
    x_brent, f_brent = MinimizationAlgorithms.brent_minimization(quadratic, 0.0, 1.0)
    print(f"  Brent minimization: minimum at {x_brent:.6f}, value {f_brent:.6e}")
    
    # Test with scipy minimize
    def multi_dim(x):
        return (x[0] - 0.3) ** 2 + (x[1] - 0.7) ** 2
    
    x_scipt, f_scipy = MinimizationAlgorithms.scipy_minimize(multi_dim, [0.0, 0.0])
    print(f"  SciPy minimize: minimum at {x_scipy}, value {f_scipy:.6e}")
    
    return True

def test_config_handling():
    """Test JSON configuration handling"""
    print("\nTesting Configuration Handling...")
    from python.utils import Utils
    
    # Create sample configuration
    sample_config = Utils.create_sample_config()
    config_file = "test_config.json"
    
    try:
        # Save configuration
        Utils.save_json_config(config_file, sample_config)
        print(f"  Saved sample configuration to {config_file}")
        
        # Load configuration
        loaded_config = Utils.load_json_config(config_file)
        print(f"  Loaded configuration with {len(loaded_config['parameters'])} parameters")
        
        # Test IVParamFitter with config
        from python.iv_param_fitter import IVParamFitter
        fitter = IVParamFitter(config_file=config_file)
        print(f"  Created IVParamFitter with configuration")
        print(f"  Data file: {fitter.data_file}")
        print(f"  Parameters to fit: {sum(fitter.to_fit.values())}")
        
        return True
    finally:
        # Clean up
        if os.path.exists(config_file):
            os.remove(config_file)

def test_complete_workflow():
    """Test complete workflow with minimal setup"""
    print("\nTesting Complete Workflow...")
    from python.iv_param_fitter import IVParamFitter
    
    # Create minimal test data
    test_data_file = "test_iv_data.txt"
    test_V = np.linspace(0, 2e-6, 50)
    test_I = np.where(test_V < 1e-6, test_V * 100, test_V * 50 + 50e-6) / 1e9  # Convert to Amps
    
    try:
        test_data = np.column_stack((test_V, test_I))
        np.savetxt(test_data_file, test_data, delimiter='\t')
        
        # Create fitter
        fitter = IVParamFitter()  # Uses default parameters
        fitter.data_file = test_data_file
        
        # Load data
        fitter.load_experiment_data(test_data_file)
        print(f"  Loaded {len(fitter.Iexp)} experimental data points")
        
        # Compute CEB properties (this will take some time)
        print("  Computing CEB properties (may take 30-60 seconds)...")
        points_computed = fitter.compute_ceb_properties()
        print(f"  Computed {points_computed} CEB model points")
        
        # Single fit iteration to test fitting
        print("  Testing single parameter fit...")
        fitter.par['beta'] = 0.1  # Reset to different value
        
        def objective(param_value):
            return fitter(param_value, 'beta')
        
        from python.minimization import MinimizationAlgorithms
        x_min, f_min = MinimizationAlgorithms.golden_section_minimization(
            objective, 0.05 * fitter.par['beta'], 2.0 * fitter.par['beta']
        )
        
        print(f"  Fitted beta: {x_min:.6e}, chi-square: {f_min:.6e}")
        
        return True
        
    except Exception as e:
        print(f"  Complete workflow test failed: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        # Clean up
        for file in [test_data_file, 'Te.txt', 'Noise.txt', 'NEP.txt', 'G.txt', 
                    'fitparameters_new.txt', 'converg.txt']:
            if os.path.exists(file):
                os.remove(file)

def main():
    print("=== Testing Python CEB Implementation ===\n")
    
    tests = [
        ("Utility Functions", test_utils),
        ("Numerical Model", test_numerical_model),
        ("Minimization Algorithms", test_minimization),
        ("Configuration Handling", test_config_handling),
    ]
    
    # Run quick tests
    results = []
    for test_name, test_func in tests:
        try:
            result = test_func()
            results.append((test_name, result))
        except Exception as e:
            print(f"FAILED: {test_name} - {e}")
            results.append((test_name, False))
    
    # Ask about complete workflow test
    print("\nNote: Complete workflow test will take 30-60 seconds")
    print("      (computes full CEB physics model)")
    
    results.append(("Complete Workflow", test_complete_workflow()))
    
    # Summary
    print("\n=== Test Results Summary ===")
    passed = sum(1 for _, result in results if result)
    total = len(results)
    
    for test_name, result in results:
        status = "PASSED" if result else "FAILED"
        print(f"  {test_name}: {status}")
    
    print(f"\nTotal: {passed}/{total} tests passed")
    return 0 if passed == total else 1

if __name__ == "__main__":
    sys.exit(main())