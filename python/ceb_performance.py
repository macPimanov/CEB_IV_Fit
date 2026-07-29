"""
High-performance CEB computation library using C++ backend with threaded implementation.

This module provides a drop-in replacement for the compute_ceb_properties method
that uses a threaded C++ implementation for significant speedup.
"""

import os
from pathlib import Path

try:
    from .ceb_bindings import compute_ceb_properties_threaded, get_amp_constants
    HAS_CPP_BACKEND = True
except ImportError as e:
    print(f"Warning: C++ backend not available ({e}), falling back to pure Python")
    HAS_CPP_BACKEND = False


def is_cpp_backend_available():
    """Check if C++ threaded backend is available."""
    return HAS_CPP_BACKEND


def compute_ceb_properties(fitter_instance, use_cpp=True):
    """
    Compute CEB properties using the fastest available backend.
    
    Args:
        fitter_instance: IVParamFitter instance containing all parameters
        use_cpp: If True, prefer C++ backend; if False, use Python backend
        
    Returns:
        dictionary with 'Inum', 'Vnum', and 'time_spent' keys
    """
    if use_cpp and HAS_CPP_BACKEND:
        return _compute_cpp(fitter_instance)
    else:
        return _compute_python(fitter_instance)


def _compute_cpp(fitter_instance):
    """Use C++ threaded backend for computation."""
    # Prepare parameters dictionary for C++ library
    params = fitter_instance.par.copy()
    
    # Get amplifier noise parameters from the fitter instance
    ampnoise = {
        'voltage_noise': fitter_instance.amp_constants['voltage_noise'],
        'current_noise': fitter_instance.amp_constants['current_noise']
    }
    
    # Call C++ threaded function
    result = compute_ceb_properties_threaded(params, amp_noise=ampnoise)
    
    # Store results in the IVParamFitter instance
    fitter_instance.Inum = result['Inum']
    fitter_instance.Vnum = result['Vnum']
    
    return result


def _compute_python(fitter_instance):
    """Fallback to pure Python implementation."""
    import time
    start_time = time.time()
    
    # Call the original Python compute_ceb_properties method
    voltage_steps = fitter_instance._original_compute_ceb_properties()
    
    time_spent = time.time() - start_time
    
    return {
        'Inum': fitter_instance.Inum,
        'Vnum': fitter_instance.Vnum,
        'time_spent': time_spent,
        'error_code': 0
    }


def patch_fitter(fitter_instance):
    """
    Patch an IVParamFitter instance to use high-performance backend.
    
    Args:
        fitter_instance: IVParamFitter instance to patch
        
    Returns:
        The patched instance (same object)
    """
    if not HAS_CPP_BACKEND:
        print("C++ backend not available, using pure Python implementation")
        return fitter_instance
    
    # Store the original method
    if not hasattr(fitter_instance, '_original_compute_ceb_properties'):
        fitter_instance._original_compute_ceb_properties = fitter_instance.compute_ceb_properties
    
    # Replace the compute method
    fitter_instance.compute_ceb_properties = lambda: _compute_patch(fitter_instance)
    
    return fitter_instance


def _compute_patch(fitter_instance):
    """Patched compute method that dispatches to C++ backend."""
    result = _compute_cpp(fitter_instance)
    return len(result['Inum']) + 1  # Return voltage steps count