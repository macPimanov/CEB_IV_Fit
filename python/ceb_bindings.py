"""
ctypes bindings for the CEB (Cooling Electron Bolometer) C++ library.

This module provides a high-performance Python interface to the threaded CEB computation
library, which significantly speeds up the compute_ceb_properties calculation.
"""

import ctypes
import os
import numpy as np
from pathlib import Path

# Get the library path
current_dir = Path(__file__).parent
lib_name = "libceb_library"
lib_path = None

# Try different library extensions based on platform
if os.name == 'nt':  # Windows
    lib_path = current_dir / f"{lib_name}.dll"
elif os.name == 'posix':
    lib_path = current_dir / f"{lib_name}.so"
else:
    raise ImportError(f"Unsupported platform: {os.name}")

if not lib_path.exists():
    raise ImportError(f"CEB library not found at: {lib_path}")

# Load the library using WinDLL on Windows
if os.name == 'nt':
    lib = ctypes.WinDLL(str(lib_path.absolute()))
else:
    lib = ctypes.CDLL(str(lib_path.absolute()))


# Define C structures matching the C library
class DoubleArray(ctypes.Structure):
    _fields_ = [
        ("data", ctypes.POINTER(ctypes.c_double)),
        ("array_size", ctypes.c_size_t)
    ]


class CEBParameters(ctypes.Structure):
    _fields_ = [
        # Physical parameters
        ("M", ctypes.c_double),          # number of bolometers in series
        ("MP", ctypes.c_double),         # number of bolometers in parallel
        ("Pbg", ctypes.c_double),        # incoming power [pW]
        ("beta", ctypes.c_double),       # returning power ratio
        ("TephPOW", ctypes.c_double),    # exponent for Te-ph, 7, 6, or 5
        ("Vol", ctypes.c_double),        # volume of absorber [um³]
        ("Z", ctypes.c_double),         # heat exchange in normal metal
        ("Tc", ctypes.c_double),         # critical temperature [K]
        ("Rn", ctypes.c_double),        # normal resistance for 1 bolometer [Ohm]
        ("Rleak", ctypes.c_double),     # leakage resistance per 1 bolometer [Ohm]
        ("Wt", ctypes.c_double),        # transparency of barrier
        ("tm", ctypes.c_double),        # depairing energy
        ("ii", ctypes.c_double),        # coefficient for Andreev current
        ("Ra", ctypes.c_double),        # normal resistance of 1 absorber [Ohm]
        ("Tp", ctypes.c_double),        # phonon temperature [K]
        ("F", ctypes.c_double),          # main frequency [GHz]
        ("dF", ctypes.c_double),        # bandwidth [GHz]
        ("dVFinVg", ctypes.c_double),   # voltage range end [Vg units]
        ("dVStartVg", ctypes.c_double), # voltage range start [Vg units]
        ("dV", ctypes.c_double),        # voltage step [V]
        ("voltage_noise", ctypes.c_double),  # [V/sqrt(Hz)]
        ("current_noise", ctypes.c_double),  # [A/sqrt(Hz)]
    ]


class CEBResult(ctypes.Structure):
    _fields_ = [
        ("Inum", DoubleArray),  # Numerical current values
        ("Vnum", DoubleArray),  # Numerical voltage values
        ("time_spent", ctypes.c_double), # Computation time in seconds
        ("error_code", ctypes.c_int),    # 0 for success, non-zero for error
        ("error_message", ctypes.c_char * 256)
    ]


# Set up function signature for compute_ceb_properties_threaded
lib.compute_ceb_properties_threaded.argtypes = [
    ctypes.POINTER(CEBParameters),
    ctypes.POINTER(CEBResult)
]
lib.compute_ceb_properties_threaded.restype = None

# Set up function signature for free_ceb_result
lib.free_ceb_result.argtypes = [ctypes.POINTER(CEBResult)]
lib.free_ceb_result.restype = None


def compute_ceb_properties_threaded(params_dict, amp_noise=None):
    """
    Compute CEB properties using the threaded C++ library.
    
    Args:
        params_dict: Dictionary containing physical parameters
        amp_noise: Dictionary with 'voltage_noise' and 'current_noise' keys (optional)
        
    Returns:
        Dictionary containing:
            - 'Inum': numpy array of numerical current values  
            - 'Vnum': numpy array of numerical voltage values
            - 'time_spent': computation time in seconds
            - 'error_code': 0 for success
            - 'error_message': error message if error occurred
        
    Raises:
        RuntimeError: If the C++ function returns an error
    """
    # Default amplifier noise (AD745)
    if amp_noise is None:
        amp_noise = {
            'voltage_noise': 1.0e-8,    # V/sqrt(Hz)
            'current_noise': 6.9e-15     # A/sqrt(Hz)
        }
    
    # Create CEBParameters structure from dictionary
    ceb_params = CEBParameters(
        M=float(params_dict.get('M', 2)),
        MP=float(params_dict.get('MP', 1)),
        Pbg=float(params_dict.get('Pbg', 0.0)),
        beta=float(params_dict.get('beta', 0.111)),
        TephPOW=float(params_dict.get('TephPOW', 5.0)),
        Vol=float(params_dict.get('Vol', 0.02)),
        Z=float(params_dict.get('Z', 0.5)),
        Tc=float(params_dict.get('Tc', 1.18)),
        Rn=float(params_dict.get('Rn', 11500.0)),
        Rleak=float(params_dict.get('Rleak', 40000000.0)),
        Wt=float(params_dict.get('Wt', 0.0001)),
        tm=float(params_dict.get('tm', 1.0)),
        ii=float(params_dict.get('ii', 0.0)),
        Ra=float(params_dict.get('Ra', 200.0)),
        Tp=float(params_dict.get('Tp', 0.19)),
        F=float(params_dict.get('F', 14.2)),
        dF=float(params_dict.get('dF', 0.1)),
        dVFinVg=float(params_dict.get('dVFinVg', 1.1)),
        dVStartVg=float(params_dict.get('dVStartVg', 0.0)),
        dV=float(params_dict.get('dV', 2e-6)),
        voltage_noise=float(amp_noise['voltage_noise']),
        current_noise=float(amp_noise['current_noise'])
    )
    
    # Create CEBResult structure
    result = CEBResult()
    
    # Call the C function
    lib.compute_ceb_properties_threaded(ctypes.byref(ceb_params), ctypes.byref(result))
    
    # Check for errors
    if result.error_code != 0:
        error_msg = result.error_message.decode('utf-8').rstrip('\x00')
        lib.free_ceb_result(ctypes.byref(result))
        raise RuntimeError(f"CEB computation failed (code {result.error_code}): {error_msg}")
    
    # Convert C arrays to numpy arrays
    inum_array = np.ctypeslib.as_array(
        ctypes.cast(result.Inum.data, ctypes.POINTER(ctypes.c_double)),
        shape=(result.Inum.array_size,)
    ).copy()
    
    vnum_array = np.ctypeslib.as_array(
        ctypes.cast(result.Vnum.data, ctypes.POINTER(ctypes.c_double)),
        shape=(result.Vnum.array_size,)
    ).copy()
    
    # Free C side memory
    lib.free_ceb_result(ctypes.byref(result))
    
    # Return results as Python dictionary
    return {
        'Inum': inum_array,
        'Vnum': vnum_array,
        'time_spent': result.time_spent,
        'error_code': result.error_code
    }


def get_amp_constants(amp_type='AD745'):
    """
    Get amplifier noise constants for different amplifier types.
    
    Args:
        amp_type: String identifier for amplifier type
        
    Returns:
        Dictionary with 'voltage_noise' and 'current_noise' keys
    """
    amplifiers = {
        'AD745': {'voltage_noise': 1.0e-8, 'current_noise': 6.9e-15},
        'OPA111': {'voltage_noise': 8.0e-9, 'current_noise': 0.8e-15},
        'AD797': {'voltage_noise': 0.9e-9, 'current_noise': 2000.0e-15},
        'IFN146': {'voltage_noise': 1.1e-9, 'current_noise': 0.3e-15},
        'OPA1641': {'voltage_noise': 5.1e-9, 'current_noise': 0.8e-15},
    }
    return amplifiers.get(amp_type, amplifiers['AD745'])