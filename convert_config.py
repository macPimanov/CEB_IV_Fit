#!/usr/bin/env python3
"""
Convert old startparams.txt format to new JSON configuration format.
This helps users migrate from the C++ implementation to Python version.
"""

import sys
import json
import argparse
from pathlib import Path

def convert_txt_to_json(txt_file, json_output_file):
    """
    Convert old startparams.txt format to JSON configuration.
    
    Old format (startparams.txt):
        parameter_name value fit_flag
        e.g., "beta 0.111 1"
    
    New format (JSON):
        {
            "parameters": {
                "beta": {"value": 0.111, "vary": true, "min": ..., "max": ...}
            }
        }
    """
    
    # Default parameter bounds
    default_bounds = {
        'Pbg': {'min': 0.0, 'max': 1.0},
        'beta': {'min': 0.0, 'max': 1.0},
        'TephPOW': {'min': 1.0, 'max': 7.0},
        'Vol': {'min': 0.001, 'max': 0.1},
        'Z': {'min': 0.1, 'max': 2.0},
        'Tc': {'min': 1.0, 'max': 1.5},
        'Rn': {'min': 5000.0, 'max': 20000.0},
        'Rleak': {'min': 1e6, 'max': 1e8},
        'Wt': {'min': 0.0, 'max': 0.001},
        'tm': {'min': 0.5, 'max': 2.0},
        'ii': {'min': 0.0, 'max': 1.0},
        'Ra': {'min': 50.0, 'max': 1000.0},
        'M': {'min': 1, 'max': 10},
        'MP': {'min': 1, 'max': 10},
        'Tp': {'min': 0.1, 'max': 0.4},
        'F': {'min': 1.0, 'max': 100.0},
        'dF': {'min': 0.01, 'max': 1.0},
        'dVFinVg': {'min': 0.5, 'max': 2.0},
        'dVStartVg': {'min': 0.0, 'max': 0.5},
        'dV': {'min': 1e-6, 'max': 10e-6}
    }
    
    # Read and parse the old format
    parameters = {}
    
    with open(txt_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            parts = line.split()
            if len(parts) < 3:
                print(f"Warning: Line {line_num} has insufficient data, skipping: {line}")
                continue
            
            param_name = parts[0]
            try:
                param_value = float(parts[1])
                fit_flag = int(parts[2])
                vary = bool(fit_flag)
            except (ValueError, IndexError) as e:
                print(f"Warning: Line {line_num} has invalid data format, skipping: {line}")
                continue
            
            # Get bounds if available, otherwise use defaults
            bounds = default_bounds.get(param_name, {'min': 0.0, 'max': 1.0})
            
            parameters[param_name] = {
                'value': param_value,
                'vary': vary,
                'min': bounds['min'],
                'max': bounds['max']
            }
    
    # Create full configuration structure
    config = {
        'data_file': 'rename_me.txt',  # Default, user should update
        'amp_type': 'AD745',
        'parameters': parameters
    }
    
    # Write the JSON configuration
    with open(json_output_file, 'w') as f:
        json.dump(config, f, indent=4)
    
    print(f"Successfully converted {txt_file} to {json_output_file}")
    print(f"Converted {len(parameters)} parameters")
    
    # Print summary
    varying_params = [name for name, param in parameters.items() if param['vary']]
    print(f"Parameters to fit: {len(varying_params)}")
    if varying_params:
        print(f"  Varying parameters: {', '.join(varying_params)}")
    
    return config

def main():
    parser = argparse.ArgumentParser(description='Convert startparams.txt to JSON configuration format')
    parser.add_argument('input_file', type=str, help='Input startparams.txt file')
    parser.add_argument('output_file', type=str, nargs='?', 
                       help='Output JSON file (default: config.json)')
    
    args = parser.parse_args()
    
    if not Path(args.input_file).exists():
        print(f"Error: Input file '{args.input_file}' not found")
        return 1
    
    # Default output file if not specified
    if args.output_file is None:
        args.output_file = 'config.json'
    
    try:
        convert_txt_to_json(args.input_file, args.output_file)
        print(f"\nIMPORTANT: Please update the 'data_file' and 'amp_type' fields in {args.output_file} to match your experimental data file.")
        return 0
    except Exception as e:
        print(f"Error during conversion: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())