#!/usr/bin/env python3
import sys
import os
import time
import json
import argparse

# Add parent directory to path to import from python subdirectory
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from python.iv_param_fitter import IVParamFitter
from python.utils import Utils

def main():
    parser = argparse.ArgumentParser(
        description='CEB IV Curve Fitting - Python Implementation',
        epilog='Example: python main.py --display --method lmfit --runs 2'
    )
    parser.add_argument('--config', type=str, default='config.json',
                       help='Path to JSON configuration file')
    parser.add_argument('--data', type=str, 
                       help='Path to experimental data file (overrides config)')
    parser.add_argument('--method', type=str, choices=['golden', 'lmfit'], default='golden',
                       help='Fitting method to use')
    parser.add_argument('--runs', type=int, default=3,
                       help='Number of fitting runs')
    parser.add_argument('--remove-offset', action='store_true',
                       help='Remove voltage offset from data')
    parser.add_argument('--save-config', type=str,
                       help='Save current configuration to specified file')
    parser.add_argument('--display', action='store_true',
                       help='Display real-time fitting progress with matplotlib (requires matplotlib and display environment)')
    
    args = parser.parse_args()
    
    start_time = time.time()
    
    try:
        # Check if config file exists
        if not os.path.exists(args.config):
            print(f"Config file '{args.config}' not found. Creating default config...")
            sample_config = Utils.create_sample_config()
            Utils.save_json_config(args.config, sample_config)
            print(f"Default config created at '{args.config}'. Please edit and run again.")
            return 0
        
        # Initialize fitter
        print("Initializing IV Parameter Fitter...")
        fitter = IVParamFitter(config_file=args.config, display=args.display)
        
        # Override data file if specified
        if args.data:
            fitter.data_file = args.data
        
        # Save config if requested
        if args.save_config:
            fitter.save_config(args.save_config)
            print(f"Configuration saved to '{args.save_config}'")
        
        # Load experimental data
        print(f"Loading experimental data from '{fitter.data_file}'...")
        data_points = fitter.load_experiment_data(fitter.data_file, args.remove_offset)
        print(f"Loaded {data_points} data points")
        
        # Compute CEB properties
        print("Computing CEB properties...")
        fitter.compute_ceb_properties()
        
        # Resample data
        Irex, Vrex = fitter.resample()
        print(f"Resampled {len(Irex)} points")
        
        # Perform fitting
        print(f"\nStarting fitting with {args.method} method ({args.runs} runs)...")
        if args.method == 'golden':
            fitter.sequential_fit(run_count=args.runs)
        else:  # lmfit
            fitter.lmfit_sequential_fit(run_count=args.runs)
        
        # Save final parameters
        final_config_file = args.config.replace('.json', '_fitted.json')
        fitter.save_config(final_config_file)
        print(f"\nFinal parameters saved to '{final_config_file}'")
        
        # Print final parameters
        print("\nFinal Parameters:")
        for param_name, value in sorted(fitter.par.items()):
            fit_status = "fit" if fitter.to_fit.get(param_name, False) else "fixed"
            print(f"  {param_name}: {value:.6e} ({fit_status})")
        
        elapsed_time = time.time() - start_time
        print(f"\nFinished. Total time: {elapsed_time:.2f} seconds")
        
        # Handle display cleanup
        if args.display:
            print("Display window will remain open until closed manually.")
            print("Press Ctrl+C or close the window to exit.")
            try:
                import matplotlib.pyplot as plt
                plt.show(block=True)
            except KeyboardInterrupt:
                print("\nInterrupted by user")
            finally:
                fitter._close_display()
        
        return 0
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        
        # Clean up display if still active
        try:
            if 'fitter' in locals() and 'args' in locals() and getattr(args, 'display', False):
                fitter._close_display()
        except:
            pass
        
        return 1

if __name__ == "__main__":
    sys.exit(main())