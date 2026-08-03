# Quick Start Guide - Python CEB IV Fitting

This guide will get you up and running with the Python implementation of the CEB IV curve fitting software.

## Step 1: Installation

### Prerequisites
- Python 3.7 or higher
- pip package manager

### Install Dependencies
```bash
# Navigate to the project directory
cd C:\work\sci\ceb_iv_fit

# Install required packages
pip install -r requirements.txt
```

### Optional: Install as Package
```bash
pip install -e .
```
This installs the package with command-line shortcuts:
- `ceb-fit` instead of `python main.py`
- `ceb-verify` for verification
- `ceb-convert` for configuration conversion

## Step 2: Prepare Configuration

### Option A: Use Default Configuration
```bash
# Default config.json is provided
python main.py --config config.json
```

### Option B: Convert Existing C++ Configuration
If you have `startparams.txt` from the C++ version:
```bash
python convert_config.py startparams.txt my_config.json
```

### Option C: Create Custom Configuration
Edit `config.json` to match your needs:

```json
{
    "data_file": "path/to/your/experimental_data.txt",
    "amp_type": "AD745",
    "parameters": {
        "Pbg": {"value": 0.0, "vary": false, "min": 0.0, "max": 1.0},
        "beta": {"value": 0.111, "vary": true, "min": 0.0, "max": 1.0},
        "Tp": {"value": 0.19, "vary": true, "min": 0.1, "max": 0.4}
        // ... other parameters
    }
}
```

## Step 3: Prepare Your Data

### Data Format
Your experimental data should be a two-column text file:
```
# Voltage (V)    Current (A)
0.000000e+00    0.000000e+00
2.000000e-06    1.500000e-07
4.000000e-06    3.100000e-07
// ... more rows
```

### Column Format
- **First column**: Voltage in Volts
- **Second column**: Current in Amperes
- **Delimiter**: Tab-separated (default) or use custom delimiter

### File Locations
Place your data file in the root directory or update the path in `config.json`.

## Step 4: Run the Fitting

### Basic Run
```bash
python main.py --config config.json
```

### Advanced Options

#### Use Different Fitting Method
```bash
# Use LMFIT (recommended for better convergence)
python main.py --config config.json --method lmfit --runs 5
```

#### Specify Custom Data File
```bash
python main.py --config config.json --data my_experiment.txt
```

#### Remove Voltage Offset
```bash
python main.py --config config.json --remove-offset
```

#### Save Final Configuration
```bash
python main.py --config config.json --save-config fitted_results.json
```

## Step 5: Understand the Output

### Console Output
```
Initializing IV Parameter Fitter...
Parameters loaded from config:
  beta = 0.111000, to fit = True
  Z = 0.500000, to fit = True
  Tp = 0.190000, to fit = True
  ...

Loading experimental data from 'data.txt'...
Loaded 1000 data points

Computing CEB properties...
  1/499: V:1.234567e-04    I:1.234567e-09    Sv:1.234567e+06    Te:1.234567e-01
  2/499: V:2.345678e-04    I:2.345678e-09    Sv:2.345678e+06    Te:2.345678e-01
  ...

Starting fitting with golden method (3 runs)...
SeqFit run 0
  beta: 1.110000e-01 -> 1.123456e-01, fmin = 1.234567e-04
  Z: 5.000000e-01 -> 4.876543e-01, fmin = 1.123456e-04
  Tp: 1.900000e-01 -> 1.923456e-01, fmin = 1.012345e-04

Final parameters saved to 'config_fitted.json'

Final Parameters:
  beta: 1.123456e-01 (fit)
  Z: 4.876543e-01 (fit)
  Tp: 1.923456e-01 (fit)
  ...

Finished. Total time: 72.45 seconds
```

### Generated Files

1. **`Te.txt`**: Electron temperature vs voltage
   - Columns: Voltage, Current, Iqp, Iand, V/Rleak, Te, Ts, DeltaT, Peph, Pand, Pleak, Pabs, Pcool

2. **`Noise.txt`**: Noise characteristics
   - Columns: Voltage, NOISEep, NOISEs, NOISEa, NOISE, NOISEph, NOISE²-NOISEph²

3. **`NEP.txt`**: Noise equivalent power
   - Columns: Voltage, Current, NEPeph, NEPs, NEPa, NEP, NEPph, Sv, NEP²-NEPph²

4. **`G.txt`**: Thermal conductance
   - Columns: Voltage, Ge, Gnis

5. **`fitparameters_new.txt`**: Fitting results and parameters
   - Appended with timestamp for each run

6. **`config_fitted.json`**: Updated parameter values
   - For use in future runs

## Step 6: Troubleshooting

### Issue: Import Errors
```bash
# Error: Module not found
# Solution: Ensure you're in the project directory
cd C:\work\sci\ceb_iv_fit
python main.py
```

### Issue: Configuration Not Found
```bash
# Error: config.json not found
# Solution: Use default or create custom config
python main.py  # Will create default config.json
```

### Issue: Data Format Errors
```bash
# Error: Unable to read experimental data
# Solution: Check file format - should be tab-separated text with two columns
# Use this to convert:
awk '{print $1 "\t" $2}' data.txt > data_fixed.txt
```

### Issue: Performance Too Slow
```bash
# Solution: Reduce computational parameters in config.json
"dV": 4e-6,  # Increase voltage step
"dVFinVg": 0.9,  # Reduce voltage range
```

### Issue: Fitting Not Converging
```bash
# Solution: Try different methods
python main.py --config config.json --method lmfit --runs 10

# Or adjust parameter bounds in config.json
"beta": {"min": 0.05, "max": 0.2}  # Narrow the search range
```

## Step 7: Verification

### Run Tests
```bash
# Basic implementation tests
python test_implementation.py

# Verification against C++ version
python verify_implementation.py --full-verification
```

### Compare Results
```bash
# If you have C++ output
python verify_implementation.py --compare Te.txt Te_cpp.txt
```

## Step 8: Integration Examples

### Use as Python Module
```python
from python.iv_param_fitter import IVParamFitter

# Initialize and run
fitter = IVParamFitter(config_file='config.json')
fitter.load_experiment_data('data.txt', remove_offset=False)
fitter.compute_ceb_properties()

# Single parameter fit
fitter(['beta'], ['Tp'])  # Fit only these parameters

# Save results
fitter.save_config('results.json')
```

### Batch Processing
```bash
# Process multiple datasets
for data_file in data/*.txt; do
    python main.py --config config.json --data "$data_file" \
                --save-config "results_$(basename $data_file .txt).json"
done
```

## Common Workflows

### Research Workflow
```bash
# 1. Generate initial fit
python main.py --config config.json --method golden --runs 3

# 2. Refine with better method
python main.py --config config_fitted.json --method lmfit --runs 5

# 3. Verify quality
python verify_implementation.py --full-verification

# 4. Analyze results
import matplotlib.pyplot as plt
te_data = np.loadtxt('Te.txt', delimiter='\t')
plt.plot(te_data[:, 0], te_data[:, 1])  # Voltage vs Current
plt.savefig('iv_curve.png')
```

### Production Workflow
```bash
# Automated processing with error handling
#!/bin/bash
for file in experiments/week_*.txt; do
    if python main.py --config config.json --data "$file" \
                    --save-config "final_$(basename $file .txt).json"; then
        echo "SUCCESS: $file processed"
    else
        echo "FAILED: $file processing"
    fi
done
```

## Getting Help

### Command Line Help
```bash
python main.py --help
```

### Documentation
- See `README_PYTHON.md` for comprehensive documentation
- See `MIGRATION_GUIDE.md` for C++ to Python migration
- See code docstrings for API documentation

### Issues and Questions
1. Check generated log files (`fitparameters_new.txt`, `converg.txt`)
2. Run verification tests (`python test_implementation.py`)
3. Review parameter convergence in console output

## Next Steps

1. **Customize Parameters**: Edit `config.json` for your specific experiment
2. **Optimize Performance**: Adjust computational parameters for your hardware
3. **Automate Workflows**: Create scripts for batch processing
4. **Visualization**: Integrate with matplotlib for result visualization
5. **Analysis**: Use pandas/NumPy for post-processing results

Congratulations! You're now ready to use the Python CEB IV fitting software.