# IV Curve Fitting - Real-time Display Usage Guide

This document provides usage examples for the new `--display` feature that provides real-time visualization during fitting.

## Requirements

- **matplotlib**: Required for display functionality
- **Display environment**: X11 server on Linux/macOS, or native display on Windows
- If these requirements are not met, the display feature is automatically disabled

## Command-line Usage

### Basic Display Example
```bash
python main.py --display
```

### Custom Display Example
```bash
python main.py --display --method lmfit --runs 2 --config my_config.json
```

### Display with Different Methods
```bash
# Using golden section method
python main.py --display --method golden --runs 3

# Using lmfit method  
python main.py --display --method lmfit --runs 2
```

## Display Features

When you run the fitting with `--display`, a matplotlib window will open showing:

### Visual Elements
- **Dual display window**: Shows two side-by-side plots with identical data
  - **Left plot**: Log Y-scale IV curve
  - **Right plot**: Linear Y-scale IV curve
- **Blue circles with lines**: Experimental data points (from Irex, Vrex)
- **Red solid line**: Numerical fit curve (from Inum, Vnum)  
- **Grid overlay**: For easier reading of values
- **Axis labels**: Voltage (V) and Current (A)

### Real-time Information
- **Evaluation counter**: Shows total number of objective function evaluations
- **Chi-squared value**: Current χ² metric showing fit quality
- **Dynamic updates**: Display refreshes after each objective function evaluation

### Display Behavior
- The window opens with proper axis limits automatically
- Red curve updates in real-time as parameters are optimized
- Blue data points remain static (reference experimental data)
- Window remains open after fitting completes
- Use Ctrl+C or close the window to exit

## Programmatic Usage

You can also use the display feature programmatically:

```python
from python.iv_param_fitter import IVParamFitter

# Enable display during fitting
fitter = IVParamFitter(config_file='config.json', display=True)

# Load data and perform fitting
fitter.load_experiment_data('experiment_data.txt')
fitter.compute_ceb_properties()

# Both fitting methods support display
fitter.sequential_fit(run_count=3)
# or
fitter.lmfit_sequential_fit(run_count=2)

# Display will close automatically or use:
# fitter._close_display()
```

## Error Handling

The display feature includes robust error handling:

1. **Missing matplotlib**: Automatically disables display with informative message
2. **No display environment**: Checks for DISPLAY environment variable (Unix) or native display (Windows)
3. **Update failures**: Non-critical, fitting continues even if display updates fail
4. **Cleanup**: Properly closes display resources even if errors occur

## Performance Considerations

- Display updates happen after each objective function evaluation
- Updates are minimal (data plotting) to avoid significant performance impact
- Large computation tasks still benefit from C++ backend when enabled
- Display can be disabled completely if needed for maximum performance

## Troubleshooting

### Display doesn't appear
- Check that matplotlib is installed: `pip install matplotlib`
- Verify display environment is available
- Check console messages for error information

### Display is slow
- Display updates are optimized for minimal overhead
- Consider disabling display for very large datasets
- The display feature is optional and non-critical for fitting

### Window closes immediately
- The display window should remain open until you close it manually
- If using automated scripts, consider using `--display` for interactive sessions only

## Examples

### Quick test of display functionality
```bash
# Run a single fitting run to test display
python main.py --display --runs 1 --method golden
```

### Full fitting with display monitoring
```bash
# Comprehensive fitting with real-time visualization
python main.py --display --method lmfit --runs 3 --data experiment_data.txt
```

### Comparison of fitting methods with display
```bash
# Test golden section method
python main.py --display --method golden --runs 2

# Test lmfit method  
python main.py --display --method lmfit --runs 2

# Compare results and fitting behavior
```

## Benefits

The display feature provides several key advantages:

1. **Real-time feedback**: Watch the fit converge in real-time
2. **Parameter impact**: See how parameter changes affect the fit
3. **Quality assessment**: Monitor χ² convergence during fitting
4. **Debugging aid**: Identify potential issues early
5. **Educational value**: Understand the fitting process visually
6. **Performance insight**: See optimization progress and evaluation count