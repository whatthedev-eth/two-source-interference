# Two-source Wave Interference Simulation

A.k.a. "Makin' Waves in Cairo!"

Simulated here is the wave interference pattern produced by two wave sources with identical wavelength `lambda`. 
- A "top-view" 2-D intensity plot is created to show the pattern inside a square plot of side length 1000 (arbitrary units). 
- The origin is at the center of the left edge of the plot. 
- The sources are placed a distance `d` apart, on the left edge of the plot symmetrically about the origin.


## User Inputs

Only these *integer* inputs for now:
- `num_pts` = number of points along each axis of the square grid, where 2 <= `num_pts` <= 25 (max value limited by allowed n_steps on testnet)
- `lambda` >= 1
- `d` >= 0

After input, both `lambda` and `d` are scaled up to be FP (fixed point) values where FP scale is `SCALE_FP` = 10**20.

## Cairo files

**intensity_plot.cairo** contains:
- `intensity_plot_arr`- The view function which accepts user inputs and initiates all calculations
- Other functions to calculate x coordinates, y coordinates, intensities at each position (x, y), and fill arrays for each of these

**wave_physics.cairo** contains functions for wave physics calculations of: 
- Wave function
- Sum of wave functions
- Intensity of combined wave

**math.cairo** contains math functions to use with fixed point quantities: 
- Square root, multiplication, division
- Distance between two points
- Shift any angle `theta` to equivalent angle `theta_shifted` within range -pi to pi
- Cosine approximation using Taylor series, requires -pi <= angle <= pi

**structs.cairo** contains struct definitions for: 
- Wave parameters common to both sources
- Wave parameters for an individual source
- Combined parameters for two waves

**constants.cairo** contains: 
- FP math constants
- Other math constants
- Physical wave parameters (except those that are user inputs)
- Plot parameters (except those that are user inputs)


## Tests

**test_intensity_plot.py**
- Contains input parameters which are fed into both Cairo and Python calculations
- Contains other constants and parameters (that should match those in **constants.cairo**) for Python calculations 
- Contains Python calculations of intensity array, to compare to Cairo calculations done by calling `intensity_plot_arr` function from **intensity_plot.cairo** using same input parameters
- Dumps `intensity_plot_arr` return to **test_intensity_plot.json** to be used by Jupyter notebook
- Prints for comparison the values (as FP) from the intensity arrays found with Cairo and Python
- To run: `pytest -s tests/test_intensity_plot.py`


**test_intensity_plot.ipynb**
- Jupyter notebook which plots two intensity arrays:
    - Array calculated in Cairo by **intensity_plot.cairo**
    - Array calculated in Python by **test_intensity_plot.py**
- Gets data from .json files created in **test_intensity_plot.py**
- To run:
    - Install [Jupyter Notebook](https://jupyter.org/install)
    - `jupyter notebook`, follow instructions to open notebook in browser
    - Open file **test_intensity_plot.ipynb**


**test_math.py**
- Compares Python calculations to those found by calling `theta_shifter` and `cosine_8th` functions from **math.cairo**
- Note:
    - `n` = 5 is used to include 5 terms (to 8th order) in Python cosine approximation, to match approximation in `cosine_8th` function in **math.cairo**
    - `num_tests` = an odd integer number of tests to run
    - Values of `theta` to test are incremented by pi/4 in range:
        *-(pi/4)(num_tests-1)/2 <= `theta` <= +(pi/4)(num_tests-1)/2*
- To run: `pytest -s tests/test_math.py`