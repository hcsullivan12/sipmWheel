# sipmwheel
This repository holds the reconstruction, characterization, and simulation software for the SiPM wheel. The objective is to reconstruct the incident position of light on an acrylic disk using data from SiPMs positioned around the periphery.

# About the code
## config/WheelConfiguration<Process>.txt
This file sets the configuration for the algorithms (characterization, reconstruction, or simulation). Among the settings are number of SiPMs, disk radius, sipm biases, gains, output paths, and disk properties. This text file interface was used for user-friendliness. Be sure to keep the format. The data must be stored in specific format based on the process being executed.

The waveform analysis is based on waveform text files written from a Teledyne Lecroy oscilloscope. After storing the data files, the waveform is smoothed out (configurable) using a moving average. Subsequently, hit candidates are found on each waveform. There are various parameters to tune hit finding and smoothing in the config file.

## Characterization
Characterization takes all the hits and creates amplitude distributions for each bias voltage and sipm. Peaks are found in the distribution and the gain is extrapolated from fitting. Results are output to the characteriationOutputFile in config. Create a data directory with subdirectories for each sipm and subsubdirectories for each bias voltage.
Example:

```
ls data/
sipm1  sipm2  sipm3  sipm4
ls data/sipm1
73.5  74.0
```

The algorithm will then read from each file in these directories. Make sure to list characterization under "process" in config. 

## Reconstruction
Reconstruction takes all the hits from each sipm for a specific event and performs a maximum likelihood estimate for the light position. First, a grid search algorithm is implented to find the approximate location, and a refinement attempt is made using a NR method. Results are output to the recoOutputFile in config. You must list "reconstruction" as the process and the sipm gains, bias, and breakdowns in config.  

## Monte Carlo
The sipmwheel Monte Carlo simulates photon propogation from a configurable position on the disk assuming a Gaussian light beam profile, uniform TPB emission, bulk and surface absorption, and surface roughness. 

# Running the code
Once the configuration has been set, run the following commands:

```
mkdir build
cd build
cmake ..
make
```
And one of following for the chosen process:
```
./reconstruct <path_to_config_file>
./characterize <path_to_config_file>
./simulate <path_to_config_file>
```
