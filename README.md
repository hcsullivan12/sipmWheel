# sipmwheel
This repository holds the reconstruction code for the SiPM wheel. The goal is to reconstruct the incident position of light on an acrylic disk which has SiPMs positioned around the periphery.

# About the code
## config/WheelConfiguration.txt
This file sets the configuration for the algorithms (characterization or reconstruction). Among the settings are number of SiPMs, disk radius, exact algorithm to run, sipm biases, gains, output paths, and others. This text file interface was used for user-friendliness. Be sure to keep the format. The data must be stored in specific format based on the process being executed.

After storing the data files to read from, the waveform is smoothed out using a moving average. Subsequently, hit candidates are found on each waveform. There are various parameters to tune hit finding and smoothing in the config file.

## Characterization
Characterization takes all the hits and creates an amplitude distribution for each bias voltage and sipm. Peaks are found in the distribution and the gain is fit. Results are output to the characteriationOutputFile in config. Create a data directory with subdirectories for each sipm and subsubdirectories for each bias voltage.
Example:

```
ls data/
sipm1  sipm2  sipm3  sipm4
ls data/sipm1
73.5  74.0
```

The algorithm will then read from each file in these directories. Make sure to list characterization under "process" in config. 

## Reconstruction
Reconstruction takes all the hits from each sipm for a specific event and performs a maximum likelihood estimate for the light position, disk attenuation length, and original light yield. Results are output to the recoOutputFile in config. You must also list "reconstruction" as the process and the sipm gains, bias, and breakdowns in config.  

# Running the code
Once the configuration has been set, run the following commands:

```
mkdir build
cd build
cmake ..
make
./sipmwheel <path_to_config_file>
```

Note: The characteriation and reconstruction plots are updated throughout the algorithm. So if you're starting a new run, make sure to delete it!
