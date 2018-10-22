//
// File: firePhotons.cxx
//
// Author: Hunter Sullivan
//
// Decription: Simulation code for sipmwheel
// 

// Includes
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include "Utilities.h"

// Preprocessing variables
#ifdef VERSION
  #define sipmWheelVersion VERSION
#endif

void ReadConfigFile(wheel::Configuration& config);
void OutputConfigInfo(wheel::Configuration& config);
void RecordPosition(wheel::Configuration& config, const std::string& value);
void Simulate(wheel::Configuration& config);

int main(int argc, char **argv)
{
  // First check to make sure the path to data is given
  if (argc < 2)
  {
    std::cerr << "\nUsage: " << argv[0] << " pathToConfigFile\n" << std::endl;
    exit(1);
  }

  // First set the configuration
  std::cout << "\nReading configuration file... " << std::endl;
  wheel::Configuration myConfig;
  myConfig.pathToConfig = argv[1];
  ReadConfigFile(myConfig);
  OutputConfigInfo(myConfig);

  // Make sure we've configured simulation 
  if      (myConfig.process == "simulate") Simulate(myConfig);
  else    { std::cout << "Error. Must specify the 'simulate' in configuration file!\n" << std::endl; exit(1); }

  return 0;
}

void OutputConfigInfo(wheel::Configuration& config)
{
  // Hello there!
  std::cout << std::setfill('-') << std::setw(80) << "-" << std::setfill(' ')  << std::endl;
  std::cout << "               SiPM Wheel Simulation Code               "      << std::endl;
  std::cout << "                 Version: " << sipmWheelVersion                << std::endl;
  std::cout << "         Author: Hunter Sullivan (UT Arlington)         "      << std::endl;
  std::cout                                                                    << std::endl;
  std::cout << "SiPM Wheel Configuration:\n";
 
  std::cout << "Process               " << config.process                        << std::endl
            << "SimulateOutputPath    " << config.simulateOutputPath             << std::endl
            << "NumerOfSiPMs          " << config.nSiPMs                         << std::endl
            << "DiskRadius            " << config.diskRadius                     << std::endl
            << "DiskThickness         " << config.diskThickness                  << std::endl
            << "SourcePosition        " << config.sourcePosition[0] << " " << config.sourcePosition[1] << std::endl
            << "IndexRefractionDisk   " << config.indexRefractionDisk            << std::endl
            << "IndexRefractionMedium " << config.indexRefractionMedium          << std::endl;
  std::cout << std::setfill('-') << std::setw(80) << "-" << std::setfill(' ') << std::endl;
  std::cout << std::endl;
}

void ReadConfigFile(wheel::Configuration& config)
{
  // Open the file
  std::ifstream file(config.pathToConfig.c_str());
  if (!file.is_open()) {
    std::cout << "Cannot open file: " << config.pathToConfig << std::endl;
    exit(1);
  }
  
  std::string header, value;
  while(std::getline(file, header, '=')) 
  {
    std::getline(file, value);

    if      (header == "process")                 config.process = value; 
    else if (header == "simulateOutputPath")      config.simulateOutputPath = value;
    else if (header == "sourcePosition")          RecordPosition(config, value); 
    else if (header == "diskRadius")              config.diskRadius = std::stof(value);
    else if (header == "diskThickness")           config.diskThickness = std::stof(value);
    else if (header == "tpbEmissionPeak")         config.tpbEmissionPeak = std::stof(value); 
    else if (header == "indexRefractionDisk")     config.indexRefractionDisk = std::stof(value);
    else if (header == "indexRefractionMedium")   config.indexRefractionMedium = std::stof(value);
    else    { std::cout << "Cannot identify " << header << std::endl; exit(1); }
  }
  
  // Place a series of safety nets here;
  if (config.process != "simulate")    { std::cout << "Error. Please specify 'simulate' as the process in config.\n" << std::endl; exit(1); }
  if (config.diskRadius == 0)          { std::cout << "Error. Disk radius == 0!\n." << std::endl; exit(1); }
  if (config.diskThickness == 0)       { std::cout << "Error. Disk thickness == 0!\n." << std::endl; exit(1); }
  if (config.tpbEmissionPeak == 0)     { std::cout << "Error. TPB emission peak == 0!\n." << std::endl; exit(1); }
  if (config.indexRefractionDisk == 0) { std::cout << "Error. Index of refraction disk == 0!\n." << std::endl; exit(1); }
  if (config.indexRefractionMedium == 0) { std::cout << "Error. Index of refraction medium == 0!\n." << std::endl; exit(1); }
  if (config.sourcePosition.size() != 0) { std::cout << "Error. Please list sourcePosition as r,theta in configuration\n" << std::endl; exit(1); }
}

void RecordPosition(wheel::Configuration& config, const std::string& value)
{ 
  std::stringstream linestream(value);
  std::string rTheta;
  while(std::getline(linestream, rTheta, ',')) config.sourcePosition.push_back(std::stof(rTheta));
}

void Simulate(wheel::Configuration& config)
{

}

