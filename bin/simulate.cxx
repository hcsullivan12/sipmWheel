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
#include "Simulator.h"
#include "ConfigReader.h"

// Preprocessing variables
#ifdef VERSION
  #define sipmWheelVersion VERSION
#endif

void ReadConfigFile(wheel::Configuration& config);
void OutputConfigInfo(const wheel::Configuration& config);
void RecordPosition(wheel::Configuration& config, const std::string& value);
void Simulate(const wheel::Configuration& config);

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
  
  wheel::ConfigReader configReader;
  configReader.ReadFile(myConfig);
   if (myConfig.process != "simulate") 
  {
    std::cout << "Error. Please list 'simulate' under 'process' in configuration file." << std::endl;
    exit(1); 
  }
  configReader.OutputConfigInfo(myConfig);

  // Start simulation
  Simulate(myConfig);

  return 0;
}

void Simulate(const wheel::Configuration& config)
{
  // Start the simulation
  wheel::Simulator simulator(config);
  simulator.Simulate(config);
}

