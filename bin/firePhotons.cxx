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
  ReadConfigFile(myConfig);
  OutputConfigInfo(myConfig);

  // Make sure we've configured simulation 
  if      (myConfig.process == "simulate") Simulate(myConfig);
  else    { std::cout << "Error. Must specify the 'simulate' in configuration file!\n" << std::endl; exit(1); }

  return 0;
}

void OutputConfigInfo(const wheel::Configuration& config)
{
  // Hello there!
  std::cout << std::setfill('-') << std::setw(80) << "-" << std::setfill(' ')  << std::endl;
  std::cout << "              SiPM Wheel Simulation Code                "      << std::endl;
  std::cout << "                 Version: " << sipmWheelVersion                << std::endl;
  std::cout << "         Author: Hunter Sullivan (UT Arlington)         "      << std::endl;
  std::cout                                                                    << std::endl;
  std::cout << "SiPM Wheel Configuration:\n";
 
  std::cout << "Process                " << config.process                        << std::endl
            << "SimulateOutputPath     " << config.simulateOutputPath             << std::endl
            << "NPhotonsToLaunch       " << config.nPhotonsToLaunch               << std::endl
            << "NumberOfSiPMs          " << config.nSiPMs                         << std::endl
            << "N0                     " << config.N0                             << std::endl
            << "Reconstruct            " << config.reconstruct                    << std::endl
            << "SiPMArea               " << config.sipmArea                       << std::endl
            << "DiskRadius             " << config.diskRadius                     << std::endl
            << "DiskThickness          " << config.diskThickness                  << std::endl
            << "SourcePosition         " << config.sourcePosition[0] << " " << config.sourcePosition[1] << std::endl
            << "SourceSigma            " << config.sourceSigma                    << std::endl
            << "SmearSigma             " << config.smearSigma                     << std::endl
            << "IndexRefractionDisk    " << config.indexRefractionDisk            << std::endl
            << "IndexRefractionEnv     " << config.indexRefractionEnv             << std::endl
            << "PermittivityDisk       " << config.permittivityDisk               << std::endl
            << "PermittivityEnv        " << config.permittivityEnv                << std::endl
            << "BulkAbsorption         " << config.bulkAbsorption                 << std::endl
            << "SurfaceAbsorptionCoeff " << config.surfaceAbsorptionCoeff         << std::endl
            << "TerminationThreshold   " << config.terminationThreshold           << std::endl
            << "MaxIterations          " << config.maxIterations                  << std::endl
            << "nVoxels                " << config.nVoxels                        << std::endl
            << "attenuationLength      " << config.attenuationLength              << std::endl
            << "recoOutputPath         " << config.recoOutputPath                 << std::endl;
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
    else if (header == "nSiPMs")                  config.nSiPMs = stoi(value);
    else if (header == "N0")                      config.N0     = stoi(value);
    else if (header == "reconstruct")             value == "true" ? config.reconstruct = true : config.reconstruct = false;
    else if (header == "sipmArea")                config.sipmArea = stof(value);
    else if (header == "sourcePosition")          RecordPosition(config, value);
    else if (header == "sourceSigma")             config.sourceSigma = std::stof(value);
    else if (header == "smearSigma")              config.smearSigma  = std::stof(value);
    else if (header == "diskRadius")              config.diskRadius = std::stof(value);
    else if (header == "diskThickness")           config.diskThickness = std::stof(value);
    else if (header == "tpbEmissionPeak")         config.tpbEmissionPeak = std::stof(value); 
    else if (header == "indexRefractionDisk")     config.indexRefractionDisk = std::stof(value);
    else if (header == "indexRefractionEnv")      config.indexRefractionEnv  = std::stof(value);
    else if (header == "permittivityDisk")        config.permittivityDisk    = std::stof(value);
    else if (header == "permittivityEnv")         config.permittivityEnv     = std::stof(value);
    else if (header == "bulkAttenuation")         config.bulkAttenuation     = std::stof(value);
    else if (header == "bulkAbsorption")          config.bulkAbsorption      = std::stof(value);
    else if (header == "surfaceAbsorptionCoeff")  config.surfaceAbsorptionCoeff = std::stof(value);
    else if (header == "nPhotonsToLaunch")        config.nPhotonsToLaunch      = std::stoi(value);
    else if (header == "terminationThreshold")    config.terminationThreshold  = std::stof(value);
    else if (header == "maxIterations")           config.maxIterations         = std::stoi(value);
    else if (header == "nVoxels")                 config.nVoxels               = std::stoi(value);
    else if (header == "attenuationLength")       config.attenuationLength     = std::stof(value);
    else if (header == "recoOutputPath")          config.recoOutputPath        = value;
    else    { std::cout << "Cannot identify " << header << std::endl; exit(1); }
  }
  
  // Place a series of safety nets here
  // TODO: May need to change these comparisons!
  //       Not good to compare floats!
  if (config.nSiPMs <= 1)              { std::cout << "Error. Please specify more than 1 sipm (even number)\n" << std::endl; exit(1); }
  if (config.process != "simulate")    { std::cout << "Error. Please specify 'simulate' as the process in config.\n" << std::endl; exit(1); }
  if (config.diskRadius <= 0)          { std::cout << "Error. Disk radius == 0!\n" << std::endl; exit(1); }
  if (config.diskThickness <= 0)       { std::cout << "Error. Disk thickness == 0!\n" << std::endl; exit(1); }
  if (config.tpbEmissionPeak <= 0)     { std::cout << "Error. TPB emission peak == 0!\n" << std::endl; exit(1); }
  if (config.indexRefractionDisk <= 0) { std::cout << "Error. Index of refraction disk == 0!\n" << std::endl; exit(1); }
  if (config.indexRefractionEnv <= 0)  { std::cout << "Error. Index of refraction medium == 0!\n" << std::endl; exit(1); }
  if (config.permittivityDisk <= 0)    { std::cout << "Error. Permittivity of disk <= 0!\n" << std::endl; exit(1); }
  if (config.permittivityEnv <= 0)     { std::cout << "Error. Permittivity of environment <= 0!\n" << std::endl; exit(1); }
  if (config.sourcePosition.size() < 2) { std::cout << "Error. Please list sourcePosition as r,theta in configuration\n" << std::endl; exit(1); }
  if (config.sourceSigma == 0)         { std::cout << "Error. Source sigma == 0!\n" << std::endl; exit(1); }
  if (config.smearSigma == 0)          { std::cout << "Error. Smear sigma == 0!\n" << std::endl; exit(1); }
  if (config.nPhotonsToLaunch <= 0) { std::cout << "Error. Please specify number of photons to simulate.\n" << std::endl; exit(1); }
  if (config.bulkAttenuation <= 0)  { std::cout << "Error. Please specify bulkAttenuation > 0\n" << std::endl; exit(1); }
  if (config.bulkAbsorption < 0) { std::cout << "Error. Please specify a bulkAbsorption > 0!\n" << std::endl; exit(1); }
  if (config.surfaceAbsorptionCoeff < 0 || config.surfaceAbsorptionCoeff > 1) { std::cout << "Error. Please specify a 0 < surfaceAbsorptionCoeff < 1.\n" << std::endl; exit(1); }
  if (config.terminationThreshold <= 0) { std::cout << "Error. Please specify termination threshold > 0.\n" << std::endl; exit(1); } 
}

void RecordPosition(wheel::Configuration& config, const std::string& value)
{ 
  std::stringstream linestream(value);
  std::string rTheta;
  while(std::getline(linestream, rTheta, ',')) config.sourcePosition.push_back(std::stof(rTheta));
}

void Simulate(const wheel::Configuration& config)
{
  // Start the simulation
  wheel::Simulator simulator(config);
  simulator.Simulate(config);
}

