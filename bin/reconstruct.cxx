//
// File: reconstruct.cxx
//
// Author: Hunter Sullivan
//
// Decription: Reconstruction software for sipmwheel.
//

// Includes
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <experimental/filesystem>
#include "TFile.h"
#include "Utilities.h"
#include "FileReader.h"
#include "Reconstructor.h"
#include "ConfigReader.h"

// Prototypes
void PrintTheFiles(const wheel::SiPMToFilesMap& sipmToFilesMap);
void GetTheRecoFiles(wheel::SiPMToFilesMap& map, const wheel::Configuration& config);
void Reconstruct(const wheel::Configuration& myConfig);
void FillSiPMInfo(wheel::SiPMInfoMap& sipmInfoMap, const wheel::Configuration& config);
void InitializeOutputFiles(const wheel::Configuration& config);

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
  if (myConfig.process != "reco") 
  {
    std::cout << "Error. Please list 'reco' under 'process' in configuration file." << std::endl;
    exit(1); 
  }
  configReader.OutputConfigInfo(myConfig);

  // Initialize output files
  InitializeOutputFiles(myConfig);

  // Start reconstuction 
  Reconstruct(myConfig);
  return 0;
}

void InitializeOutputFiles(const wheel::Configuration& config)
{
  // Reconstruction output file
  TFile f(config.recoOutputPath.c_str(), "RECREATE");
  f.Close();
  
  // Waveform files
  TFile f2(config.rawWaveformPath.c_str(), "RECREATE");
  f2.Close();
  TFile f3(config.modWaveformPath.c_str(), "RECREATE");
  f3.Close();
}

void Reconstruct(const wheel::Configuration& myConfig)
{
  // Fill the sipm info first
  wheel::SiPMInfoMap sipmInfoMap;
  FillSiPMInfo(sipmInfoMap, myConfig);

  std::cout << "SiPM Information: \n";
  for (const auto& sipm : sipmInfoMap)
  {
    std::cout << "SiPM " << sipm.first << "  Bias = " << sipm.second.bias << "  Gain = " << sipm.second.gain << "  Breakdown = " << sipm.second.breakdown << std::endl;
  }

  // First we need to get the files to read data from
  std::cout << "\nGetting the files from: " << myConfig.pathToData << std::endl;
  wheel::SiPMToFilesMap sipmToFilesMap;
  GetTheRecoFiles(sipmToFilesMap, myConfig);

  // Make sure there is the same amount of data for each sipm
  const unsigned& nFiles = sipmToFilesMap.find(1)->second.size();
  for (const auto& sipm : sipmToFilesMap)
  {
    if (sipm.second.size() != nFiles) { std::cout << "Error. Different data sizes for SiPMs\n" << std::endl; exit(1); }
  }

  // Option to print the collected files
  if (myConfig.printFiles) PrintTheFiles(sipmToFilesMap);

  // Main workhourse here
  // Start reco
  for (unsigned trigger = 1; trigger <= nFiles; trigger++)
  {
    std::cout << "\n" << std::setfill('*') << std::setw(80) << "" << std::setfill(' ')  << std::endl; 
    std::cout << "\nInitializing event " << trigger << "..." << std::endl;
    // Now read the files
    std::cout << "Reading the files... " << std::endl;
    wheel::SiPMToTriggerMap sipmToTriggerMap;
    wheel::FileReader fr;
    fr.ReadFiles(sipmToTriggerMap, sipmToFilesMap, trigger, sipmInfoMap, myConfig);
    // Make sure we just have 1 file/sipm here
    for (const auto& sipm : sipmToTriggerMap)
    {
      if (sipm.second.size() != 1) { std::cout << "Error. There is not 1 file specified for each sipm\n." << std::endl; exit(1); }
    } 
 
    // Start reconstruction
    wheel::Reconstructor Reconstructor;
    Reconstructor.Reconstruct(sipmToTriggerMap, sipmInfoMap, myConfig, trigger);
  }
}

void FillSiPMInfo(wheel::SiPMInfoMap& sipmInfoMap, const wheel::Configuration& config)
{
  // Loop over sipms
  for (unsigned sipm = 1; sipm <= config.nSiPMs; sipm++)
  {
    wheel::SiPMInfo sipmInfo;
    sipmInfo.gain      = config.gains.find(sipm)->second;
    sipmInfo.breakdown = config.breakdowns.find(sipm)->second; 
    sipmInfo.bias      = *config.biases.begin();

    sipmInfoMap.emplace(sipm, sipmInfo);
  }
}

void GetTheRecoFiles(wheel::SiPMToFilesMap& map, const wheel::Configuration& config)
{
  namespace stdfs = std::experimental::filesystem;

  for (unsigned sipm = 1; sipm <= config.nSiPMs; sipm++)
  {
    // What is the directory
    std::string sipmDir = config.pathToData + "/sipm" + std::to_string(sipm);
    std::set<std::string> sipmData;
    // Create a temp map for storing these
    std::set<std::string> tempMap; 
     
    // Only get the bias directory listed in config
    for(const auto& biasDirIter : stdfs::directory_iterator(sipmDir.c_str()))
    {
      // First get the bias
      std::string bias;
      std::string path = biasDirIter.path();
      for (std::string::iterator ch = path.end(); *ch != '/'; ch--) bias = *ch + bias;
      const float thisBias = std::stof(bias);

      if (thisBias != *config.biases.begin()) continue;
      
      // Loop over files?
      //std::set<std::string> files;
      for (const auto& fileIter : stdfs::directory_iterator(biasDirIter.path().c_str()))
      {
        // Store these into our set
        tempMap.insert(fileIter.path());
      }
      //tempMap.insert(files);
    }
    // Insert into our map
    map.emplace(sipm, tempMap);
  }
}

void PrintTheFiles(const wheel::SiPMToFilesMap& map)
{
  std::cout << "\nReading the following files: " << std::endl;
  for (const auto& sipmData : map)
  {
    std::cout << sipmData.second.size() << " files for SiPM " << sipmData.first << ":" << std::endl;
    for (const auto& file : sipmData.second)
    {
      std::cout << file << std::endl;
    }
    std::cout << std::endl;
  }
}
