//
// File: characterize.cxx
//
// Author: Hunter Sullivan
//
// Decription: Characterization software for sipmwheel.
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
#include "Characterizer.h"
#include "ConfigReader.h"

// Prototypes
void PrintTheFiles(const wheel::SiPMToBiasTriggerMap& sipmToBiasTriggerMap);
void GetTheCharacterizationFiles(wheel::SiPMToBiasTriggerMap& map, const wheel::Configuration& config);
void Characterize(const wheel::Configuration& myConfig);
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
  if (myConfig.process != "characterize")
  { 
    std::cout << "Error. Please list 'characterize' under 'process' in configuration file." << std::endl; 
    exit(1); 
  }
  configReader.OutputConfigInfo(myConfig);

  // Initialize output files
  InitializeOutputFiles(myConfig);
 
  // Start characterization
  Characterize(myConfig);
  return 0;
}

void InitializeOutputFiles(const wheel::Configuration& config)
{
  if (config.process == "characterize") 
  {
    TFile f1(config.characterizeOutputPath.c_str(), "RECREATE");
    f1.Close();
  }
  if (config.process == "reco")
  {
    TFile f(config.recoOutputPath.c_str(), "RECREATE");
    f.Close();
  }
  
  // Waveform files
  TFile f2(config.rawWaveformPath.c_str(), "RECREATE");
  f2.Close();
  TFile f3(config.modWaveformPath.c_str(), "RECREATE");
  f3.Close();
}

void Characterize(const wheel::Configuration& myConfig)
{
  // Start a timer
  clock_t start = clock();

  // First we need to get the files to read data from
  std::cout << "Getting the files from: " << myConfig.pathToData << std::endl;
  wheel::SiPMToBiasTriggerMap sipmToBiasTriggerMap;
  GetTheCharacterizationFiles(sipmToBiasTriggerMap, myConfig);
  // Option to print the collected files
  if (myConfig.printFiles) PrintTheFiles(sipmToBiasTriggerMap);
  // Now read the files
  std::cout << "Reading the files... " << std::endl;
  // Loop over the sipms
  wheel::Characterizer ch;
  ch.Initialize(myConfig);
  for (auto& sipm : sipmToBiasTriggerMap)
  {
    wheel::SiPMToTriggerMap sipmToTriggerMap;
    wheel::FileReader fr;
    fr.ReadFiles(sipmToTriggerMap, sipm.second, sipm.first, myConfig);

    // Characterize
    wheel::SiPMInfoMap sipmInfoMap;
    ch.Characterize(sipmInfoMap, sipmToTriggerMap, myConfig);
  }
  // Save the plots
  ch.SaveCharacterizationPlots(myConfig);

  clock_t end = clock();
  double duration = ((double) (end - start)) / CLOCKS_PER_SEC;
  std::cout << "Run time of " << duration << " s" << std::endl;
}

void GetTheCharacterizationFiles(wheel::SiPMToBiasTriggerMap& map, const wheel::Configuration& config)
{
  namespace stdfs = std::experimental::filesystem;

  for (unsigned sipm = 1; sipm <= config.nSiPMs; sipm++)
  {
    // What is the directory
    std::string sipmDir = config.pathToData + "/sipm" + std::to_string(sipm);
    std::set<std::string> sipmData;
    // Create a temp map for storing these
    std::map<float, std::set<std::string>> tempMap; 
    
    // Loop over the bias directories
    for(const auto& biasDirIter : stdfs::directory_iterator(sipmDir.c_str()))
    {
      // First get the bias
      std::string bias;
      std::string path = biasDirIter.path();
      for (std::string::iterator ch = path.end(); *ch != '/'; ch--) bias = *ch + bias;
      const float thisBias = std::stof(bias);

      // Loop over files?
      std::set<std::string> files;
      for (const auto& fileIter : stdfs::directory_iterator(biasDirIter.path().c_str()))
      {
        // Store these into our set
        files.insert(fileIter.path());
      }
      tempMap.emplace(thisBias, files);
    }

    // Insert into our map
    map.emplace(sipm, tempMap);
  }
}

void PrintTheFiles(const wheel::SiPMToBiasTriggerMap& map)
{
  std::cout << "\nReading the following files: " << std::endl;
  for (const auto& sipmData : map)
  {
    std::cout << sipmData.second.size() << " biases for SiPM " << sipmData.first << ":" << std::endl;
    for (const auto& bias : sipmData.second)
    {
      std::cout << "Bias = " << bias.first << std::endl;
      for (const auto& file : bias.second) std::cout << file << std::endl;
    }
    std::cout << std::endl;
  }
}
