// 
// File: ConfigReader.cxx
//
// Author: Hunter Sullivan
//
// Discription: Structure to read the configuration file.
//

#include "ConfigReader.h"
#include <sstream>

// Preprocessing variables
#ifdef VERSION
#define sipmWheelVersion VERSION
#endif

namespace wheel {

ConfigReader::ConfigReader()
{}

ConfigReader::~ConfigReader()
{}

void ConfigReader::ReadFile(Configuration& config)
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

    if      (header == "pathToData")              config.pathToData      = value;
    else if (header == "rawWaveformPath")         config.rawWaveformPath = value;
    else if (header == "modWaveformPath")         config.modWaveformPath = value;
    else if (header == "recoOutputPath")          config.recoOutputPath  = value;
    else if (header == "printFiles")              value == "true" ? config.printFiles       = true : config.printFiles       = false;
    else if (header == "baselineSubtract")        value == "true" ? config.baselineSubtract = true : config.baselineSubtract = false;
    else if (header == "saveRawWaveforms")        value == "true" ? config.saveRawWaveforms = true : config.saveRawWaveforms = false;
    else if (header == "saveModWaveforms")        value == "true" ? config.saveModWaveforms = true : config.saveModWaveforms = false;
    else if (header == "nSiPMs")                  config.nSiPMs          = std::stoi(value);
    else if (header == "smoothWaveform")          value == "true" ? config.smoothWaveform = true : config.smoothWaveform = false;
    else if (header == "smaRange")                config.smaRange        = std::stoi(value);
    else if (header == "hitFinderSearch")         config.hitFinderSearch = std::stoi(value);
    else if (header == "resolution")              config.resolution      = std::stoi(value);
    else if (header == "minimumHitAmp")           config.minimumHitAmp   = std::stof(value);
    else if (header == "process")                 value == "characterize" ? config.process  = "characterize" : config.process = "reco"; // default to reco
    else if (header == "characterizeAmpThr")      config.characterizeAmpThr      = std::stof(value);
    else if (header == "characterizeAmpSig")      config.characterizeAmpSig      = std::stof(value);
    else if (header == "characterizeAmpFitRange") config.characterizeAmpFitRange = std::stof(value);
    else if (header == "characterizeOutputPath")  config.characterizeOutputPath  = value;
    else if (header == "nFilesCharacterize")      config.nFilesCharacterize      = std::stoi(value);
    else if (header == "nBiases")                 config.nBiases                 = std::stoi(value);
    else if (header == "biases")                  RecordBiases(config, value); 
    else if (header == "gains")                   RecordGains(config.gains, config, value);
    else if (header == "breakdowns")              RecordGains(config.breakdowns, config, value);
    else if (header == "diskRadius")              config.diskRadius = std::stof(value);
    else if (header == "attenuationLength")       config.attenuationLength = std::stof(value); 
    else if (header == "nVoxels")                 config.nVoxels    = std::stoi(value);
    else if (header == "maxIterations")           config.maxIterations = std::stoi(value);
    else    { std::cout << "Cannot identify " << header << " listed in configuration file." << std::endl; exit(1); }
  }
  
  // Place a series of safety nets here;
  if (config.nBiases != config.biases.size())                                { std::cout << "Error. Number of biases does not match in config.\n" << std::endl; exit(1); }
  if (config.process == "reco" && config.nSiPMs != config.gains.size())      { std::cout << "Error. Number of gains does not match number of SiPMs\n." << std::endl; exit(1); }
  if (config.process == "reco" && config.nSiPMs != config.breakdowns.size()) { std::cout << "Error. Number of breakdowns does not match number of SiPMs\n." << std::endl; exit(1); }
  if (config.process == "reco" && config.biases.size() != 1)                 { std::cout << "Error. Can only specify one bias for reconstruction\n." << std::endl; exit(1); }
}

void ConfigReader::OutputConfigInfo(const Configuration& config)
{
  // Hello there!
  std::cout << std::setfill('-') << std::setw(80) << "-" << std::setfill(' ')  << std::endl;
  std::cout << "     SiPM Wheel Characterization and Analysis Code      "      << std::endl;
  std::cout << "                 Version: " << sipmWheelVersion                << std::endl;
  std::cout << "         Author: Hunter Sullivan (UT Arlington)         "      << std::endl;
  std::cout                                                                    << std::endl;
  std::cout << "SiPM Wheel Configuration:\n";
  if (config.process == "reco") {
  std::cout << "Process             " << config.process                        << std::endl
            << "PathToData          " << config.pathToData                     << std::endl
            << "SaveRawWaveforms    " << config.saveRawWaveforms               << std::endl
            << "SaveModWaveforms    " << config.saveModWaveforms               << std::endl
            << "RawWaveformsPath    " << config.rawWaveformPath                << std::endl
            << "ModWaveformsPath    " << config.modWaveformPath                << std::endl
            << "RecoOutputPath      " << config.recoOutputPath                 << std::endl
            << "BaselineSubtract    " << config.baselineSubtract               << std::endl
            << "SmoothWaveform      " << config.smoothWaveform                 << std::endl
            << "SMARange            " << config.smaRange                       << std::endl
            << "WaveformResolution  " << config.resolution                     << std::endl
            << "HitFinderSearch     " << config.hitFinderSearch                << std::endl
            << "MinimumHitAmp       " << config.minimumHitAmp                  << std::endl
            << "NumerOfSiPMs        " << config.nSiPMs                         << std::endl
            << "Biases              ";for(const auto& bias : config.biases)     std::cout << bias        << "  ";
                                                                                std::cout << std::endl; std::cout 
            << "Gains               ";for(const auto& sipm : config.gains)      std::cout << sipm.second << "  ";
                                                                                std::cout << std::endl; std::cout  
            << "Breakdowns          ";for(const auto& sipm : config.breakdowns) std::cout << sipm.second << "  ";
                                                                                std::cout << std::endl; std::cout 
            << "nVoxels             " << config.nVoxels                        << std::endl
            << "maxIterations       " << config.maxIterations                  << std::endl
            << "AttenuationLength   " << config.attenuationLength              << std::endl
            << "DiskRadius          " << config.diskRadius                     << std::endl;
  std::cout << std::setfill('-') << std::setw(80) << "-" << std::setfill(' ') << std::endl;
  std::cout << std::endl;
  return;
  }
  std::cout << "Process             " << config.process                        << std::endl
            << "PathToData          " << config.pathToData                     << std::endl
            << "OutputPath          " << config.characterizeOutputPath         << std::endl
            << "NFiles/SiPM         " << config.nFilesCharacterize             << std::endl
            << "SaveRawWaveforms    " << config.saveRawWaveforms               << std::endl
            << "SaveModWaveforms    " << config.saveModWaveforms               << std::endl
            << "RawWaveformsPath    " << config.rawWaveformPath                << std::endl
            << "ModWaveformsPath    " << config.modWaveformPath                << std::endl
            << "BaselineSubtract    " << config.baselineSubtract               << std::endl
            << "SmoothWaveform      " << config.smoothWaveform                 << std::endl
            << "SMARange            " << config.smaRange                       << std::endl
            << "WaveformResolution  " << config.resolution                     << std::endl
            << "HitFinderSearch     " << config.hitFinderSearch                << std::endl
            << "MinimumHitAmp       " << config.minimumHitAmp                  << std::endl
            << "NumerOfSiPMs        " << config.nSiPMs                         << std::endl
            << "AmpThreshold        " << config.characterizeAmpThr             << std::endl
            << "AmpSigma            " << config.characterizeAmpSig             << std::endl
            << "AmpFitRange         " << config.characterizeAmpFitRange        << std::endl;
  std::cout << std::setfill('-') << std::setw(80) << "-" << std::setfill(' ') << std::endl;	
  std::cout << std::endl;
}

void ConfigReader::RecordGains(std::map<unsigned, float>& map, const Configuration& config, const std::string& value)
{
  // Make sure n sipms has been defined
  if (!config.nSiPMs) { std::cout << "Error. Please specify number of sipms before listing gains and breakdowns.\n" << std::endl; exit(1); }
  
  std::stringstream linestream(value);
  std::string gain;
  unsigned sipm = 1;
  while(std::getline(linestream, gain, ',')) { map.emplace(sipm, std::stof(gain)); sipm++; }
}


void ConfigReader::RecordBiases(Configuration& config, const std::string& value)
{
  // Make sure n biases has been defined
  if (!config.nBiases) { std::cout << "Error. Please specify number of biases before listing biases.\n" << std::endl; exit(1); }
  
  std::stringstream linestream(value);
  std::string bias;
  while(std::getline(linestream, bias, ',')) config.biases.insert(std::stof(bias));
}

}

