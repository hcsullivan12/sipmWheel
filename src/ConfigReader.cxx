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

    if      (header == "process")                 config.process = value; 
    else if (header == "simulateOutputPath")      config.recoOutputPath = value;
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
    else if (header == "pathToData")              config.pathToData      = value;
    else if (header == "rawWaveformPath")         config.rawWaveformPath = value;
    else if (header == "modWaveformPath")         config.modWaveformPath = value;
    else if (header == "recoOutputPath")          config.recoOutputPath  = value;
    else if (header == "printFiles")              value == "true" ? config.printFiles       = true : config.printFiles       = false;
    else if (header == "baselineSubtract")        value == "true" ? config.baselineSubtract = true : config.baselineSubtract = false;
    else if (header == "saveRawWaveforms")        value == "true" ? config.saveRawWaveforms = true : config.saveRawWaveforms = false;
    else if (header == "saveModWaveforms")        value == "true" ? config.saveModWaveforms = true : config.saveModWaveforms = false;
    else if (header == "smoothWaveform")          value == "true" ? config.smoothWaveform = true : config.smoothWaveform = false;
    else if (header == "smaRange")                config.smaRange        = std::stoi(value);
    else if (header == "hitFinderSearch")         config.hitFinderSearch = std::stoi(value);
    else if (header == "resolution")              config.resolution      = std::stoi(value);
    else if (header == "minimumHitAmp")           config.minimumHitAmp   = std::stof(value);
    else if (header == "characterizeAmpThr")      config.characterizeAmpThr      = std::stof(value);
    else if (header == "characterizeAmpSig")      config.characterizeAmpSig      = std::stof(value);
    else if (header == "characterizeAmpFitRange") config.characterizeAmpFitRange = std::stof(value);
    else if (header == "characterizeOutputPath")  config.characterizeOutputPath  = value;
    else if (header == "nFilesCharacterize")      config.nFilesCharacterize      = std::stoi(value);
    else if (header == "nBiases")                 config.nBiases                 = std::stoi(value);
    else if (header == "biases")                  RecordBiases(config, value); 
    else if (header == "gains")                   RecordGains(config.gains, config, value);
    else if (header == "breakdowns")              RecordGains(config.breakdowns, config, value);
    else if (header == "maxIterations")           config.maxIterations = std::stoi(value);
    else    { std::cout << "Cannot identify " << header << " listed in configuration file." << std::endl; exit(1); }
  }
  
  // Place a series of safety nets here;
  if (config.process == "characterize" && config.nBiases != config.biases.size()) { std::cout << "Error. Number of biases does not match in config.\n" << std::endl; exit(1); }
  if (config.process == "reco" && config.nBiases != config.biases.size()) { std::cout << "Error. Number of biases does not match in config.\n" << std::endl; exit(1); }
  if (config.process == "reco" && config.nSiPMs != config.gains.size())      { std::cout << "Error. Number of gains does not match number of SiPMs\n." << std::endl; exit(1); }
  if (config.process == "reco" && config.nSiPMs != config.breakdowns.size()) { std::cout << "Error. Number of breakdowns does not match number of SiPMs\n." << std::endl; exit(1); }
  if (config.process == "reco" && config.biases.size() != 1)                 { std::cout << "Error. Can only specify one bias for reconstruction\n." << std::endl; exit(1); }
  if (config.nSiPMs <= 1)              { std::cout << "Error. Must specify an even number of sipms!\n" << std::endl; exit(1); }
  if (config.diskRadius <= 0)          { std::cout << "Error. Disk radius == 0!\n" << std::endl; exit(1); }
  if (config.process == "simulate" && config.diskThickness <= 0)       { std::cout << "Error. Disk thickness == 0!\n" << std::endl; exit(1); }
  if (config.process == "simulate" && config.tpbEmissionPeak <= 0)     { std::cout << "Error. TPB emission peak == 0!\n" << std::endl; exit(1); }
  if (config.process == "simulate" && config.indexRefractionDisk <= 0) { std::cout << "Error. Index of refraction disk == 0!\n" << std::endl; exit(1); }
  if (config.process == "simulate" && config.indexRefractionEnv <= 0)  { std::cout << "Error. Index of refraction medium == 0!\n" << std::endl; exit(1); }
  if (config.process == "simulate" && config.permittivityDisk <= 0)    { std::cout << "Error. Permittivity of disk <= 0!\n" << std::endl; exit(1); }
  if (config.process == "simulate" && config.permittivityEnv <= 0)     { std::cout << "Error. Permittivity of environment <= 0!\n" << std::endl; exit(1); }
  if (config.process == "simulate" && config.sourcePosition.size() < 2) { std::cout << "Error. Please list sourcePosition as r,theta in configuration\n" << std::endl; exit(1); }
  if (config.process == "simulate" && config.sourceSigma == 0)         { std::cout << "Error. Source sigma == 0!\n" << std::endl; exit(1); }
  if (config.process == "simulate" && config.smearSigma == 0)          { std::cout << "Error. Smear sigma == 0!\n" << std::endl; exit(1); }
  if (config.process == "simulate" && config.nPhotonsToLaunch <= 0) { std::cout << "Error. Please specify number of photons to simulate.\n" << std::endl; exit(1); }
  if (config.process == "simulate" && config.bulkAttenuation <= 0)  { std::cout << "Error. Please specify bulkAttenuation > 0\n" << std::endl; exit(1); }
  if (config.process == "simulate" && config.bulkAbsorption < 0) { std::cout << "Error. Please specify a bulkAbsorption > 0!\n" << std::endl; exit(1); }
  if (config.process == "simulate" && (config.surfaceAbsorptionCoeff < 0 || config.surfaceAbsorptionCoeff > 1)) { std::cout << "Error. Please specify a 0 < surfaceAbsorptionCoeff < 1.\n" << std::endl; exit(1); }
  if (config.process == "simulate" && config.terminationThreshold <= 0) { std::cout << "Error. Please specify termination threshold > 0.\n" << std::endl; exit(1); } 
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
  if (config.process == "reco") 
  {
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
  else if (config.process == "characterize")
  {
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
  else
  {
    std::cout << "Process                " << config.process                        << std::endl
              << "SimulateOutputPath     " << config.recoOutputPath                 << std::endl
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
              << "attenuationLength      " << config.attenuationLength              << std::endl;
    std::cout << std::setfill('-') << std::setw(80) << "-" << std::setfill(' ') << std::endl;
    std::cout << std::endl;  
  }
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

void ConfigReader::RecordPosition(Configuration& config, const std::string& value)
{ 
  std::stringstream linestream(value);
  std::string rTheta;
  while(std::getline(linestream, rTheta, ',')) config.sourcePosition.push_back(std::stof(rTheta));
}

}

