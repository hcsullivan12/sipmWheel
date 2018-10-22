//
// File: main.cxx
//
// Author: Hunter Sullivan
//
// Decription: Analysis code for sipmwheel
//             This analyzes waveforms seen by each sipm for each "trigger"
// 
// Future work: Implent this into the live event display.
//

// Includes
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <dirent.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <experimental/filesystem>
#include <map>
#include <list>
#include "TFile.h"
#include "TF2.h"
#include "TH2.h"
#include "TGraphPolar.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "Utilities.h"
#include "FileReader.h"
#include "Analyzer.h"
#include "Characterizer.h"

// Preprocessing variables
#ifdef VERSION
  #define sipmWheelVersion VERSION
#endif

// Prototypes
void PrintTheFiles(const wheel::SiPMToFilesMap& sipmToFilesMap);
void PrintTheFiles(const wheel::SiPMToBiasTriggerMap& sipmToBiasTriggerMap);
void GetTheCharacterizationFiles(wheel::SiPMToBiasTriggerMap& map, const wheel::Configuration& config);
void GetTheFiles(wheel::SiPMToFilesMap& map, const wheel::Configuration& config);
void ReadConfigFile(wheel::Configuration& config);
void SaveWaveforms(wheel::FileReader& fr, const wheel::Configuration& config);
void Characterize(const wheel::Configuration& myConfig);
void Reco(const wheel::Configuration& myConfig);
void RecordBiases(wheel::Configuration& config, const std::string& value);
void RecordGains(std::map<unsigned, float>& map, const wheel::Configuration& config, const std::string& value);
void FillSiPMInfo(wheel::SiPMInfoMap& sipmInfoMap, const wheel::Configuration& config);
void OutputConfigInfo(wheel::Configuration& config);
void SaveCharacterizationPlots(std::map<unsigned, std::vector<TH1D>>         ampDists,
                               std::map<unsigned, std::vector<TGraphErrors>> ampPeaks,
                               const wheel::Configuration&                   config);

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

  //const std::string process = argv[2];
  if      (myConfig.process == "characterize") Characterize(myConfig);
  else if (myConfig.process == "reco")         Reco(myConfig);
  else    { std::cout << "Error. Must choose characterize or reco.\n" << myConfig.process << std::endl; exit(1); }

  return 0;
}

void OutputConfigInfo(wheel::Configuration& config)
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
            << "RawWaveformsFile    " << config.rawWaveformPath                << std::endl
            << "ModWaveformsFile    " << config.modWaveformPath                << std::endl
            << "RecoOutputFile      " << config.recoOutputFile                 << std::endl
            << "BaselineSubtract    " << config.baselineSubtract               << std::endl
            << "SMARange            " << config.smaRange                       << std::endl
            << "WaveformResolution  " << config.resolution                     << std::endl
            << "HitSigmaThreshold   " << config.hitSigma                       << std::endl
            << "HitFinderSearch     " << config.hitFinderSearch                << std::endl
            << "MinimumHitAmp       " << config.minimumHitAmp                  << std::endl
            << "NumerOfSiPMs        " << config.nSiPMs                         << std::endl
            << "Biases              ";for(const auto& bias : config.biases)     std::cout << bias        << "  ";
                                                                                std::cout << std::endl; std::cout 
            << "Gains               ";for(const auto& sipm : config.gains)      std::cout << sipm.second << "  ";
                                                                                std::cout << std::endl; std::cout  
            << "Breakdowns          ";for(const auto& sipm : config.breakdowns) std::cout << sipm.second << "  ";
                                                                                std::cout << std::endl; std::cout 
            << "ThetaBinSize        " << config.thetaBinSize                   << std::endl
            << "RadiusBinSize       " << config.radiusBinSize                  << std::endl
            << "AttenLengthBinSize  " << config.attenuationLengthBinSize       << std::endl
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
            << "OutputFile          " << config.characterizeOutputFile         << std::endl
            << "NFiles/SiPM         " << config.nFilesCharacterize             << std::endl
            << "SaveRawWaveforms    " << config.saveRawWaveforms               << std::endl
            << "SaveModWaveforms    " << config.saveModWaveforms               << std::endl
            << "RawWaveformsFile    " << config.rawWaveformPath                << std::endl
            << "ModWaveformsFile    " << config.modWaveformPath                << std::endl
            << "BaselineSubtract    " << config.baselineSubtract               << std::endl
            << "SMARange            " << config.smaRange                       << std::endl
            << "WaveformResolution  " << config.resolution                     << std::endl
            << "HitSigmaThreshold   " << config.hitSigma                       << std::endl
            << "HitFinderSearch     " << config.hitFinderSearch                << std::endl
            << "MinimumHitAmp       " << config.minimumHitAmp                  << std::endl
            << "NumerOfSiPMs        " << config.nSiPMs                         << std::endl
            << "AmpThreshold        " << config.characterizeAmpThr             << std::endl
            << "AmpSigma            " << config.characterizeAmpSig             << std::endl
            << "AmpFitRange         " << config.characterizeAmpFitRange        << std::endl;
  std::cout << std::setfill('-') << std::setw(80) << "-" << std::setfill(' ') << std::endl;	
  std::cout << std::endl;
}

void Reco(const wheel::Configuration& myConfig)
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
  GetTheFiles(sipmToFilesMap, myConfig);

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
    // Option to output a few waveforms
    if (trigger == 1) SaveWaveforms(fr, myConfig);
  
    wheel::Analyzer analyzer;
    analyzer.Reconstruct(sipmToTriggerMap, sipmInfoMap, myConfig, trigger);
  }
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
  for (auto& sipm : sipmToBiasTriggerMap)
  {
    wheel::SiPMToTriggerMap sipmToTriggerMap;
    wheel::FileReader fr;
    fr.ReadFiles(sipmToTriggerMap, sipm.second, sipm.first, myConfig);
    // Option to output graphs
    if (sipm.first == 1) SaveWaveforms(fr, myConfig);
    // Characterize
    wheel::SiPMInfoMap sipmInfoMap;
    ch.Characterize(sipmInfoMap, sipmToTriggerMap, myConfig);
  }
  // Save the plots
  SaveCharacterizationPlots(ch.GetAmpDists(), ch.GetAmpPeaks(), myConfig);

  clock_t end = clock();
  double duration = ((double) (end - start)) / CLOCKS_PER_SEC;
  std::cout << "Run time of " << duration << " s" << std::endl;
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

    if      (header == "pathToData")              config.pathToData      = value;
    else if (header == "rawWaveformPath")         config.rawWaveformPath = value;
    else if (header == "modWaveformPath")         config.modWaveformPath = value;
    else if (header == "recoOutputFile")          config.recoOutputFile  = value;
    else if (header == "printFiles")              value == "true" ? config.printFiles       = true : config.printFiles       = false;
    else if (header == "baselineSubtract")        value == "true" ? config.baselineSubtract = true : config.baselineSubtract = false;
    else if (header == "saveRawWaveforms")        value == "true" ? config.saveRawWaveforms = true : config.saveRawWaveforms = false;
    else if (header == "saveModWaveforms")        value == "true" ? config.saveModWaveforms = true : config.saveModWaveforms = false;
    else if (header == "nSiPMs")                  config.nSiPMs          = std::stoi(value);
    else if (header == "smaRange")                config.smaRange        = std::stoi(value);
    else if (header == "hitFinderSearch")         config.hitFinderSearch = std::stoi(value);
    else if (header == "resolution")              config.resolution      = std::stoi(value);
    else if (header == "hitSigma")                config.hitSigma        = std::stof(value);
    else if (header == "minimumHitAmp")           config.minimumHitAmp   = std::stof(value);
    else if (header == "process")                 value == "characterize" ? config.process  = "characterize" : config.process = "reco"; // default to reco
    else if (header == "characterizeAmpThr")      config.characterizeAmpThr      = std::stof(value);
    else if (header == "characterizeAmpSig")      config.characterizeAmpSig      = std::stof(value);
    else if (header == "characterizeAmpFitRange") config.characterizeAmpFitRange = std::stof(value);
    else if (header == "characterizeOutputFile")  config.characterizeOutputFile  = value;
    else if (header == "nFilesCharacterize")      config.nFilesCharacterize      = std::stoi(value);
    else if (header == "nBiases")                 config.nBiases                 = std::stoi(value);
    else if (header == "biases")                  RecordBiases(config, value); 
    else if (header == "gains")                   RecordGains(config.gains, config, value);
    else if (header == "breakdowns")              RecordGains(config.breakdowns, config, value);
    else if (header == "thetaBinSize")            config.thetaBinSize  = std::stof(value);
    else if (header == "radiusBinSize")           config.radiusBinSize = std::stof(value);
    else if (header == "attenuationLengthBinSize") config.attenuationLengthBinSize = std::stof(value);
    else if (header == "diskRadius")              config.diskRadius = std::stof(value);
    else if (header == "attenuationLength")       config.attenuationLength = std::stof(value); 
    else if (header == "nVoxels")                 config.nVoxels    = std::stoi(value);
    else if (header == "maxIterations")           config.maxIterations = std::stoi(value);
    else    { std::cout << "Cannot identify " << header << std::endl; exit(1); }
  }
  
  // Place a series of safety nets here;
  if (config.nBiases != config.biases.size())                                { std::cout << "Error. Number of biases does not match in config.\n" << std::endl; exit(1); }
  if (config.process == "reco" && config.nSiPMs != config.gains.size())      { std::cout << "Error. Number of gains does not match number of SiPMs\n." << std::endl; exit(1); }
  if (config.process == "reco" && config.nSiPMs != config.breakdowns.size()) { std::cout << "Error. Number of breakdowns does not match number of SiPMs\n." << std::endl; exit(1); }
  if (config.process == "reco" && config.biases.size() != 1)                 { std::cout << "Error. Can only specify one bias for reconstruction\n." << std::endl; exit(1); }
}

void RecordGains(std::map<unsigned, float>& map, const wheel::Configuration& config, const std::string& value)
{
  // Make sure n sipms has been defined
  if (!config.nSiPMs) { std::cout << "Error. Please specify number of sipms before listing gains and breakdowns.\n" << std::endl; exit(1); }
  
  std::stringstream linestream(value);
  std::string gain;
  unsigned sipm = 1;
  while(std::getline(linestream, gain, ',')) { map.emplace(sipm, std::stof(gain)); sipm++; }
}


void RecordBiases(wheel::Configuration& config, const std::string& value)
{
  // Make sure n biases has been defined
  if (!config.nBiases) { std::cout << "Error. Please specify number of biases before listing biases.\n" << std::endl; exit(1); }
  
  std::stringstream linestream(value);
  std::string bias;
  while(std::getline(linestream, bias, ',')) config.biases.insert(std::stof(bias));
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

void GetTheFiles(wheel::SiPMToFilesMap& map, const wheel::Configuration& config)
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

void SaveCharacterizationPlots(std::map<unsigned, std::vector<TH1D>>         ampDists, 
                               std::map<unsigned, std::vector<TGraphErrors>> ampPeaks, 
                               const wheel::Configuration&                   config)
{
  // Canvases for out histos
  TCanvas masterAmpDist("masterAmpDist",  "All Amplitude Distributions", 1000, 1000);
  TCanvas masterGainPlot("masterGainPlot", "All Gains", 1000, 1000); 

  // Divide the canvases
  masterAmpDist.Divide(config.nBiases,config.nSiPMs);
  masterGainPlot.Divide(config.nBiases,config.nSiPMs);

  // Loop over sipms
  unsigned ampCounter(1);
  unsigned gainCounter(1);
  for (unsigned sipm = 1; sipm <= config.nSiPMs; sipm++)
  {
    // Get the dists for this sipm
    //std::vector<TH1D>& dists = ampDists.find(sipm)->second;
    // Draw the amp dists and gain plot
    for (auto& dist : ampDists.find(sipm)->second)
    {
      // Amplitude dist
      masterAmpDist.cd(ampCounter);
      gStyle->SetOptStat(0);
      dist.GetXaxis()->SetTitle("Integral/a.u.");
      dist.Draw();
      masterAmpDist.Update();
      dist.GetXaxis()->SetTitle("Area/a.u."); 
      masterAmpDist.Modified();
      ampCounter++;
    }
    for (auto& peaks : ampPeaks.find(sipm)->second)
    {
      // Gain
      masterGainPlot.cd(gainCounter);
      peaks.SetMarkerStyle(20);
      peaks.SetMarkerColor(1);
      peaks.Draw("AP");
      gainCounter++;
      peaks.GetXaxis()->SetLimits(0, 5);
      peaks.SetMinimum(0);      
    }
  }

  TFile f(config.characterizeOutputFile.c_str(), "RECREATE");
  masterAmpDist.Write();
  masterGainPlot.Write();
  f.Close();
};

void SaveWaveforms(wheel::FileReader& fr, const wheel::Configuration& config)
{
  // First output raw
  if (config.saveRawWaveforms)
  {
    std::cout << "Outputing waveforms... " << std::endl;
    // Create a file to write to
    TFile f(config.rawWaveformPath.c_str(), "RECREATE");
  
    unsigned counter = 0;
    for (auto& g : fr.GetRawGraphs())
    {
      TCanvas c(std::to_string(counter).c_str(), std::to_string(counter).c_str(), 500, 500);
      //g->GetXaxis()->SetRangeUser(1900, 2100);
      //g->Fit(fits[counter], "R+", "", 1988, 2100);
      g.Draw("");
      // Write this to file
      c.Write();
      counter++;
    }
  f.Close();
  }
  // Now output modified
  if (config.saveModWaveforms)
  {
    std::cout << "Outputing modified waveforms... " << std::endl;
    // Create a file to write to
    TFile f(config.modWaveformPath.c_str(), "RECREATE");
  
    unsigned counter = 0;
    for (auto& g : fr.GetModGraphs())
    {
      TCanvas c(std::to_string(counter).c_str(), std::to_string(counter).c_str(), 500, 500);
      //g->GetXaxis()->SetRangeUser(1900, 2100);
      //g->Fit(fits[counter], "R+", "", 1988, 2100);
      g.Draw("");
      for (auto& m : fr.GetMarkers()[counter])
      {
        m.first.Draw("same");
        m.second.Draw("same");
      }
      // Write this to file
      c.Write();
      counter++;
    }
  f.Close();
  }
}
