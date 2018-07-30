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

// Prototypes
void PrintTheFiles(const wheel::SiPMToFilesMap& sipmToFilesMap);
void PrintTheFiles(const wheel::SiPMToBiasTriggerMap& sipmToBiasTriggerMap);
void GetTheCharacterizationFiles(wheel::SiPMToBiasTriggerMap& map, const wheel::Configuration& config);
void GetTheFiles(wheel::SiPMToFilesMap& map, const wheel::Configuration& config);
void ReadConfigFile(wheel::Configuration& config);
void SaveWaveforms(std::vector<TGraph>& graphs, std::vector<std::vector<TMarker>>& markers, const wheel::Configuration& config);
void Characterize(const wheel::Configuration& myConfig);
void Reco(const wheel::Configuration& myConfig);
void RecordBiases(wheel::Configuration& config, const std::string& value);
void RecordGains(std::map<unsigned, float>& map, const wheel::Configuration& config, const std::string& value);
void FillSiPMInfo(wheel::SiPMInfoMap& sipmInfoMap, const wheel::Configuration& config);
void MakeRecoPlots(wheel::Analyzer& analyzer, const wheel::Configuration& config, const unsigned& trigger);
void OutputConfigInfo(wheel::Configuration& config);
void SaveCharacterizationPlots(std::multimap<unsigned, std::vector<TH1D>>         ampDists,
                               std::multimap<unsigned, std::vector<TGraphErrors>> ampPeaks,
                               const wheel::Configuration&                        config);

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
  std::cout << "         Author: Hunter Sullivan (UT Arlington)         "      << std::endl;
  std::cout                                                                    << std::endl;
  std::cout << "SiPM Wheel Configuration:\n";
  if (config.process == "reco") {
  std::cout << "Process             " << config.process                        << std::endl
            << "PathToData          " << config.pathToData                     << std::endl
            << "SaveWaveforms       " << config.saveWaveforms                  << std::endl
            << "WaveformsFile       " << config.outputPath                     << std::endl
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
            << "DiskRadius          " << config.diskRadius                     << std::endl;
  std::cout << std::setfill('-') << std::setw(80) << "-" << std::setfill(' ') << std::endl;
  std::cout << std::endl;
  return;
  }
  std::cout << "Process             " << config.process                        << std::endl
            << "PathToData          " << config.pathToData                     << std::endl
            << "OutputFile          " << config.characterizeOutputFile         << std::endl
            << "SaveWaveforms       " << config.saveWaveforms                  << std::endl
            << "WaveformsFile       " << config.outputPath                     << std::endl
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
  // First we need to get the files to read data from
  std::cout << "\nGetting the files from: " << myConfig.pathToData << std::endl;
  wheel::SiPMToFilesMap sipmToFilesMap;
  GetTheFiles(sipmToFilesMap, myConfig);

  // Option to print the collected files
  if (myConfig.printFiles) PrintTheFiles(sipmToFilesMap);

  // Now read the files
  std::cout << "\nReading the files... " << std::endl;
  wheel::SiPMToTriggerMap sipmToTriggerMap;
  wheel::FileReader fr;
  fr.ReadFiles(sipmToTriggerMap, sipmToFilesMap, myConfig);
 
  // Option to output graphs
  if (myConfig.saveWaveforms) SaveWaveforms(fr.GetGraphs(), fr.GetMarkers(), myConfig);

  // Fill the sipm info first
  wheel::SiPMInfoMap sipmInfoMap;
  FillSiPMInfo(sipmInfoMap, myConfig);

  // Start the main work horse after filling sipm info
  // To reduce overhead, run reco and make plots before looping
  // Make sure there is the same amount of data for each sipm
  const unsigned& nFiles = sipmToTriggerMap.find(1)->second.size();
  for (const auto& sipm : sipmToTriggerMap)
  {
    if (sipm.second.size() != nFiles) { std::cout << "Error. Different data sizes for SiPMs\n" << std::endl; exit(1); }
  }
  // Loop over the triggers
  for (unsigned trigger = 1; trigger <= sipmToTriggerMap.find(1)->second.size(); trigger++)
  {
    std::cout << std::setfill('*') << std::setw(80) << "-" << std::setfill(' ')  << std::endl; 
    std::cout << "\nInitializing event " << trigger << "..." << std::endl;
    wheel::Analyzer analyzer;
    analyzer.RunReco(sipmToTriggerMap, sipmInfoMap, myConfig, trigger);
    // Now make the plots
    MakeRecoPlots(analyzer, myConfig, trigger);
  }
}

void Characterize(const wheel::Configuration& myConfig)
{
  // First we need to get the files to read data from
  std::cout << "Getting the files from: " << myConfig.pathToData << std::endl;
  wheel::SiPMToBiasTriggerMap sipmToBiasTriggerMap;
  GetTheCharacterizationFiles(sipmToBiasTriggerMap, myConfig);

  // Option to print the collected files
  if (myConfig.printFiles) PrintTheFiles(sipmToBiasTriggerMap);

  // Now read the files
  std::cout << "Reading the files... " << std::endl;
  wheel::SiPMToTriggerMap sipmToTriggerMap;
  wheel::FileReader fr;
  fr.ReadFiles(sipmToTriggerMap, sipmToBiasTriggerMap, myConfig);
 
  // Option to output graphs
  if (myConfig.saveWaveforms) SaveWaveforms(fr.GetGraphs(), fr.GetMarkers(), myConfig);

  wheel::Characterizer ch;
  wheel::SiPMInfoMap sipmInfoMap;
  ch.Characterize(sipmInfoMap, sipmToTriggerMap, myConfig);

  // Save the plots
  SaveCharacterizationPlots(ch.GetAmpDists(), ch.GetAmpPeaks(), myConfig);
}

void FillSiPMInfo(wheel::SiPMInfoMap& sipmInfoMap, const wheel::Configuration& config)
{
  // Loop over sipms
  for (unsigned sipm = 1; sipm <= config.nSiPMs; sipm++)
  {
    wheel::SiPMInfo sipmInfo;
    sipmInfo.gain      = config.gains.find(sipm)->second;
    sipmInfo.breakdown = config.breakdowns.find(sipm)->second; 

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

    if      (header == "pathToData")              config.pathToData     = value;
    else if (header == "outputPath")              config.outputPath     = value;
    else if (header == "recoOutputFile")          config.recoOutputFile = value;
    else if (header == "printFiles")              value == "true" ? config.printFiles       = true : config.printFiles       = false;
    else if (header == "baselineSubtract")        value == "true" ? config.baselineSubtract = true : config.baselineSubtract = false;
    else if (header == "saveWaveforms")           value == "true" ? config.saveWaveforms    = true : config.saveWaveforms    = false;
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
    else if (header == "nBiases")                 config.nBiases                 = std::stoi(value);
    else if (header == "biases")                  RecordBiases(config, value); 
    else if (header == "gains")                   RecordGains(config.gains, config, value);
    else if (header == "breakdowns")              RecordGains(config.breakdowns, config, value);
    else if (header == "thetaBinSize")            config.thetaBinSize  = std::stof(value);
    else if (header == "radiusBinSize")           config.radiusBinSize = std::stof(value);
    else if (header == "attenuationLengthBinSize") config.attenuationLengthBinSize = std::stof(value);
    else if (header == "diskRadius")              config.diskRadius = std::stof(value);
    else    { std::cout << "Cannot identify " << header << std::endl; exit(1); }
  }
  
  // Place a series of safety nets here;
  if (config.nBiases != config.biases.size())     { std::cout << "Error. Number of biases does not match in config.\n" << std::endl; exit(1); }
  if (config.nSiPMs  != config.gains.size())      { std::cout << "Error. Number of gains does not match number of SiPMs\n." << std::endl; exit(1); }
  if (config.nSiPMs  != config.breakdowns.size()) { std::cout << "Error. Number of breakdowns does not match number of SiPMs\n." << std::endl; exit(1); }
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
     
    // Loop over the bias directories
    for(const auto& biasDirIter : stdfs::directory_iterator(sipmDir.c_str()))
    {
      // First get the bias
      std::string bias;
      std::string path = biasDirIter.path();
      for (std::string::iterator ch = path.end(); *ch != '/'; ch--) bias = *ch + bias;
      const float thisBias = std::stof(bias);

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

void SaveCharacterizationPlots(std::multimap<unsigned, std::vector<TH1D>>         ampDists, 
                               std::multimap<unsigned, std::vector<TGraphErrors>> ampPeaks, 
                               const wheel::Configuration&                        config)
{
  // Canvases for out histos
  TCanvas masterAmpDist("masterAmpDist",  "All Amplitude Distributions", 1000, 1000);
  TCanvas masterGainPlot("masterGainPlot", "All Gains", 1000, 1000); 

  // Divide the canvases
  masterAmpDist.Divide(config.nBiases,config.nSiPMs);
  masterGainPlot.Divide(config.nBiases,config.nSiPMs);

  // Loop over sipms
  unsigned counter(1);
  for (unsigned sipm = 1; sipm <= config.nSiPMs; sipm++)
  {
    // Get the dists for this sipm
    //std::vector<TH1D>& dists = ampDists.find(sipm)->second;
    // Draw the amp dists and gain plot
    for (auto& dist : ampDists.find(sipm)->second)
    {
      // Amplitude dist
      masterAmpDist.cd(counter);
      gStyle->SetOptStat(0);
      dist.Draw("apl");      
      counter++;
    }
    counter = 1;
    for (auto& peaks : ampPeaks.find(sipm)->second)
    {
      // Gain
      masterGainPlot.cd(counter);
      peaks.SetMarkerStyle(20);
      peaks.SetMarkerColor(1);
      peaks.Draw("AP");
      counter++;
      peaks.GetXaxis()->SetLimits(0, 5);
      peaks.SetMinimum(0);      
    }
  }

  TFile f(config.characterizeOutputFile.c_str(), "RECREATE");
  masterAmpDist.Write();
  masterGainPlot.Write();
  f.Close();
};

void SaveWaveforms(std::vector<TGraph>& graphs, std::vector<std::vector<TMarker>>& markers, const wheel::Configuration& config)
{
  std::cout << "\nOutputing waveforms... \n" << std::endl;
  // Create a file to write to
  TFile f(config.outputPath.c_str(), "RECREATE");
  
  unsigned counter = 0;
  for (auto& g : graphs)
  {
    // Set a safety net here
    if (counter >= 10 || counter > markers.size()) break;

    TCanvas c(std::to_string(counter).c_str(), std::to_string(counter).c_str(), 500, 500);
    //g->GetXaxis()->SetRangeUser(1900, 2100);
    //g->Fit(fits[counter], "R+", "", 1988, 2100);
    g.Draw("");
    for (auto& m : markers[counter])
    {
      m.Draw("same");
    }
    // Write this to file
    c.Write();
    counter++;
  }
 
  f.Close();;
}

Double_t disk(Double_t* x, Double_t* par)
{
  Double_t f = x[0]*x[0] + x[1]*x[1];
  return f;
}

void MakeRecoPlots(wheel::Analyzer& analyzer, const wheel::Configuration& config, const unsigned& trigger)
{
  // Get the maps needed for the plots
  auto accumulatorMap = analyzer.GetAccumulatorMap();
  auto mleParameters  = analyzer.GetMLEParameters();

  // Make the comparison plot 
  std::string name = "data_mle_trigger" + std::to_string(trigger);
  TCanvas c1(name.c_str(), name.c_str(), 1000, 1000);
  // Get the counts lambda_m using the mlestimates for r, theta, and N0
  std::vector<unsigned> prediction;
  prediction.reserve(config.nSiPMs);
  for (int m = 1; m <= config.nSiPMs; m++) 
  {
    float lambda_m = analyzer.ComputeLambda( mleParameters[0].second.find("radius")->second, mleParameters[0].second.find("theta")->second, 
                                             mleParameters[0].first, m, mleParameters[0].second.find("attenuationLength")->second);
    prediction[m - 1] = static_cast<unsigned>(lambda_m);
  }

  name = "pred_trigger" + std::to_string(trigger);
  TH1D pred(name.c_str(), name.c_str(), config.nSiPMs, 0, config.nSiPMs);
  std::cout << std::endl;
  std::cout << "Bin comparison:\n";
  for (int posBin = 1; posBin <= config.nSiPMs; posBin++) 
  {
    std::cout << "Pred Bin " << posBin << " :  " << prediction[posBin - 1] << " p.e." << std::endl;
    pred.SetBinContent( posBin, prediction[posBin - 1] );
  }

  pred.SetFillStyle(3001);
  pred.SetFillColor(kRed);
  pred.SetLineWidth(3);
  pred.SetLineColor(kRed);
  pred.GetXaxis()->SetTitle("SiPM Position");
  pred.GetYaxis()->SetTitle("p.e");
  pred.SetTitle("Estimator for SiPM Wheel");
  gStyle->SetOptStat(0);
  pred.SetMaximum(30);
  pred.SetMinimum(0);
  pred.Draw();

  name = "dataHisto_trigger" + std::to_string(trigger);
  TH1D dataHisto(name.c_str(), name.c_str(), config.nSiPMs, 0, config.nSiPMs);
  for (int posBin = 1; posBin <= config.nSiPMs; posBin++) 
  {
    std::cout << "Data Bin " << posBin << " :  " << analyzer.GetData().find(posBin)->second << " p.e." << std::endl;
    dataHisto.SetBinContent( posBin, analyzer.GetData().find(posBin)->second );
  }
  dataHisto.SetMarkerStyle(21);
  dataHisto.SetMarkerSize(2);
  dataHisto.Draw("same P");

  TLegend leg(0.1,0.6,0.3,0.7);
  leg.AddEntry(&pred, "Estimator", "f");
  leg.AddEntry(&dataHisto, "Data", "p");
  leg.Draw("same");

  // Make the reco plot
  name = "reco_trigger" + std::to_string(trigger);
  TCanvas c2(name.c_str(), name.c_str(), 1000, 1000);
  TF2 fDisk("fDisk", disk, -2*config.diskRadius, 2*config.diskRadius, -2*config.diskRadius, 2*config.diskRadius, 0);
  fDisk.SetContour(1);
  fDisk.SetContourLevel(0, config.diskRadius*config.diskRadius);
  fDisk.SetLineColor(1);
  fDisk.SetLineWidth(4);
  //fDisk.GetHistogram()->GetXaxis()->SetTitle("x/cm");
  //fDisk.GetHistogram()->GetYaxis()->SetTitle("y/cm");
  fDisk.SetTitle("Reconstructed Position of Point Source");
  double x = mleParameters[0].second.find("radius")->second*TMath::Cos(mleParameters[0].second.find("theta")->second*(TMath::Pi()/180) );
  double y = mleParameters[0].second.find("radius")->second*TMath::Sin(mleParameters[0].second.find("theta")->second*(TMath::Pi()/180) );
  // Compute errors
  //std::vector<float> eX, eY;
  //ComputeErrors(eX, eY, mleParameters[0], config);
  TMarker point(x, y, 8);
  point.SetMarkerColor(4);
  point.SetMarkerSize(2);
  fDisk.Draw();
  fDisk.GetHistogram()->GetXaxis()->SetTitle("x/cm");
  fDisk.GetHistogram()->GetYaxis()->SetTitle("y/cm");
  c2.Modified();
  point.Draw("same");

  // Draw y = 0
  TF1 line0("line0", "0", -config.diskRadius, config.diskRadius);
  line0.SetLineStyle(3);
  line0.SetLineColor(1);
  line0.SetLineWidth(3);
  line0.Draw("same");

  // Draw y = x
  TF1 line1("line1", "x", -(sqrt(2)/2)*config.diskRadius, (sqrt(2)/2)*config.diskRadius);
  line1.SetLineStyle(3);
  line1.SetLineColor(1);
  line1.SetLineWidth(3);
  line1.Draw("same");

  // Draw y = -x
  TF1 line2("line2", "-x", -(sqrt(2)/2)*config.diskRadius, (sqrt(2)/2)*config.diskRadius);
  line2.SetLineStyle(3);
  line2.SetLineColor(1);
  line2.SetLineWidth(3);
  line2.Draw("same");

  // Draw x = 0
  TF1 line3("line3", "100000*x", -0.00001*config.diskRadius, 0.00001*config.diskRadius);
  line3.SetLineStyle(3);
  line3.SetLineColor(1);
  line3.SetLineWidth(3);
  line3.Draw("same"); 

  // Now I want to draw a likelihood heat map 
  name = "heatMap_trigger" + std::to_string(trigger);
  TCanvas c3(name.c_str(), "Likelihood Heat Map", 800, 800);
  unsigned rBins     = config.diskRadius/config.radiusBinSize;
  unsigned thetaBins = 360/config.thetaBinSize;
  
  TH2F lhHist("lhHist", "Likelihood Heat Map", thetaBins, 0, 360, rBins, 0, config.diskRadius);
  //TH2F lhHist("lhHist", "lhHist", 2*rBins, -config.diskRadius, config.diskRadius, 2*rBins, -config.diskRadius, config.diskRadius);
  // Loop through the accumulator map to find all r, theta, and likelihood 
  // associated with the mle attenuation length
  const float& attenL = mleParameters[0].second.find("attenuationLength")->second;
  for (const auto& index : accumulatorMap)
  {
    /*// Convert the r and theta to x and y
    double x = index.second.find("radius")->second*TMath::Cos(index.second.find("theta")->second*(TMath::Pi()/180) );
    double y = index.second.find("radius")->second*TMath::Sin(index.second.find("theta")->second*(TMath::Pi()/180) );

    lhHist.SetBinContent(lhHist.GetXaxis()->FindBin(x), lhHist.GetYaxis()->FindBin(y), TMath::Exp(index.second.find("likelihood")->second));*/
    /*std::cout << "Radius     = " << index.second.find("radius")->second     << "  "
              << "Theta      = " << index.second.find("theta")->second      << "  "
              << "Likelihood = " << index.second.find("likelihood")->second << std::endl;    */
    lhHist.SetBinContent(lhHist.GetXaxis()->FindBin(index.second.find("theta")->second),
                         lhHist.GetYaxis()->FindBin(index.second.find("radius")->second), 
                         TMath::Exp(index.second.find("likelihood")->second));
  }

  //lhHist.Draw("colz");

  c3.SetTheta(90.0);
  c3.SetPhi(0.0);
  c3.SetTitle("Likelihood Heat Map");
  gStyle->SetPalette(kBird);
  lhHist.Draw("lego2 pol");

  // Draw the mle for r, theta
  /*Double_t theta[1]   = {mleParameters[0].second.find("theta")->second*(TMath::Pi()/180)};
  Double_t radius[1]  = {mleParameters[0].second.find("radius")->second};
  Double_t etheta[1]  = {(config.thetaBinSize/2)*(TMath::Pi()/180)};
  Double_t eradius[1] = {config.radiusBinSize/2};*/
  
  /*std::cout << theta[0]*(180/TMath::Pi()) << " " << radius[0] << " " << etheta[0]*(180/TMath::Pi()) << " " << eradius[0] << std::endl;*/
  
  // This does not work for some reason, not plotting marker
  /* TGraphPolar p(1, theta, radius, etheta, eradius);
  p.SetMarkerStyle(8);
  p.SetMarkerSize(5);
  p.SetMarkerColor(2);
  p.SetLineColor(1);
  p.SetLineWidth(3);
  p.Draw("AOPE");
  c2.Update();
  p.SetMinRadial(0);
  p.SetMaxRadial(config.diskRadius);
  p.GetPolargram()->SetToRadian();*/

  /*TMarker m(theta[0], radius[0], 8);
  m.SetMarkerColor(2);
  m.SetMarkerSize(1.5);
  m.Draw("same");*/

  // Write these to a file
  std::cout << "\nSaving plots to reco output file... \n" << std::endl;
  TFile f(config.recoOutputFile.c_str(), "UPDATE");
  c1.Write();
  c2.Write();
  c3.Write();
}
