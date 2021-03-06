// File: FileReader.cxx
//
// Author: Hunter Sullivan
//
// Discription: Structure to read the input waveform files.
//

#include "FileReader.h"
#include "TCanvas.h"
#include "TFile.h"

namespace wheel {

FileReader::FileReader()
{}

FileReader::~FileReader()
{}

/*
 * This method is for characterization
 * 
 **/
void FileReader::ReadFiles(SiPMToTriggerMap&    sipmToTriggerMap, 
                           const BiasToFileMap& biasMap, 
                           const unsigned&      sipm, 
                           const Configuration& config)
{
  // Create a trigger list for this sipm
  // If this is characterization, this is 
  // just a list of the different biases
  // Place a safety net here
  if (config.process == "characterize") 
  {
    if (config.biases.size() != biasMap.size()) { std::cout << "Error. The files do not match the biases listed in config.\n"; std::cout << std::endl; exit(1); } 
  }

  std::cout << "\nFinding hits for SiPM " << sipm << "..." << std::endl;
  std::vector<HitCandidateVec> triggerList;
  triggerList.reserve(biasMap.begin()->second.size());

  // Loop over the biases for this sipm
  for (const auto& bias : biasMap)
  {
    std::cout << "Bias = " << bias.first << "\n";
    unsigned counter = 0;
    // Loop over the files for this bias and sipm
    for (const auto& file : bias.second)
    {
      // Just to make sure we didnt get stuck
      if (counter % 100 == 0) std::cout << "Trigger #" << counter << std::endl;
      counter++;
      // Option to only do so many
      if (counter > config.nFilesCharacterize) break;
      // Create a temp hitVec
      HitCandidateVec hitCandVec;
      // Now read this file
      ReadFile(hitCandVec, file, bias.first, sipm, counter, config);
      // Safety net to protect against division settings
      triggerList.emplace_back(hitCandVec);
    }
  }
  // There will only be one sipm here
  sipmToTriggerMap.emplace(sipm, triggerList);
}

/*
 * This method is for reco
 *
 */
void FileReader::ReadFiles(SiPMToTriggerMap& sipmToTriggerMap, 
                           const SiPMToFilesMap& map, 
                           const unsigned& trigger, 
                           const SiPMInfoMap& sipmInfoMap, 
                           const Configuration& config)
{
  // Loop over the sipms
  for (const auto& sipm : map)
  {
    // Create a hitCandVec for this sipm
    HitCandidateVec hitCandVec;
    std::vector<HitCandidateVec> triggerList;
    // There should only be one bias specified
    // Now read this file
    const std::string& file = *std::next(map.find(sipm.first)->second.begin(), trigger - 1);
    ReadFile(hitCandVec, file, sipmInfoMap.find(sipm.first)->second.bias, sipm.first, trigger, config);
    // Safety net to protect against division settings
    triggerList.emplace_back(hitCandVec);
    
    sipmToTriggerMap.emplace(sipm.first, triggerList);
  }
}

/*
 *  The method extracts the waveform from the data file
 *
 */
void FileReader::ReadFile(HitCandidateVec& hitCandidateVec, const std::string& filename, const float& bias, const unsigned& channel, const unsigned& trigger, const Configuration& config)
{
  // Open the file
  std::ifstream file(filename.c_str());
  if (!file.is_open()) 
  {
    std::cout << "Cannot open file: " << filename << std::endl;
    return;
  }
  std::string word = "";
  std::vector<float> signal;

  // Skip first text in file
  while (word != "Time,Ampl") file >> word;
  
  std::getline(file, word); std::getline(file, word);
  int counter = 0;
  int push = 0;
  
  while(!file.eof()) 
  {
    // Data starts
    std::string yTemp;
    std::string xTemp;
    std::getline(file, xTemp, ',');
    std::getline(file, yTemp);
    counter++;

    signal.push_back( -1*atof(yTemp.c_str()) );
  }
  
  // Analyze this waveform
  Analyze(signal, hitCandidateVec, filename, bias, channel, trigger, config);
}

/*
 *  This method analyzes each waveform:
 *    1) Smoothing alg
 *    2) HitFinder alg
 *
 */
void FileReader::Analyze(std::vector<float>&  signal, 
                         HitCandidateVec&     hitCandidateVec, 
                         const std::string&   filename, 
                         const float&         bias, 
                         const unsigned&      channel, 
                         const unsigned&      trigger,
                         const Configuration& config)
{
  // Raw waveform: Only allow to store a few
  if (config.process == "reco" && config.saveRawWaveforms)
  {
    TGraph g(signal.size());
    std::string name = "Trigger" + std::to_string(trigger) + "_Channel" + std::to_string(channel);
    g.SetNameTitle(name.c_str(), name.c_str());

    unsigned tick(0);
    for (const auto& amp : signal)
    {
      g.SetPoint(tick, tick, amp);
      tick++;
    }
 
    TCanvas c(name.c_str(), name.c_str(), 500, 500);
    g.Draw();   

    // Add to our output file
    TFile f(config.rawWaveformPath.c_str(), "UPDATE");
    g.Write();
    f.Close();
  }

  // Smooth the waveform
  WaveformAlg waveformAlg; 
  if (config.smoothWaveform) waveformAlg.SmoothWaveform(signal, config); 

  // Let's do the hit finding now
  HitCandidateVec      hitCandVec;
  MergeHitCandidateVec mergedHitsVec;
  waveformAlg.FindHitCandidates(signal, 0, channel, bias, hitCandVec, config);
  
  // Append the hits
  hitCandidateVec.insert(hitCandidateVec.end(), hitCandVec.begin(), hitCandVec.end());

  // Make the markers for the graphs
  if (config.process == "reco" && config.saveModWaveforms) MakeTheMarkers(signal, channel, trigger, hitCandVec, config);
}

void FileReader::MakeTheMarkers(const std::vector<float>& signal, const unsigned& channel, const unsigned& trigger, const HitCandidateVec& hitCandVec, const Configuration& config)
{
  // Make the markers first
  std::list<TMarker> mks;
  for (const auto& hit : hitCandVec)
  {
    //std::cout << "Hit start = " << hit.startTickAmp << "  hit stop = " << hit.stopTickAmp << std::endl;
    TMarker mMax(hit.stopTick, hit.stopTickAmp, 23);
    TMarker mMin(hit.startTick, hit.startTickAmp, 23);
    TMarker mPeak(hit.hitPeakTick, hit.hitPeak, 23);
    mMax.SetMarkerColor(4);
    mMax.SetMarkerSize(1.5);
    mMin.SetMarkerColor(4);
    mMin.SetMarkerSize(1.5); 
    mPeak.SetMarkerColor(2);
    mPeak.SetMarkerSize(1.5);
    mks.emplace_back(mMax);
    mks.emplace_back(mMin);
    mks.emplace_back(mPeak);
  } 
  
  // Now write the waveforms
  TGraph g(signal.size());
  std::string name = "Trigger" + std::to_string(trigger) + "_Channel" + std::to_string(channel);
  g.SetNameTitle(name.c_str(), name.c_str());
  unsigned tick(0);
  for (const auto& amp : signal)
  {
    g.SetPoint(tick, tick, amp);
    tick++;
  }

  TCanvas c(name.c_str(), name.c_str(), 500, 500);
  g.Draw("");
  for (auto& m : mks) m.Draw("same");
    
  // Add to our output file
  TFile f(config.modWaveformPath.c_str(), "UPDATE");
  c.Write();
  f.Close();
}
}

