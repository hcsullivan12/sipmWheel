// 
// File: FileReader.cxx
//
// Author: Hunter Sullivan
//
// Discription: Structure to read the input waveform files.
//

#include "FileReader.h"

namespace wheel {

FileReader::FileReader()
{}

FileReader::~FileReader()
{}

void FileReader::ReadFiles(SiPMToTriggerMap& sipmToTriggerMap, const SiPMToBiasTriggerMap& map, const Configuration& config)
{
  // Loop over the sipms
  unsigned channel = 1;
  for (const auto& sipm : map)
  {
    // Create a trigger list for this sipm
    // If this is characterization, this is 
    // just a list of the different biases
    // Place a safety net here
    if (config.process == "characterize") 
    {
      if (config.biases.size() != sipm.second.size()) { std::cout << "Error. The files do not match the biases listed in config.\n"; std::cout << std::endl; exit(1); } 
    }
    std::vector<HitCandidateVec> triggerList;
    unsigned totalSize = 0;
    for (const auto& bias : sipm.second) totalSize += bias.second.size();
    triggerList.reserve(sipm.second.size());

    // Loop over the biases for this sipm
    for (const auto& bias : sipm.second)
    {
      // Loop over the files for this bias and sipm
      for (const auto& file : bias.second)
      {
        // Create a temp hitVec
        HitCandidateVec hitCandVec;
        // Now read this file
        
        ReadFile(hitCandVec, file, bias.first, sipm.first, config);
        // Safety net to protect against division settings
        if (hitCandVec.size() < 20) triggerList.emplace_back(hitCandVec);
        triggerList.emplace_back(hitCandVec);
      }
    }
    sipmToTriggerMap.emplace(sipm.first, triggerList);
    channel++;
  }

}

void FileReader::ReadFiles(SiPMToTriggerMap& sipmToTriggerMap, const SiPMToFilesMap& map, const Configuration& config)
{
  // Loop over the sipms
  unsigned channel = 1;
  for (const auto& sipm : map)
  {
    // Create a trigger list for this sipm
    std::vector<HitCandidateVec> triggerList;
    triggerList.reserve(sipm.second.size());

    // Since the files and biases are ordered
    // should be able to iterate through bias set
    std::set<float>::iterator bias = config.biases.begin();
    // Loop over the files for this sipm
    for (const auto& file : sipm.second)
    {
      // Create a temp hitVec
      HitCandidateVec hitCandVec;
      // Now read this file
      ReadFile(hitCandVec, file, *bias, channel, config);
      // Safety net to protect against division settings
      if (hitCandVec.size() < 20) triggerList.emplace_back(hitCandVec);
      //bias++;
    }
    sipmToTriggerMap.emplace(channel, triggerList);
    channel++;
  }
}

void FileReader::ReadFile(HitCandidateVec& hitCandidateVec, const std::string& filename, const float& bias, const unsigned& channel, const Configuration& config)
{
  // Open the file
  std::ifstream file(filename.c_str());
  if (!file.is_open()) {
    std::cout << "Cannot open file: " << filename << std::endl;
    return;
  }

  std::string word = "";
  std::vector<float> signal;

  // Skip first text
  while (word != "Time,Ampl") {
    file >> word;
  }

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

    signal.push_back( atof(yTemp.c_str()) );
  }
 
  WaveformAlg waveformAlg; 
  waveformAlg.SmoothWaveform2(signal, config);
  std::vector<float> waveform;
  counter = 0;
  for(const auto& amp : signal)
  {
    if (counter % config.resolution != 0) 
    {
      counter++;
      continue;
    }
    counter++;
    waveform.push_back(amp);
  }

  // Let's do the hit finding now
  HitCandidateVec hitCandVec;
  waveformAlg.FindHitCandidates(waveform, 0, channel, bias, hitCandVec, config);
  // Safety net to protect against division settings
  if (hitCandVec.size() > 20) return;


  // Only allow to store a few
  if (config.saveWaveforms && waveforms.size() < 10)
  {
    TGraph g(waveform.size());
    g.SetNameTitle(filename.c_str(), filename.c_str());
    unsigned sample = 1;
    for (const auto& amp : waveform)
    {
      g.SetPoint(sample, sample - 1, amp);
      sample++;
    }
    waveforms.push_back(g);
  }

  std::cout << "Found: " << hitCandVec.size() << " hits. " << std::endl;
  // Append
  hitCandidateVec.insert(hitCandidateVec.end(), hitCandVec.begin(), hitCandVec.end());

  if (config.saveWaveforms && waveforms.size() <= 10) MakeTheMarkers(hitCandVec);

  //fits.push_back(FitTheHits(waveform, hitCandVec));

  /*for (const auto& hit : hitCandVec) 
  {
    std::cout << "\nHIT...\n";
    std::cout << "Start tick: " << hit.startTick  << std::endl
	      << "Stop  tick: " << hit.stopTick   << std::endl
	      << "Centertick: " << hit.hitCenter  << std::endl
	      << "Hit height: " << hit.hitHeight  << std::endl;
  } */
}

void FileReader::MakeTheMarkers(const HitCandidateVec& hitCandVec)
{
  std::vector<TMarker> mks;
  for (const auto& hit : hitCandVec)
  {
    TMarker m(hit.hitCenter, hit.hitHeight, 8);
    m.SetMarkerColor(2);
    m.SetMarkerSize(1.5);
    mks.push_back(m);
  } 
  markers.push_back(mks);
}
}

