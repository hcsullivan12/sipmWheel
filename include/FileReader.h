// 
// File: FileReader.h
//
// Author: Hunter Sullivan
//
// Discription: Structure to read the input waveform files.
//

#ifndef FILEREADER_H
#define FILEREADER_H

#include "Utilities.h"
#include "WaveformAlg.h"
#include <fstream>
#include <iostream>
#include "TGraph.h"
#include "TMarker.h"

namespace wheel {

class FileReader {

public:
  FileReader();
  ~FileReader();

  void Analyze(std::vector<float>&  waveform, 
               HitCandidateVec&     hitCandidateVec, 
               const std::string&   filename, 
               const float&         bias, 
               const unsigned&      channel, 
               const Configuration& config);
  // For reco  
  void ReadFiles(SiPMToTriggerMap&           sipmToTriggerMap, 
                 const SiPMToFilesMap&       map, 
                 const unsigned&             trigger, 
                 const SiPMInfoMap&          sipmInfoMap, 
                 const wheel::Configuration& config);
  // For characterization
  void ReadFiles(SiPMToTriggerMap&           sipmToTriggerMap, 
                 const BiasToFileMap&        biasMap, 
                 const unsigned&             sipm, 
                 const wheel::Configuration& config);

 
  std::vector<TGraph>&  GetRawGraphs()  { return rawWaveforms; }
  std::vector<TGraph>&  GetModGraphs()  { return modWaveforms; }
  wheel::MarkerPairVec& GetMarkers()    { return markers; }
  
private:

  void ReadFile(HitCandidateVec& hitCandidateVec, const std::string& filename, const float& bias, const unsigned& channel, const wheel::Configuration& config);
  void MakeTheMarkers(const HitCandidateVec& hitCandVec);

  std::vector<TGraph>  rawWaveforms;
  std::vector<TGraph>  modWaveforms;
  MarkerPairVec        markers; 

};
}

#endif
