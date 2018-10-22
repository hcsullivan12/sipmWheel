// 
// File: Characterizer.h
//
// Author: Hunter Sullivan
//
// Discription: Structure to characterize sipms.
//

#ifndef CHARACTERIZER_H
#define CHARACTERIZER_H

#include "Utilities.h"
#include <iostream>
#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TStyle.h"
#include "TGraphErrors.h"

namespace wheel {

class Characterizer {

public:
  Characterizer();
  ~Characterizer();
  
  void Characterize(SiPMInfoMap& sipmInfoMap, const SiPMToTriggerMap& sipmToTriggerMap, const Configuration& config);
  void SaveCharacterizationPlots(const wheel::Configuration& config);
  void Initialize(const Configuration& config);

private:

  void         MakeHistograms(const unsigned& sipm, const std::vector<HitCandidateVec>& hitCandVec, const Configuration& config);
  TGraphErrors FitGain(TH1D& hs, const unsigned& sipm, const unsigned& nBias, const Configuration& config);

  std::map<unsigned, std::vector<TH1D>>         m_ampDists;
  std::map<unsigned, std::vector<TGraphErrors*>> m_ampPeaks;
  SiPMGains m_sipmGains;

};
}

#endif
