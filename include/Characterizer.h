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

  std::multimap<unsigned, std::vector<TH1D>>&         GetAmpDists() { return ampDists; }
  std::multimap<unsigned, std::vector<TGraphErrors>>& GetAmpPeaks() { return ampPeaks; }
    
private:

  void         MakeHistograms(const unsigned& sipm, const std::vector<HitCandidateVec>& hitCandVec, const Configuration& config);
  TGraphErrors FitGain(TH1D& hs, const unsigned& sipm, SiPMGains& sipmGains, const unsigned& nBias, const Configuration& config);

  std::multimap<unsigned, std::vector<TH1D>> ampDists;
  std::multimap<unsigned, std::vector<TGraphErrors>> ampPeaks;

};
}

#endif
