// 
// File: Analyzer.h
//
// Author: Hunter Sullivan
//
// Description: Structure to perfom reconstruction on sipm data.
//

#ifndef ANALYZER_H
#define ANALYZER_H

#include "Utilities.h"
#include "WaveformAlg.h"
#include "Rtypes.h"
#include "TThread.h"
#include <iostream>

namespace wheel {

class Analyzer {

public:
  Analyzer();
  ~Analyzer();
  
  void  RunReco(SiPMToTriggerMap& sipmToTriggerMap, const SiPMInfoMap& sipmInfoMap, const Configuration& config, const unsigned& trigger);
  void  WheelReco(SiPMToTriggerMap& sipmToTriggerMap, const SiPMInfoMap& sipmInfoMap, const unsigned& trigger);
  float ComputeLambda(const float& r, const float& theta, const unsigned& N0, const unsigned& m, const float& attenuationLength);

  AccumulatorMap&            GetAccumulatorMap() { return p_accumulatorMap; }
  AccumulatorMap&            GetMLEParameters()  { return p_mleParameters; }
  std::map<unsigned, float>& GetData()           { return p_data; }
  
private:

  unsigned CountPhotons(SiPMToTriggerMap& sipmToTriggerMap, const SiPMInfoMap& sipmInfoMap, const unsigned& trigger);
  float    ComputeLikelihood(const float& r, const float& theta, const unsigned& N0, const float& attenuationLength);
  void     Handle(const unsigned& N0);
  void     FillAccumulatorMap();

  float                     p_thetaBinSize;
  float                     p_radiusBinSize;
  float                     p_attenuationLengthBinSize;
  float                     p_likelihood;
  float                     p_N0;
  float                     p_attenuationLength;
  float                     p_radius;
  float                     p_theta;
  float                     p_beta;
  float                     p_diskRadius;
  unsigned                  p_nSiPMs;
  std::map<unsigned, float> p_data;  
  AccumulatorMap            p_accumulatorMap;
  AccumulatorMap            p_mleParameters;
};
}

#endif
