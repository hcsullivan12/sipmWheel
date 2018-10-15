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

  AccumulatorMap&            GetAccumulatorMap() { return m_accumulatorMap; }
  AccumulatorMap&            GetMLEParameters()  { return m_mleParameters; }
  std::map<unsigned, float>& GetData()           { return m_data; }
  
private:

  std::pair<unsigned, unsigned> InitData(SiPMToTriggerMap& sipmToTriggerMap, const SiPMInfoMap& sipmInfoMap, const unsigned& trigger);
  float    ComputeLikelihood(const float& r, const float& theta, const unsigned& N0, const float& attenuationLength);
  void     Handle(const unsigned& N0);
  void     FillAccumulatorMap();

  float                     m_thetaBinSize;
  float                     m_radiusBinSize;
  float                     m_attenuationLengthBinSize;
  float                     m_likelihood;
  float                     m_N0;
  float                     m_attenuationLength;
  float                     m_radius;
  float                     m_theta;
  float                     m_beta;
  float                     m_diskRadius;
  unsigned                  m_nSiPMs;
  std::map<unsigned, float> m_data;  
  AccumulatorMap            m_accumulatorMap;
  AccumulatorMap            m_mleParameters;
};
}

#endif
