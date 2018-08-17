// 
// File: WaveformAlg.h
//
// Author: Hunter Sullivan
//
// Discription: Structure to perform various algos on the input waveforms.
//

#ifndef WAVEFORMALG_H
#define WAVEFORMALG_H

#include "Utilities.h"
#include <fstream>
#include <iostream>

namespace wheel {

class WaveformAlg {

public:
  WaveformAlg();
  ~WaveformAlg();
  
  void SmoothWaveform(std::vector<float>& signal, const wheel::Configuration& config);
  void SmoothWaveform2(std::vector<float>& signal, const wheel::Configuration& config);
  void FindHitCandidates(std::vector<float>&         waveform,
                         size_t                      roiStartTick,
                         size_t                      channel,
                         const float&                bias,
    		                 wheel::HitCandidateVec&     hitCandVec,
                         const wheel::Configuration& config); 

private:

  std::pair<float, float> ComputeNoise(std::vector<float>& signal, float& baselineInter, const wheel::Configuration& config);
  void                    FindHitCandidates(std::vector<float>::const_iterator startItr,
                                            std::vector<float>::const_iterator stopItr,
                                            const std::pair<float, float>&     noiseParameters,
                                            const float&                       bias,
                                            size_t                             roiStartTick,
		                                        HitCandidateVec&                   hitCandVec,
                                            const wheel::Configuration&        config);
};
}

#endif
