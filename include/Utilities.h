// 
// File: Utilities.h
//
// Author: Hunter Sullivan
//
// Detail: Structure to hold configuration settings and defining some useful
//         containers.
//

#ifndef UTILITIES_H
#define UTILITIES_H

#include <vector>
#include <map>
#include <list>
#include <set>
#include "TMarker.h"

namespace wheel {

struct HitCandidate
{
  size_t channel;
  size_t startTick;
  size_t stopTick;
  float  bias;
  size_t nPhotons;
  float  hitBase;
  float  hitCenter;
  float  hitHeight;
};

struct SiPMInfo 
{
  float    gain;
  float    breakdown;
  float    bias;
};


using HitCandidateVec      = std::vector<HitCandidate>;
using SiPMToHitCandVecMap  = std::map<unsigned, HitCandidateVec>;
using SiPMToFilesMap       = std::map<unsigned, std::set<std::string>>;
using SiPMToTriggerMap     = std::map<unsigned, std::vector<HitCandidateVec>>;
using SiPMToBiasTriggerMap = std::map<unsigned, std::map<float, std::set<std::string>>>;  
using SiPMGains            = std::multimap<unsigned, float>;
using SiPMInfoMap          = std::map<unsigned, SiPMInfo>;
using BinIndex             = std::pair<unsigned, std::map<std::string, float>>;
using AccumulatorMap       = std::vector<BinIndex>;  // N0, r, theta, likelihood
using DataList             = std::list<std::map<unsigned, float>>;
using AccumulatorMapList   = std::list<AccumulatorMap>;
using BiasToFileMap        = std::map<float, std::set<std::string>>;
using MarkerPairVec        = std::vector<std::vector<std::pair<TMarker,TMarker>>>;

struct Configuration 
{
  std::string      pathToData;
  std::string      rawWaveformPath;
  std::string      modWaveformPath;
  std::string      pathToConfig;
  std::string      recoOutputFile;
  bool             printFiles;
  bool             baselineSubtract;
  bool             saveRawWaveforms;
  bool             saveModWaveforms;
  unsigned         nSiPMs;
  unsigned         smaRange;
  unsigned         resolution;
  float            hitSigma;
  unsigned         hitFinderSearch;
  float            minimumHitAmp;
  std::string      process;
  float            characterizeAmpThr;
  float            characterizeAmpSig;
  float            characterizeAmpFitRange;
  std::string      characterizeOutputFile;
  unsigned         nFilesCharacterize;
  unsigned                  nBiases;
  std::set<float>           biases;
  std::map<unsigned, float> gains;
  std::map<unsigned, float> breakdowns;
  float                     thetaBinSize;
  float                     radiusBinSize;
  float                     attenuationLengthBinSize;
  float                     diskRadius;
};
}

#endif
