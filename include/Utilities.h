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
  float  startTickAmp;
  float  stopTickAmp;
  float  bias;
  size_t nPhotons;
  float  hitBase;
  float  hitPeak;
  float  hitPeakTick;
  float  hitAmplitude;
  float  hitIntegral;
};

struct SiPMInfo 
{
  float    gain;
  float    breakdown;
  float    bias;
};


using HitCandidateVec      = std::vector<HitCandidate>;
using MergeHitCandidateVec = std::vector<HitCandidateVec>;
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
  std::string      pathToData;                         ///< Path to csv data files
  std::string      rawWaveformPath;                    ///< Output path for raw waveforms
  std::string      modWaveformPath;                    ///< Output path for modified waveforms
  std::string      pathToConfig;                       ///< Path to configuration file
  std::string      recoOutputFile;                     ///< Output path for reconstruction results
  bool             printFiles;                         ///< Option to print out data files
  bool             baselineSubtract;                   ///< Option to baseline subtract
  bool             saveRawWaveforms;                   ///< Option to save raw waveforms
  bool             saveModWaveforms;                   ///< Option to save modified waveforms
  unsigned         nSiPMs;                             ///< Number of sipms 
  unsigned         smaRange;                           ///< Range around center for applying average on waveform
  unsigned         resolution;                         ///< Resolution to keep on waveform
  float            hitSigma;                           ///< 
  unsigned         hitFinderSearch;                    
  float            minimumHitAmp;                      ///< Hit threshold
  std::string      process;                            ///< Process: Reconstruction or characterization
  float            characterizeAmpThr;
  float            characterizeAmpSig;
  float            characterizeAmpFitRange;
  std::string      characterizeOutputFile;
  unsigned         nFilesCharacterize;                 ///< Number of files to use from data
  unsigned                  nBiases;                   ///< Number of biases used
  std::set<float>           biases;                    ///< The bias values (in V)
  std::map<unsigned, float> gains;                     ///< SiPM gains (in mV/p.e./O.V.)
  std::map<unsigned, float> breakdowns;                ///< SiPM breakdowns (in V)
  float                     thetaBinSize;              
  float                     radiusBinSize;
  float                     attenuationLengthBinSize;
  float                     diskRadius;                ///< Disk radius for reconstruction (in cm)
  float                     attenuationLength;         ///< Attenuation length for disk (in cm)
  unsigned                  nVoxels;                   ///< Number of voxels to segment the disk geometry (=NxN)
  unsigned                  maxIterations;
};
}

#endif
