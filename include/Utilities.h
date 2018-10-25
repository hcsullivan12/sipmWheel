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
using SiPMGains            = std::vector<std::vector<std::vector<float>>>;
using SiPMInfoMap          = std::map<unsigned, SiPMInfo>;
using BinIndex             = std::pair<unsigned, std::map<std::string, float>>;
using AccumulatorMap       = std::vector<BinIndex>;  // N0, r, theta, likelihood
using DataList             = std::list<std::map<unsigned, float>>;
using AccumulatorMapList   = std::list<AccumulatorMap>;
using BiasToFileMap        = std::map<float, std::set<std::string>>;
using MarkerListVec        = std::vector<std::list<TMarker>>;

struct Configuration 
{
  // Name of the process
  std::string               process;                   ///< Process: Reconstruction or characterization
  // File handling information
  std::string               pathToData;                ///< Path to csv data files
  std::string               rawWaveformPath;           ///< Output path for raw waveforms
  std::string               modWaveformPath;           ///< Output path for modified waveforms
  std::string               pathToConfig;              ///< Path to configuration file
  std::string               recoOutputPath;            ///< Output path for reconstruction results
  std::string               characterizeOutputPath;    ///< Output path for characterization results
  std::string               simulateOutputPath;        ///< Output path for simulation results
  // Output options and signal processing
  bool                      printFiles;                ///< Option to print out data files
  bool                      baselineSubtract;          ///< Option to baseline subtract
  bool                      saveRawWaveforms;          ///< Option to save raw waveforms
  bool                      saveModWaveforms;          ///< Option to save modified waveforms
  bool                      smoothWaveform;            ///< Option to perform a running average
  unsigned                  smaRange;                  ///< Range around center for applying average on waveform
  unsigned                  resolution;                ///< Resolution to keep on waveform
  unsigned                  hitFinderSearch;           ///< Number of samples to skip in hit finding
  float                     minimumHitAmp;             ///< Hit threshold
  // Reconstruction parameters
  unsigned                  nVoxels;                   ///< Number of voxels to segment the disk geometry (=NxN)
  unsigned                  maxIterations;             ///< Max number of iterations for NR correction
  // Characterization parameters
  float                     characterizeAmpThr;        ///< Threshold for peak fitting
  float                     characterizeAmpSig;        ///< Spread for peak finding
  float                     characterizeAmpFitRange;   ///< Spread for fitting
  unsigned                  nFilesCharacterize;        ///< Number of files to use from data
  // Disk and sipm information
  unsigned                  nBiases;                   ///< Number of biases used
  std::set<float>           biases;                    ///< The bias values (in V)
  std::map<unsigned, float> gains;                     ///< SiPM gains (in mV/p.e./O.V.)
  std::map<unsigned, float> breakdowns;                ///< SiPM breakdowns (in V)
  float                     diskRadius;                ///< Disk radius (in cm)
  float                     diskThickness;             ///< Thickness of the disk (in cm)
  unsigned                  nSiPMs;                    ///< Number of sipms 
  float                     sipmArea;                  ///< Area of sipms
  float                     attenuationLength;         ///< Attenuation length for disk (in cm)
  // Simulation information
  std::vector<float>        sourcePosition;            ///< Light source position. radius from center and theta measured counterclockwise from sipm1 (cm, deg)
  float                     sourceSigma;               ///< Spread in light position
  float                     smearSigma;                ///< Spread in incidence angle
  float                     tpbEmissionPeak;           ///< TPB emission peak wavelength (nm)
  float                     indexRefractionDisk;       ///< Index of refraction for the disk
  float                     indexRefractionEnv;        ///< Index of refraction for the sorrounding environment
  unsigned                  nPhotonsToLaunch;          ///< Number of photons to simulate
  float                     terminationThreshold;      ///< Threshold for photon termination
  float                     bulkAttenuation;           ///< Bulk attenuation length
  float                     bulkAbsorption;            ///< Bulk absorption length
  float                     surfaceAbsorptionCoeff;    ///< Surface absorption coeffecient
};
}

#endif
