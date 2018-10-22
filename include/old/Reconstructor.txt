// 
// File: Reconstructor.h
//
// Author: Hunter Sullivan
//
// Description: Reconstruction algorithm for sipm data.
//              Our disk will be segmented into voxels, where the 
//              corresponding weights, c_ij, represent the probability
//              that a photon leaving pixel j is detected by sipm i.
//              The likelihood is maximized using an EM algorithm.

#ifndef Reconstructor_H
#define Reconstructor_H

#include "Utilities.h"
#include "WaveformAlg.h"
#include "Rtypes.h"
#include "TThread.h"
#include <iostream>

namespace wheel {
	
class Voxel {

public:
	Voxel(const float& x, const float& y);
	~Voxel();
	
	void SetOldIntensity(const float& intensity) { m_oldIntensity = intensity; };
	void SetNewIntensity(const float& intensity) { m_newIntensity = intensity; };
  void SetProbabilityMap(const std::map<unsigned, float>& pMap) { m_probabilityMap = pMap; }; 
	
	float X()              const { return m_x; };
	float Y()              const { return m_y; };
	float OldIntensity()   const { return m_oldIntensity; };
	float NewIntensity()   const { return m_newIntensity; };
	std::map<unsigned, float> ProbabilityMap() const { return m_probabilityMap; };

private:
    
  float m_x;                                   ///< x position that this voxel is centered on
  float m_y;                                   ///< y position that this voxel is centered on
  float m_oldIntensity;                        ///< old intensity estimate for this voxel
  float m_newIntensity;                        ///< new intensity estimate for this voxel
  std::map<unsigned, float> m_probabilityMap;  ///< probability that photon leaving this voxel reaches position "first"
};

class Reconstructor {

public:
  Reconstructor();
  ~Reconstructor();
  
  void Reconstruct(SiPMToTriggerMap& sipmToTriggerMap, const SiPMInfoMap& sipmInfoMap, const Configuration& config, const unsigned& trigger);
  
private:

  void Initialize(const Configuration& config);
  void InitVoxelList();
  void Reconstruct(SiPMToTriggerMap& sipmToTriggerMap, const SiPMInfoMap& sipmInfoMap, const unsigned& trigger);
  void FirstEstimate(const std::pair<unsigned, unsigned>& maxSiPM_counts);
  std::pair<unsigned, unsigned> InitData(SiPMToTriggerMap& sipmToTriggerMap, const SiPMInfoMap& sipmInfoMap, const unsigned& trigger);
  float ComputeWeight(const unsigned& sipm, const float& x, const float& y);
  void CheckConvergence(float& eps);
  void MakeNextEstimate(unsigned& iterator);
  void SaveResults();


    
  std::list<Voxel>          m_voxelList;         ///< Container for voxels
  float                     m_attenuationLength; ///< Attenuation length for the disk
  float                     m_beta;              ///< Angle between sipm positions
  float                     m_diskRadius;        ///< Disk radius
  unsigned                  m_nSiPMs;            ///< Number of sipms around the disk
  unsigned                  m_nVoxels;           ///< Granularity of source position
  std::string               m_recoOutputFile;    ///< Output path for results
  std::map<unsigned, float> m_data;              ///< Number of photons detected by each sipm 
};
}

#endif
