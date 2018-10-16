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

class Voxel {

public:
	Voxel(const float& x, const float& y, const float& r, const float& theta);
	~Voxel();
		
	float X()     const { return m_x; };
	float Y()     const { return m_y; };
  float R()     const { return m_r; };
  float Theta() const { return m_theta; }

private:
    
  float m_x;      ///< x position that this voxel is centered on
  float m_y;      ///< y position that this voxel is centered on
  float m_r;      ///< radius from center for this voxel
  float m_theta;  ///< angle with respect to sipm 1 (in degrees)
};

class Analyzer {

public:
  Analyzer();
  ~Analyzer();
  
  void  Reconstruct(SiPMToTriggerMap& sipmToTriggerMap, const SiPMInfoMap& sipmInfoMap, const Configuration& config, const unsigned& trigger);
  void  Reconstruct(SiPMToTriggerMap& sipmToTriggerMap, const SiPMInfoMap& sipmInfoMap, const unsigned& trigger);
  float ComputeLambda(const float& r, const float& theta, const unsigned& N0, const unsigned& m);

  AccumulatorMap            GetAccumulatorMap() const { return m_accumulatorMap; };
  AccumulatorMap            GetMLEParameters()  const { return m_mleParameters; };
  std::map<unsigned, float> GetData()           const { return m_data; };
  
private:

  void     Initialize(const Configuration& config);
  std::pair<unsigned, unsigned> InitData(SiPMToTriggerMap& sipmToTriggerMap, const SiPMInfoMap& sipmInfoMap, const unsigned& trigger);
  float    ComputeLikelihood(const float& r, const float& theta, const unsigned& N0);
  void     Handle(const unsigned& N0);
  void     FillAccumulatorMap();
  void     InitVoxelList();
  void     ComputeConfidenceIntervals();

  float                     m_mlLikelihood;
  float                     m_mlN0;
  float                     m_mlX;
  float                     m_mlY;
  float                     m_mlRadius;
  float                     m_mlTheta;
  float                     m_beta;
  float                     m_diskRadius;
  float                     m_attenuationLength;
  unsigned                  m_nSiPMs;
  unsigned                  m_nVoxels;
  std::list<Voxel>          m_voxelList;
  std::map<unsigned, float> m_data;  
  AccumulatorMap            m_accumulatorMap;
  AccumulatorMap            m_mleParameters;
};
}

#endif
