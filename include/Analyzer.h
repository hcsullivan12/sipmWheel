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
#include "TMath.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TArrayD.h"
#include "TDecompLU.h"

namespace wheel {

class Voxel {

public:
	Voxel(const float& x, const float& y, const float& r, const float& theta);
	~Voxel();
		
	float       X()     const { return m_x; };
	float       Y()     const { return m_y; };
  float       R()     const { return m_r; };
  float       Theta() const { return m_theta; };
  double      CB()    const { return m_cb; };

  void SetCB(const double& cb) { m_cb = cb; };

private:
    
  float  m_x;      ///< x position that this voxel is centered on
  float  m_y;      ///< y position that this voxel is centered on
  float  m_r;      ///< radius from center for this voxel
  float  m_theta;  ///< angle with respect to sipm 1 (in degrees)
  double m_cb;
};

class Analyzer {

public:
  Analyzer();
  ~Analyzer();
  
  void  Reconstruct(SiPMToTriggerMap& sipmToTriggerMap, const SiPMInfoMap& sipmInfoMap, const Configuration& config, const unsigned& trigger);
  void  Reconstruct(SiPMToTriggerMap& sipmToTriggerMap, const SiPMInfoMap& sipmInfoMap, const unsigned& trigger);
  float ComputeLambda(const float& r, const float& theta, const unsigned& N0, const unsigned& m);

  std::map<unsigned, unsigned> GetData()           const { return m_data; };
  
private:

  void     Initialize(const Configuration& config);
  std::pair<unsigned, unsigned> InitData(SiPMToTriggerMap& sipmToTriggerMap, const SiPMInfoMap& sipmInfoMap, const unsigned& trigger);
  double   ComputeLogLikelihood(const float& x, const float& y, const unsigned& N0);
  void     Handle(const unsigned& N0);
  void     FillAccumulatorMap();
  void     InitVoxelList();
  void     ComputeConfidenceIntervals();
  void     MakePlot();
  void     ConvertToPolar(float& r, float& thetaDeg, const float& x, const float& y);
  void     FirstEstimate(const std::pair<unsigned, unsigned>& maxSiPM_counts);
  void     RefineEstimate(unsigned& iterator);
  std::string JointCLFunc(const float& cl);
 
  double                       m_mlLogLikelihood;
  float                        m_mlN0;
  float                        m_mlX;
  float                        m_mlY;
  float                        m_mlRadius;
  float                        m_mlTheta;
  float                        m_oldGuessX;  
  float                        m_oldGuessY;  
  float                        m_oldGuessN0;
  double                       m_oldGuessMLogL;
  float                        m_newGuessX;  
  float                        m_newGuessY;
  float                        m_newGuessN0;  
  double                       m_newGuessMLogL;
  float                        m_beta;
  float                        m_diskRadius;
  float                        m_attenuationLength;
  unsigned                     m_maxIterations;
  unsigned                     m_nSiPMs;
  unsigned                     m_nVoxels;
  std::list<Voxel>             m_voxelList;
  std::map<float, float>       m_deltaX;
  std::map<float, float>       m_deltaY;
  std::map<unsigned, unsigned> m_data;  
  std::string                  m_recoOutputFile;
  std::vector<float>           m_fisherEigenvalues;
  std::vector<std::vector<float>> m_fisherEigenvectors;
  std::vector<float>           m_sigmaInverse;
  float                        m_sigmaXDiag;
  float                        m_sigmaYDiag;
  float                        m_sigmaX;
  float                        m_sigmaY;
  float                        m_correlation;
  bool                         m_doCI;
};
}

#endif
