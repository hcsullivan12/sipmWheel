// 
// File: Reconstructor.h
//
// Author: Hunter Sullivan
//
// Description: Structure to perfom reconstruction on sipm data.
//

#ifndef RECONSTRUCTOR_H
#define RECONSTRUCTOR_H

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

class Reconstructor {

public:
  Reconstructor();
  ~Reconstructor();
  
  void  Reconstruct(SiPMToTriggerMap& sipmToTriggerMap, const SiPMInfoMap& sipmInfoMap, const Configuration& config, const unsigned& trigger);
  void  Reconstruct(unsigned& N0);
  void  SetData(const std::map<unsigned, unsigned>& data) { m_data = data; };
  void  Initialize(const Configuration& config);
  void  MakePlot(const unsigned& trigger);


  const double   ML()    { return m_mlLogLikelihood; }
  const float    X()     { return m_mlX; }
  const float    Y()     { return m_mlY; }
  const float    R()     { return m_mlRadius; }
  const float    Theta() { return m_mlTheta; }
  const unsigned N0()    { return m_mlN0; }
    
private:

  std::pair<unsigned, unsigned> InitData(SiPMToTriggerMap& sipmToTriggerMap, const SiPMInfoMap& sipmInfoMap, const unsigned& trigger);
  double   ComputeLogLikelihood(const float& x, const float& y, const unsigned& N0);
  float    ComputeLambda(const float& r, const float& theta, const unsigned& N0, const unsigned& m);
  void     Handle(const unsigned& N0);
  void     InitVoxelList();
  void     ConvertToPolar(float& r, float& thetaDeg, const float& x, const float& y);
  void     RefineEstimate(unsigned& iterator);
 
  double                       m_mlLogLikelihood; //< Log likelihood for the MLE
  float                        m_mlN0;            //< MLE for N0
  float                        m_mlX;             //< MLE for x (cm)
  float                        m_mlY;             //< MLE for y (cm)
  float                        m_mlRadius;        //< MLE for r (cm)
  float                        m_mlTheta;         //< MLE for theta (deg)

  float                        m_oldGuessX;       //< old estimate for x in NR method
  float                        m_oldGuessY;       //< old estimate for y in NR method
  float                        m_oldGuessN0;      //< old estimate for N0 in NR method
  double                       m_oldGuessMLogL;   //< old estimate for log likelihood in NR method

  float                        m_beta;            
  float                        m_diskRadius;
  float                        m_attenuationLength;  //< input attenuation length
  unsigned                     m_maxIterations;      //< maximum number of iterations for NR method
  unsigned                     m_nSiPMs;             //< number of sipms 
  unsigned                     m_nVoxels;            //< number of voxels used to segment disk
  std::list<Voxel>             m_voxelList;          //< list of created voxels
  std::map<unsigned, unsigned> m_data;               //< measured counts (sipm, np.e.)
  unsigned                     m_maxCounts;          //< maximum number of p.e., used for plotting
  std::string                  m_recoOutputPath;     //< output file for plots
};
}

#endif
