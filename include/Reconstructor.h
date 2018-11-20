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
#include "LogLikelihoodCalculator.h"
#include "LogLikelihoodFunction.h"
#include <iostream>
#include "TMath.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

namespace wheel {

class LogLikelihoodCalculator;
class LogLikelihoodFunc; 

class Reconstructor {

public:
  Reconstructor();
  ~Reconstructor();
  
  void  Reconstruct(SiPMToTriggerMap&    sipmToTriggerMap, 
                    const SiPMInfoMap&   sipmInfoMap, 
                    const Configuration& config, 
                    const unsigned&      trigger);
  void  Reconstruct(unsigned& N0);
  void  SetData(const Configuration& config, 
                const std::map<unsigned, unsigned>& data) { m_data = data; m_logLikeFunc.Initialize(config, m_data); }
  void  Initialize(const Configuration& config);
  void  MakePlot(const unsigned& trigger);

  double ML() { return m_mlLogLikelihood; };
  double X() { return m_mlX; };
  double Y() { return m_mlY; };
  double Radius() { return m_mlRadius; };
  double Theta() { return m_mlTheta; };
  double N0() { return m_mlN0; };
    
private:

  void     InitData(SiPMToTriggerMap&  sipmToTriggerMap, 
                     const SiPMInfoMap& sipmInfoMap, 
                     const unsigned&    trigger);
  void     ConvertToPolar(float&       r, 
                          float&       thetaDeg, 
                          const float& x, 
                          const float& y);
 
  double                       m_mlLogLikelihood; //< Log likelihood for the MLE
  double                       m_mlN0;            //< MLE for N0
  double                       m_mlX;             //< MLE for x (cm)
  double                       m_mlY;             //< MLE for y (cm)
  double                       m_mlRadius;        //< MLE for r (cm)
  double                       m_mlTheta;         //< MLE for theta (deg)

  float                        m_beta;            
  float                        m_diskRadius;
  float                        m_attenuationLength;  //< input attenuation length
  unsigned                     m_maxIterations;      //< maximum number of iterations minimizer
  unsigned                     m_nSiPMs;             //< number of sipms 
  std::map<unsigned, unsigned> m_data;               //< measured counts (sipm, np.e.)
  unsigned                     m_totalCounts;        //< total amount of counts
  unsigned                     m_maxCounts;          //< maximum number of p.e., used for plotting
  unsigned                     m_maxSiPM;            //< sipm ID which saw maxCounts
  std::string                  m_recoOutputPath;     //< output file for plots
  LogLikelihoodFunction        m_logLikeFunc; 
  LogLikelihoodCalculator      m_logLikeCalc;
};
}

#endif
