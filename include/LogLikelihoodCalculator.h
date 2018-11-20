// 
// File: LogLikelihoodCalculator.h
//
// Author: Hunter Sullivan
//
// Description: Structure to calculate log likelihood. 
//

#ifndef LOGLIKELIHOODCALCULATOR_H
#define LOGLIKELIHOODCALCULATOR_H

#include "Utilities.h"
#include "Rtypes.h"
#include <iostream>
#include "TMath.h"
#include <map>

namespace wheel {

class LogLikelihoodCalculator 
{
public:
  LogLikelihoodCalculator();
  ~LogLikelihoodCalculator();
  
  void  Initialize(const Configuration&                config, 
                   const std::map<unsigned, unsigned>& data);

  double ComputeLogLikelihood(const double& N0, 
                              const double& r, 
                              const double& thetaDeg);
  double ComputeLambda(const double&   N0, 
                       const double&   r, 
                       const double&   thetaDeg, 
                       const unsigned& m);
  double ComputeLogFactorial(const unsigned& counts);
 
  float                        m_beta;            
  float                        m_diskRadius;
  float                        m_attenuationLength;  //< input attenuation length
  unsigned                     m_nSiPMs;             //< number of sipms 
  std::map<unsigned, unsigned> m_data;               //< measured counts (sipm, np.e.)
};
}

#endif
