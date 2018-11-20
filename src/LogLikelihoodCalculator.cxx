//
// File: LogLikelihoodCalculator.cxx
//
// Author: Hunter Sullivan
//
// Description: Structure to calculate log likelihood.
//

#include "LogLikelihoodCalculator.h"
#include "TMath.h"

namespace wheel {

LogLikelihoodCalculator::LogLikelihoodCalculator()
{}

LogLikelihoodCalculator::~LogLikelihoodCalculator()
{}

void LogLikelihoodCalculator::Initialize(const Configuration&                config, 
                                         const std::map<unsigned, unsigned>& data)
{
  // Initialize
  m_diskRadius        = config.diskRadius;
  m_attenuationLength = config.attenuationLength;
  m_beta              = 360/config.nSiPMs;
  m_nSiPMs            = config.nSiPMs;
  m_data              = data;

  if (m_diskRadius <= 0) { std::cout << "Error. Please specify disk radius!\n"; std::exit(1); }
  if (m_attenuationLength <= 0) { std::cout << "Error. Please specify an attenuation length!\n"; std::exit(1); }
  if (config.nSiPMs <= 0) { std::cout << "Error. Please specify number of SiPMs!\n"; std::exit(1); }
}

double LogLikelihoodCalculator::ComputeLogLikelihood(const double& N0, 
                                                     const double& r, 
                                                     const double& thetaDeg)
{
  // Sum over terms ----> k_m*ln(lambda_m) - lambda_m - ln(k_m!)
  double sum = 0;
  for (int sipm = 1; sipm <= m_nSiPMs; sipm++) 
  {
    double   lambda_m = ComputeLambda(N0, r, thetaDeg, sipm);
    unsigned k_m      = m_data.find(sipm)->second;
    // Since we may be dealing with large numbers, we need to handle the factorial carefully
    double term = k_m*std::log(lambda_m) - lambda_m - ComputeLogFactorial(k_m);
    sum = sum + term;

    //std::cout << lambda_m << " " << k_m << " " << term << std::endl;
  }
  return sum;
} 

double LogLikelihoodCalculator::ComputeLogFactorial(const unsigned& counts)
{
  double sum(0);
  for (unsigned n = 1; n <= counts; n++) sum = sum + std::log(n);
  return sum;
}

double LogLikelihoodCalculator::ComputeLambda(const double&   N0, 
                                              const double&   r, 
                                              const double&   thetaDeg, 
                                              const unsigned& sipm)
{ 
  // Assumptions:
  //    1) Only bulk absorption 
  //    2) Detection effeciency of sipms is 100%
  //
  // The equation I'll be using is...
  // 
  // I(r,theta) = (I_0/r_m)*cos(alpha)*exp(-r_m/attenuationLength)
  //
  //    1) 1/r_m term comes from this being a 2D wave from point-like source
  //    2) cos(alpha) comes from the projection onto the sipm normal axis
  //    3) Exponential term comes from bulk absorption
  
  // Compute angle between sipm normal and v = (x,y), and distance from (x,y) to sipm
  double radToDeg     = 180/TMath::Pi();
    
  double angleXYandSiPMRad = ((sipm - 1)*m_beta - thetaDeg)*(1/radToDeg);
  double sipmToXYSquared   = r*r + m_diskRadius*m_diskRadius - 2*r*m_diskRadius*TMath::Cos(angleXYandSiPMRad);
    
  // Safety here
  if (sipmToXYSquared < 0) { std::cout << "Error! Square root of negative number!\n"; return 0; } 
  double sipmToXY = std::sqrt(sipmToXYSquared);

  // Get the angle between sipmToXY and the normal for this sipm
  // Put a protection here
  double cosAngleSiPMToXYandSiPM(1);
  if (sipmToXY != 0) cosAngleSiPMToXYandSiPM = (-r*r + sipmToXYSquared + m_diskRadius*m_diskRadius)/(2*sipmToXY*m_diskRadius);
 
  double weight = N0*cosAngleSiPMToXYandSiPM*TMath::Exp(-sipmToXY*sipmToXY/(m_attenuationLength*m_attenuationLength));
  return weight;
}
}

