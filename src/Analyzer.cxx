//
// File: Analyzer.cxx
//
// Author: Hunter Sullivan
//
// Description: Structure to run reconstruction on sipm data.
//

#include "Analyzer.h"
#include "TMath.h"

namespace wheel {

Analyzer::Analyzer()
{}

Analyzer::~Analyzer()
{}

void Analyzer::Reconstruct(SiPMToTriggerMap& sipmToTriggerMap, const SiPMInfoMap& sipmInfoMap, const Configuration& config, const unsigned& trigger)
{
  // Initialize parameters
  Init(config);
  // Start reconstruction
  Reconstruct(sipmToTriggerMap, sipmInfoMap, trigger);
}

void Analyzer::Init(const Configuration& config)
{
  // Initialize
  m_thetaBinSize      = config.thetaBinSize;
  m_radiusBinSize     = config.radiusBinSize;
  m_diskRadius        = config.diskRadius;
  //m_attenuationLengthBinSize = config.attenuationLengthBinSize;
  m_attenuationLength = config.attenuationLength;
  m_beta              = 360/config.nSiPMs;
  m_nSiPMs            = config.nSiPMs;
  m_mlLikelihood      = std::numeric_limits<float>::lowest();
  m_mlRadius = 0; m_mlTheta = 0; m_mlN0 = 0; 
  m_data.clear();
}

void Analyzer::Reconstruct(SiPMToTriggerMap& sipmToTriggerMap, const SiPMInfoMap& sipmInfoMap, const unsigned& trigger)
{
  // Start the timer for this trigger
  clock_t start = clock();

  // First to cound number of photons in each hit
  // also N0 set to lower threshold to decrease computation time
  const auto maxSiPM_counts = InitData(sipmToTriggerMap, sipmInfoMap, trigger);
  unsigned N0 = maxSiPM_counts.second;

  std::cout << "Running MLE...\n";
  // Start main loop
  while ( N0 <= 300 ) 
  {
    Handle(N0);
    N0++;
  }
 
  // We should have the ml now
  // Loop again and fill the accumulator map
  // for mle N0 
  FillAccumulatorMap(); 
  // Maximum is at the end
  unsigned ID = m_accumulatorMap.size()-1;
  std::cout << "Max likelihood = " << m_accumulatorMap[ID].second.find("likelihood")->second        << std::endl
            << "Radius         = " << m_accumulatorMap[ID].second.find("radius")->second            << " cm\n"
            << "Theta          = " << m_accumulatorMap[ID].second.find("theta")->second             << " deg\n"
            << "N0             = " << m_accumulatorMap[ID].first                                    << " photons\n";
 
  clock_t end = clock();
  double duration = ((double) (end - start)) / CLOCKS_PER_SEC;
  std::cout << "Run time of " << duration << " s" << std::endl;

  // Add these parameters to our mle vec
  m_mleParameters.emplace_back(m_accumulatorMap[ID].first, m_accumulatorMap[ID].second);
}

void Analyzer::Handle(const unsigned& N0)
{
  // Variables used for incrementing
  // Start with r = 0, theta = 0, (N0 is what's passed) and increment by the constants defined above
  float radius(0);
  float theta(0);

  //std::cout << "N0: " << N0 << std::endl;
  while ( radius <= m_diskRadius ) 
  {
    //std::cout << "radius: " << radius << std::endl;
    while ( theta < 360 ) 
    {
      // Since there are so many iterations here,
      // let's only fill the accumulator map 
      // with the max likelihood iterations
      // So find the N0 and attenL that correspond to ML
      float likelihood = ComputeLikelihood(radius, theta, N0);
      //std::cout << likelihood << std::endl;
      if ( likelihood > m_mlLikelihood) { m_mlLikelihood = likelihood; m_mlN0 = N0; m_mlRadius = radius; m_mlTheta = theta; } 
        
      theta = theta + m_thetaBinSize;
    }
    theta = 0;
    radius = radius + m_radiusBinSize;
  }
  radius = 0;
  //std::cout << "ML " << m_likelihood << "  atten: " << m_attenuationLength << "  radius " << m_radius << "  theta " << m_theta << "  N0 " << m_N0 << std::endl;
}

void Analyzer::FillAccumulatorMap()
{
  // Variables used for incrementing
  // Start with r = 0, theta = 0, (N0 is what's passed) and increment by the constants defined above
  float radius = 0;
  float theta  = 0;

  while ( radius <= m_diskRadius ) 
  {
    while ( theta < 360 ) 
    {
      // Compute likelihood using the N0 and attenL for ml already found
      float likelihood = ComputeLikelihood(radius, theta, m_mlN0);

      std::map<std::string, float> params;
      params.emplace("radius", radius);         
      params.emplace("theta", theta); 
      params.emplace("likelihood", likelihood); 
      m_accumulatorMap.emplace_back(m_mlN0, params);

      theta = theta + m_thetaBinSize;
    }
    theta = 0;
    radius = radius + m_radiusBinSize;
  }

  // Now sort for output purposes
  std::sort(m_accumulatorMap.begin(), m_accumulatorMap.end(), 
            [](const BinIndex& biL, const BinIndex& biR) 
            {return biL.second.find("likelihood")->second < biR.second.find("likelihood")->second;});
}

std::pair<unsigned, unsigned> Analyzer::InitData(SiPMToTriggerMap& sipmToTriggerMap, const SiPMInfoMap& sipmInfoMap, const unsigned& trigger)
{ 
  // Which sipm saw the largest number of photons?
  unsigned maxSiPM;
  unsigned max(0);

  // Loop over all the hits and set the n of photons
  for (auto& sipm : sipmToTriggerMap)
  {
    // Get the gain and bdown for this sipm
    const float& thisGain = sipmInfoMap.find(sipm.first)->second.gain;
    const float& thisBD   = sipmInfoMap.find(sipm.first)->second.breakdown;
    // Total number for this sipm
    unsigned sipmCounts(0);
    for (auto& hit : sipm.second[0])
    {
      //std::cout << hit.hitHeight << "  " << hit.bias << "  " << thisBD << "\n";
      hit.nPhotons = std::round(hit.hitAmplitude/( thisGain*( hit.bias - thisBD ) ));   // hit.bias - thisBD
      sipmCounts += hit.nPhotons;
      //std::cout << hit.nPhotons << std::endl;
    }
    if (sipmCounts > max) { max = sipmCounts; maxSiPM = sipm.first; }

    // Store counts in data container
    m_data.emplace(sipm.first, sipmCounts);
    std::cout << "SiPM " << sipm.first << " --> " << sipmCounts << " p.e.\n";
  }

  return std::make_pair(maxSiPM, max);
}

float Analyzer::ComputeLikelihood(const float& r, const float& theta, const unsigned& N0)
{
  // Sum over terms ----> k_m*ln(lambda_m) - lambda_m - ln(k_m!)
  float sum = 0;
  for (int sipm = 1; sipm <= m_nSiPMs; sipm++) {
    float lambda_m = ComputeLambda(r, theta, N0, sipm);
    //std::cout << "nPhotons: " << m_data.find(m)->second << "  lambda_m " << lambda_m << "   factorial " << TMath::Factorial(m_data.find(m)->second) << "  term ";
    float term = m_data.find(sipm)->second*log(lambda_m) - lambda_m - log(TMath::Factorial(m_data.find(sipm)->second));
    //std::cout << term << std::endl;
    sum = sum + term;
  }
  return sum;
}

float Analyzer::ComputeLambda(const float& r, const float& theta, const unsigned& N0, const unsigned& sipm)
{
  // Assumptions:
  //    1) Only bulk absorption 
  //    2) Detection effeciency of sipms is 100%
  //
  // The equation I'll be using is...
  // 
  // I(r,theta) = (I_0/r_m)*cos(theta)*exp(-r_m/attenuationLength)
  //
  //    1) 1/r_m term comes from this being a 2D wave from point-like source
  //    2) cos(theta) comes from the projection onto the sipm normal axis
  //    3) Exponential term comes from bulk absorption
  
  // Compute angle between sipm normal and v = (x,y), and distance from (x,y) to sipm
  //float r                 = std::sqrt(x*x + y*y);
  //float thetaDeg          = std::abs(TMath::ASin(y/r))*180/TMath::Pi();
  // Integral 
  //float b                 = 2*m_diskRadius/m_attenuationLength;
  //float pi                = TMath::Pi();
  //float integralNorm      = m_attenuationLength*(2 - 0.467289);
  // Handle theta convention
  // 2nd quadrant
  //if (x < 0 && y > 0) thetaDeg += 90;
  // 3rd quadrant
  //if (x < 0 && y < 0) thetaDeg += 180;
  // 4th quadrant
  //if (x > 0 && y < 0) thetaDeg += 270;
  // Remember, we overlaid a square grid on our disk, so some of the voxels
  // we need to ignore
  if (r >= m_diskRadius) return 0; 
  float angleXYandSiPMRad = ( (sipm - 1)*m_beta - theta )*( TMath::Pi()/180 );
  float sipmToXYSquared   = r*r + m_diskRadius*m_diskRadius - 2*r*m_diskRadius*TMath::Cos(angleXYandSiPMRad);
    
  // Safety here
  if (sipmToXYSquared < 0) 
  {
    std::cout << "Error! Square root of negative number!\n";
    return 0;
  }
  float sipmToXY = TMath::Sqrt(sipmToXYSquared);
  // Get the angle between sipmToXY and the normal for this sipm
  // Put a protection here
  float cosAngleSiPMToXYandSiPM(1);
  if (sipmToXY != 0) cosAngleSiPMToXYandSiPM = (-r*r + sipmToXYSquared + m_diskRadius*m_diskRadius)/(2*sipmToXY*m_diskRadius);
 
  float relativeWeight = N0*cosAngleSiPMToXYandSiPM*TMath::Exp(-sipmToXY/m_attenuationLength)/sipmToXY;
  //std::cout << "Weight = " << relativeWeight << "  at sipm " << sipm << " from x = " << x << " y = " << y << "  sipmToXY = " << sipmToXY << " angleXYandSIPM = " << angleXYandSiPMRad*180/TMath::Pi() << " beta = " << m_beta << "  thetaDeg = " << thetaDeg <<std::endl;
  if (relativeWeight < 0) { std::cout << "UH OH! WEIGHT < 0!!\n"; std::exit(1); }
  return relativeWeight;
}
}
