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

void Analyzer::Handle(const unsigned& N0)
{
  //long N0 = (long) ptr;

  // Variables used for incrementing
  // Start with r = 0, theta = 0, (N0 is what's passed) and increment by the constants defined above
  float radius = 0;
  float theta  = 0;
  float attenuationLength = 5.0;

  //TThread::Lock();
  //std::cout << "N0: " << N0 << std::endl;
  while (attenuationLength <= 40)
  {
    //std::cout << "mu: " << attenuationLength << std::endl;
    while ( radius <= p_diskRadius ) 
    {
      //std::cout << "radius: " << radius << std::endl;
      while ( theta < 360 ) 
      {
        // Since there are so many iterations here,
        // let's only fill the accumulator map 
        // with the max likelihood iterations
        // So find the N0 and attenL that correspond to ML
        float likelihood = ComputeLikelihood( radius, theta, N0, attenuationLength );
        //std::cout << likelihood << std::endl;
        if ( likelihood > p_likelihood) { p_likelihood = likelihood; p_N0 = N0; p_radius = radius; p_theta = theta; p_attenuationLength = attenuationLength; } 
        
        theta = theta + p_thetaBinSize;
      }
      theta = 0;
      radius = radius + p_radiusBinSize;
    }
    radius = 0;
    attenuationLength = attenuationLength + p_attenuationLengthBinSize;
  }
  //std::cout << "ML " << p_likelihood << "  atten: " << p_attenuationLength << "  radius " << p_radius << "  theta " << p_theta << "  N0 " << p_N0 << std::endl;
}

void Analyzer::RunReco(SiPMToTriggerMap& sipmToTriggerMap, const SiPMInfoMap& sipmInfoMap, const Configuration& config)
{
  // Initialize
  p_thetaBinSize  = config.thetaBinSize;
  p_radiusBinSize = config.radiusBinSize;
  p_diskRadius    = config.diskRadius;
  p_attenuationLengthBinSize = config.attenuationLengthBinSize;
  p_beta          = 360/config.nSiPMs;
  p_nSiPMs        = config.nSiPMs;
  p_likelihood    = std::numeric_limits<float>::lowest();

 
  // Loop over different triggers
  //for (unsigned trigger = 1; trigger <= sipmToTriggerMap.find(1)->second.size(); trigger++)
  //{
    WheelReco(sipmToTriggerMap, sipmInfoMap, 1);
  //}

}

void Analyzer::FillAccumulatorMap()
{
  // Variables used for incrementing
  // Start with r = 0, theta = 0, (N0 is what's passed) and increment by the constants defined above
  float radius = 0;
  float theta  = 0;

  while ( radius <= p_diskRadius ) 
  {
    while ( theta < 360 ) 
    {
      // Compute likelihood using the N0 and attenL for ml already found
      float likelihood = ComputeLikelihood( radius, theta, p_N0, p_attenuationLength );

      std::map<std::string, float> params;
      params.emplace("radius", radius);         params.emplace("theta", theta); 
      params.emplace("likelihood", likelihood); params.emplace("attenuationLength", p_attenuationLength);
      p_accumulatorMap.emplace_back(p_N0, params);

      theta = theta + p_thetaBinSize;
    }
    theta = 0;
    radius = radius + p_radiusBinSize;
  }

  // Now sort for output purposes
  std::sort(p_accumulatorMap.begin(), p_accumulatorMap.end(), 
            [](const BinIndex& biL, const BinIndex& biR) 
            {return biL.second.find("likelihood")->second < biR.second.find("likelihood")->second;});

  /*for (const auto& id : p_accumulatorMap)
  {
    std::cout << "LH = "       << id.second.find("likelihood")->second << "  N0 = "    << id.first 
              << "  radius = " << id.second.find("radius")->second     << "  theta = " << id.second.find("theta")->second << std::endl;
  }*/
}

void Analyzer::WheelReco(SiPMToTriggerMap& sipmToTriggerMap, const SiPMInfoMap& sipmInfoMap, const unsigned& trigger)
{
   clock_t start = clock();

  //sipm_errors = dr.GetSiPMCountErrors(); 

  // First to cound number of photons in each hit
  // also N0 set to lower threshold to decrease computation time
  // and clear temp data vec for this trigger
  p_data.clear();
  unsigned N0 = CountPhotons(sipmToTriggerMap, sipmInfoMap, trigger);

  std::cout << "\nRunning MLE...\n" << std::endl;

  while ( N0 <= 1000 ) 
  {
    Handle(N0);
    N0++;
  }
 
  // We should have the ml now
  // Loop again and fill the accumulator map
  // for const N0 and attenuation length
  FillAccumulatorMap(); 

  // Now find the maximum
  float max = std::numeric_limits<float>::lowest();;
  unsigned ID(0);
  for (unsigned index = 0; index < p_accumulatorMap.size(); index++)
  {
    if (p_accumulatorMap[index].second.find("likelihood")->second > max) { max = p_accumulatorMap[index].second.find("likelihood")->second; ID = index; } 
  }
 
  std::cout << std::endl
            << "Max likelihood = " << p_accumulatorMap[ID].second.find("likelihood")->second        << std::endl
            << "Radius         = " << p_accumulatorMap[ID].second.find("radius")->second            << " cm\n"
            << "Theta          = " << p_accumulatorMap[ID].second.find("theta")->second             << " deg\n"
            << "N0             = " << p_accumulatorMap[ID].first                                    << " photons\n"
            << "Attenuation    = " << p_accumulatorMap[ID].second.find("attenuationLength")->second << " cm\n";
 
  clock_t end = clock();
  double duration = ((double) (end - start)) / CLOCKS_PER_SEC;
  std::cout << "\nRun time of " << duration << std::endl;

  // Add the parameters to our mle vec
  //std::vector<float> temp = {p_accumulatorMap[ID].second[1], p_accumulatorMap[ID].second[2], p_accumulatorMap[ID].second[3]};
  p_mleParameters.emplace_back(p_accumulatorMap[ID].first, p_accumulatorMap[ID].second);
}

unsigned Analyzer::CountPhotons(SiPMToTriggerMap& sipmToTriggerMap, const SiPMInfoMap& sipmInfoMap, const unsigned& trigger)
{
  // Get an average count to get a starting place for reco
  unsigned maxPhotons(0);
  
  // Loop over all the hits and set the n of photons
  for (auto& sipm : sipmToTriggerMap)
  {
    // Get the gain and bdown for this sipm
    const float& thisGain = sipmInfoMap.find(sipm.first)->second.gain;
    const float& thisBD   = sipmInfoMap.find(sipm.first)->second.breakdown;
    std::cout << "SiPM " << sipm.first << "  Gain = " << thisGain << "  Breakdown = " << thisBD << std::endl;
    unsigned max(0);
    for (auto& hit : sipm.second[trigger])
    {
      //std::cout << hit.hitHeight << "  " << hit.bias << "  " << thisBD << "\n";
      hit.nPhotons = std::round(hit.hitHeight/( thisGain*( hit.bias - thisBD ) ));   // hit.bias - thisBD
      //std::cout << hit.nPhotons << std::endl;
      if (hit.nPhotons > max)               max = hit.nPhotons;
      if (hit.nPhotons > maxPhotons) maxPhotons = hit.nPhotons;
    }
    // Store max in the data vec 
    p_data.emplace(sipm.first, max);
  }

  return maxPhotons;
}

float Analyzer::ComputeLikelihood(const float& r, const float& theta, const unsigned& N0, const float& attenuationLength)
{
  // Sum over terms ----> k_m*ln(lambda_m) - lambda_m - ln(k_m!)
  float sum = 0;
  for (int m = 1; m <= p_nSiPMs; m++) {
    float lambda_m = ComputeLambda( r, theta, N0, m, attenuationLength );
    //std::cout << "nPhotons: " << p_data.find(m)->second << "  lambda_m " << lambda_m << "   factorial " << TMath::Factorial(p_data.find(m)->second) << "  term ";
    float term = p_data.find(m)->second*log( lambda_m ) - lambda_m - log(TMath::Factorial(p_data.find(m)->second));
    //std::cout << term << std::endl;
    sum = sum + term;
  }
  return sum;
}

float Analyzer::ComputeLambda(const float& r, const float& theta, const unsigned& N0, const unsigned& m, const float& attenuationLength)
{
  ///lambda_m = N0e^( -r_m/attenuationLength )
  // Convert angle to radians 
  const double angleRad   = ( (m - 1)*p_beta - theta )*( TMath::Pi()/180 );
  const double r_mSquared = r*r + p_diskRadius*p_diskRadius - 2*r*p_diskRadius*TMath::Cos( angleRad );
  if (r_mSquared < 0) {
    std::cout << "Error! Square root of negative number!\n";
    return 0;
  }
  const double r_m = TMath::Sqrt( r_mSquared );
  const double lambda_m = N0*TMath::Exp( -r_m/attenuationLength );//(r_m);

  return lambda_m;
}
}
