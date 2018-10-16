//
// File: Analyzer.cxx
//
// Author: Hunter Sullivan
//
// Description: Structure to run reconstruction on sipm data.
//

#include "Analyzer.h"
#include "TMath.h"
#include "TF1.h"
#include "TF2.h"

namespace wheel {

Voxel::Voxel(const float& x, const float& y, const float& r, const float& theta)
 : m_x(x), m_y(y), m_r(r), m_theta(theta)
{}

Voxel::~Voxel()
{}

Analyzer::Analyzer()
{}

Analyzer::~Analyzer()
{}

void Analyzer::Reconstruct(SiPMToTriggerMap& sipmToTriggerMap, const SiPMInfoMap& sipmInfoMap, const Configuration& config, const unsigned& trigger)
{
  // Initialize voxels and containers
  std::cout << "Initializing reconstruction...\n";  
  Initialize(config);
  // Start reconstruction
  Reconstruct(sipmToTriggerMap, sipmInfoMap, trigger);
}

void Analyzer::Initialize(const Configuration& config)
{
  // Initialize
  m_nVoxels           = config.nVoxels;  
  m_diskRadius        = config.diskRadius;
  m_attenuationLength = config.attenuationLength;
  if (m_attenuationLength == 0) { std::cout << "Error. Please specify an attenuation length!\n"; std::exit(1); }
  m_beta              = 360/config.nSiPMs;
  m_nSiPMs            = config.nSiPMs;
  m_mlLikelihood      = std::numeric_limits<float>::lowest();
  m_mlRadius = 0; m_mlTheta = 0; m_mlN0 = 0; 
  m_data.clear();
  // Create the voxels
  InitVoxelList();
}

void Analyzer::InitVoxelList()
{
  // Since our geometry is a circle, loop through
  // x coordinate and y coordinate creating new voxels.
  // We can do this by building the first quadrant and duplicating.
  // Positions of voxels will be their centers.
  const unsigned firstQuadWidth = std::sqrt(m_nVoxels)/2;
  const float    increment      = 2*m_diskRadius/std::sqrt(m_nVoxels);
  // We will work our way out from the origin
  float x(increment/2);

  for (unsigned voxelCounterX = 1; voxelCounterX <= firstQuadWidth; voxelCounterX++)
  {
    // Reset the y position
    float y(increment/2);
    for (unsigned voxelCounterY = 1; voxelCounterY <= firstQuadWidth; voxelCounterY++)
    {
      // Compute r and theta 
      float r        = std::sqrt(x*x + y*y);
      float thetaDeg = std::abs(TMath::ASin(y/r))*(180/TMath::Pi());

      // Create our voxels
      Voxel v1( x,  y,  r,  thetaDeg);     Voxel v2(-x,  y,  r,  thetaDeg+90);
      Voxel v3(-x, -y,  r,  thetaDeg+180); Voxel v4( x, -y,  r,  thetaDeg+270);

      // Ignore if the voxel lies outside our ROI
      if (r >= m_diskRadius) continue;

      m_voxelList.emplace_back(v1); m_voxelList.emplace_back(v2);
      m_voxelList.emplace_back(v3); m_voxelList.emplace_back(v4);

      y += increment;
    }
    x += increment;
  }
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
  // We will cover the entire parameter space, with the cost of computation time
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
  std::cout << "Max likelihood = " << m_accumulatorMap.back().second.find("likelihood")->second  << std::endl
            << "X              = " << m_accumulatorMap.back().second.find("x")->second           << " cm\n"
            << "Y              = " << m_accumulatorMap.back().second.find("y")->second           << " cm\n"
            << "Radius         = " << m_accumulatorMap.back().second.find("radius")->second      << " cm\n"
            << "Theta          = " << m_accumulatorMap.back().second.find("theta")->second       << " deg\n" 
            << "N0             = " << m_accumulatorMap.back().first                              << " photons\n";

  ComputeConfidenceIntervals();
 
  clock_t end = clock();
  double duration = ((double) (end - start)) / CLOCKS_PER_SEC;
  std::cout << "Run time of " << duration << " s" << std::endl;

  // Add these parameters to our mle vec
  m_mleParameters.emplace_back(m_accumulatorMap.back().first, m_accumulatorMap.back().second);
}

void Analyzer::Handle(const unsigned& N0)
{
  for (const auto& voxel : m_voxelList)
  {
    // Likelihood for this parameter set
    float likelihood = ComputeLikelihood(voxel.R(), voxel.Theta(), N0);
    if (likelihood > m_mlLikelihood) { m_mlLikelihood = likelihood; m_mlN0 = N0; m_mlX = voxel.X(); m_mlY = voxel.Y(); m_mlRadius = voxel.R(); m_mlTheta = voxel.Theta(); } 
  }
}

float Analyzer::ComputeLikelihood(const float& r, const float& thetaDeg, const unsigned& N0)
{
  // Sum over terms ----> k_m*ln(lambda_m) - lambda_m - ln(k_m!)
  float sum = 0;
  for (int sipm = 1; sipm <= m_nSiPMs; sipm++) 
  {
    float lambda_m = ComputeLambda(r, thetaDeg, N0, sipm);
    //std::cout << "nPhotons: " << m_data.find(m)->second << "  lambda_m " << lambda_m << "   factorial " << TMath::Factorial(m_data.find(m)->second) << "  term ";
    float term = m_data.find(sipm)->second*log(lambda_m) - lambda_m - log(TMath::Factorial(m_data.find(sipm)->second));
    //std::cout << term << std::endl;
    sum = sum + term;
  }
  return sum;
}

float Analyzer::ComputeLambda(const float& r, const float& thetaDeg, const unsigned& N0, const unsigned& sipm)
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
  float radToDeg     = 180/TMath::Pi();
    
  float angleXYandSiPMRad = ((sipm - 1)*m_beta - thetaDeg)*(1/radToDeg);
  float sipmToXYSquared   = r*r + m_diskRadius*m_diskRadius - 2*r*m_diskRadius*TMath::Cos(angleXYandSiPMRad);
    
  // Safety here
  if (sipmToXYSquared < 0) { std::cout << "Error! Square root of negative number!\n"; return 0; } 
  float sipmToXY = std::sqrt(sipmToXYSquared);

  // Get the angle between sipmToXY and the normal for this sipm
  // Put a protection here
  float cosAngleSiPMToXYandSiPM(1);
  if (sipmToXY != 0) cosAngleSiPMToXYandSiPM = (-r*r + sipmToXYSquared + m_diskRadius*m_diskRadius)/(2*sipmToXY*m_diskRadius);
 
  float weight = N0*cosAngleSiPMToXYandSiPM*TMath::Exp(-sipmToXY/m_attenuationLength)/sipmToXY;
  //std::cout << "Weight = " << relativeWeight << "  at sipm " << sipm << " from x = " << x << " y = " << y << "  sipmToXY = " << sipmToXY << " angleXYandSIPM = " << angleXYandSiPMRad*180/TMath::Pi() << " beta = " << m_beta << "  thetaDeg = " << thetaDeg <<std::endl;
  if (weight < 0) { std::cout << "UH OH! WEIGHT < 0!!\n"; std::exit(1); }
  return weight;
}

void Analyzer::FillAccumulatorMap()
{
  for (const auto& voxel : m_voxelList)
  {
    // Likelihood for this parameter set
    float likelihood = ComputeLikelihood(voxel.R(), voxel.Theta(), m_mlN0);
    // Update our accumulator map
    std::map<std::string, float> params;
    params.emplace("x", voxel.X());
    params.emplace("y", voxel.Y());
    params.emplace("radius", voxel.R());         
    params.emplace("theta", voxel.Theta()); 
    params.emplace("likelihood", likelihood); 
    m_accumulatorMap.emplace_back(m_mlN0, params);
  }

  // Order according to likelihood
  std::sort(m_accumulatorMap.begin(), m_accumulatorMap.end(), 
            [](const BinIndex& biL, const BinIndex& biR) 
            {return biL.second.find("likelihood")->second < biR.second.find("likelihood")->second;});
}

void Analyzer::ComputeConfidenceIntervals()
{
  // Some notes:
  //
  // [x-, x+]  --> x +- c*1/sqrt(I)
  //
  // c = 2.576 --> 99%
  // c = 2.326 --> 98%
  // c = 1.96  --> 95%
  // c = 1.645 --> 90%
  //
  // Our equation is:
  //
  // I = sum_k [ (data[k]/lambda^2[k])(dlambda/dx)^2 - (1 - data[k]/lambda[k])(d^2lambda/dx2) ]
  //
  // Rememeber:
  // 
  // lambda ~ N0*cos(alpha)*Exp(-r/attenuationLength)/r
  //
 
  std::vector<float> sum(2,0);
  for (unsigned sipm = 1; sipm <= m_nSiPMs; sipm++)
  {
    // Define the 2D lambda function to evaluate 1st and 2nd derivatives
    // We need two of these to evaluate partials
    // Be careful here, variable is always x!
    std::string r_x              = "std::sqrt(x*x + "+std::to_string(m_mlX*m_mlX)+")";
    std::string thetaDeg_x       = "TMath::ACos("+std::to_string(m_mlX)+"/"+r_x+")*180/TMath::Pi()";
    std::string cosp_x           = "TMath::Cos( (("+std::to_string(sipm)+"-1)*"+std::to_string(m_beta)+" - "+thetaDeg_x+")*TMath::Pi()/180 )";
    std::string rBar_x           = "std::sqrt("+r_x+"*"+r_x+" + "+std::to_string(m_diskRadius)+"*"+std::to_string(m_diskRadius)+" - 2*"+r_x+"*"+std::to_string(m_diskRadius)+"*"+cosp_x+")";
    std::string cosAlpha_x       = "(("+rBar_x+"*"+rBar_x+" + "+std::to_string(m_diskRadius)+"*"+std::to_string(m_diskRadius)+" - "+r_x+"*"+r_x+")/(2*"+rBar_x+"*"+std::to_string(m_diskRadius)+"))";
    std::string lambdaXConstEquation = std::to_string(m_mlN0)+"*"+cosAlpha_x+"*TMath::Exp(-"+rBar_x+"/"+std::to_string(m_attenuationLength)+")/"+rBar_x; 

    /*std::cout << "r_x        = " << r_x        << "\n"
              << "thetaDeg_x = " << thetaDeg_x << "\n"
              << "cosp_x     = " << cosp_x     << "\n"
              << "rBar_x     = " << rBar_x     << "\n"
              << "cosAlpha_x = " << cosAlpha_x << "\n"
              << "lambdaXC   = " << lambdaXConstEquation << "\n\n";*/

    TF1 lambdaXConst("lambdaXConst", lambdaXConstEquation.c_str(), -m_diskRadius, m_diskRadius);
   
    std::string r_y            = "std::sqrt(x*x + "+std::to_string(m_mlY*m_mlY)+")";
    std::string thetaDeg_y     = "TMath::ACos(x/"+r_y+")*180/TMath::Pi()";
    std::string cosp_y         = "TMath::Cos( (("+std::to_string(sipm)+"-1)*"+std::to_string(m_beta)+" - "+thetaDeg_y+")*TMath::Pi()/180 )";
    std::string rBar_y         = "std::sqrt("+r_y+"*"+r_y+" + "+std::to_string(m_diskRadius)+"*"+std::to_string(m_diskRadius)+" - 2*"+r_y+"*"+std::to_string(m_diskRadius)+"*"+cosp_y+")";
    std::string cosAlpha_y     = "(("+rBar_y+"*"+rBar_y+" + "+std::to_string(m_diskRadius)+"*"+std::to_string(m_diskRadius)+" - "+r_y+"*"+r_y+")/(2*"+rBar_y+"*"+std::to_string(m_diskRadius)+"))";
    std::string lambdaYConstEquation = std::to_string(m_mlN0)+"*"+cosAlpha_y+"*TMath::Exp(-"+rBar_y+"/"+std::to_string(m_attenuationLength)+")/"+rBar_y; 

    TF1 lambdaYConst("lambdaYConst", lambdaYConstEquation.c_str(), -m_diskRadius, m_diskRadius);

    // Derivatives
    float dlambdadx   = lambdaYConst.Derivative(m_mlX);
    float d2lambdadx2 = lambdaYConst.Derivative2(m_mlX); 
    float dlambdady   = lambdaXConst.Derivative(m_mlY);
    float d2lambdady2 = lambdaXConst.Derivative2(m_mlY); 

    // We also need lambda in terms of x and y
    std::string r              = "std::sqrt(x*x + y*y)";
    std::string thetaDeg       = "TMath::ACos(x/"+r+")*180/TMath::Pi()";
    std::string cosp           = "TMath::Cos( (("+std::to_string(sipm)+"-1)*"+std::to_string(m_beta)+" - "+thetaDeg+")*TMath::Pi()/180 )";
    std::string rBar           = "std::sqrt("+r+"*"+r+" + "+std::to_string(m_diskRadius)+"*"+std::to_string(m_diskRadius)+" - 2*"+r+"*"+std::to_string(m_diskRadius)+"*"+cosp+")";
    std::string cosAlpha       = "(("+rBar+"*"+rBar+" + "+std::to_string(m_diskRadius)+"*"+std::to_string(m_diskRadius)+" - "+r+"*"+r+")/(2*"+rBar+"*"+std::to_string(m_diskRadius)+"))";
    std::string lambdaEquation = std::to_string(m_mlN0)+"*"+cosAlpha+"*TMath::Exp(-"+rBar+"/"+std::to_string(m_attenuationLength)+")/"+rBar; 

    TF2 lambda2D("lambda2D", lambdaEquation.c_str(), -m_diskRadius, m_diskRadius, -m_diskRadius, m_diskRadius);

    // Lambda evaluated at our ml estimate
    float lambda = lambda2D.Eval(m_mlX, m_mlY); 
  
    // Finally...
    std::cout << "Sum[0] = " << sum[0] << "  Sum[1] = " << sum[1] << "  ratio = " << (m_data[sipm-1]/lambda) << std::endl;

    sum[0] += (m_data[sipm-1]/(lambda*lambda))*(dlambdadx*dlambdadx) + d2lambdadx2*(1 - (m_data[sipm-1]/lambda));
    sum[1] += (m_data[sipm-1]/(lambda*lambda))*(dlambdady*dlambdady) + d2lambdady2*(1 - (m_data[sipm-1]/lambda));
  }

  std::cout << "Sum[0] = " << sum[0] << "  Sum[1] = " << sum[1] << std::endl;
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
}
