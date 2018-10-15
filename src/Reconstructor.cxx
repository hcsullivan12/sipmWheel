//
// File: Reconstructor.cxx
//
// Author: Hunter Sullivan
//
// Description: Reconstruction algorithm for sipm data.
//              Our disk will be segmented into voxels, where the 
//              corresponding weights, c_ij, represent the probability
//              that a photon leaving pixel j is detected by sipm i.
//              The likelihood is maximized using an EM algorithm.

#include "Reconstructor.h"
#include "TMath.h"
#include "TF2.h"
#include "TH2F.h"
#include "TFile.h"

namespace wheel {

Voxel::Voxel(const float& x, const float& y)
 : m_x(x), m_y(y)
{}

Voxel::~Voxel()
{}

Reconstructor::Reconstructor()
{}

Reconstructor::~Reconstructor()
{}

void Reconstructor::Reconstruct(SiPMToTriggerMap& sipmToTriggerMap, const SiPMInfoMap& sipmInfoMap, const Configuration& config, const unsigned& trigger)
{
  // Initialize voxels and containers
  std::cout << "Initializing reconstruction...\n";	
  Initialize(config);
  // Start reconstruction 
  std::cout << "Starting reconstruction...\n";
  Reconstruct(sipmToTriggerMap, sipmInfoMap, trigger);
  // Save the results
  SaveResults();
}

void Reconstructor::Initialize(const Configuration& config)
{
  // Useful parameters
  m_nSiPMs            = config.nSiPMs;
  m_beta              = 360.0/m_nSiPMs;
  m_diskRadius        = config.diskRadius;
  m_nVoxels           = config.nVoxels;
  m_attenuationLength = config.attenuationLength;
  m_recoOutputFile    = config.recoOutputFile;
  // Clear the data container
  m_data.clear();
  // Fill the voxel map
  InitVoxelList();
}

void Reconstructor::InitVoxelList()
{
  // This map will be filled with voxels.
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
      Voxel v1( x,  y); Voxel v2(-x,  y);
      Voxel v3(-x, -y); Voxel v4( x, -y);

      // Ignore if the voxel lies outside our ROI
      float dist = std::sqrt(x*x + y*y);
      if (dist >= m_diskRadius) continue;

      m_voxelList.emplace_back(v1); m_voxelList.emplace_back(v2);
      m_voxelList.emplace_back(v3); m_voxelList.emplace_back(v4);

      y += increment;
    }
    x += increment;
  }
  
  // Now we need to store the weights
  // First compute relative weights, then normalize to the whole disk
  std::vector<float> sums(m_voxelList.size(), 0);
  unsigned voxelCounter(0);
  for (auto& voxel : m_voxelList)
  {
    // Probability map for this voxel
    std::map<unsigned, float> probabilityMap;
    for (unsigned sipm = 1; sipm <= m_nSiPMs; sipm++)
    {
	    float weight = ComputeWeight(sipm, voxel.X(), voxel.Y());
	    sums[voxelCounter] += weight;
      probabilityMap.emplace(sipm, weight);
    }
    voxel.SetProbabilityMap(probabilityMap);
    voxelCounter++;
  }
  voxelCounter = 0;
  for (auto& voxel : m_voxelList)
  {
    // Probability map for this voxel
    auto probabilityMap = voxel.ProbabilityMap();
    for (unsigned sipm = 1; sipm <= m_nSiPMs; sipm++)
    {
      probabilityMap.find(sipm)->second = probabilityMap.find(sipm)->second/sums[voxelCounter];
    }
    voxel.SetProbabilityMap(probabilityMap);
    voxelCounter++;
  }
}

float Reconstructor::ComputeWeight(const unsigned& sipm, const float& x, const float& y)
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
  float r                 = std::sqrt(x*x + y*y);
  float thetaDeg          = std::abs(TMath::ASin(y/r))*180/TMath::Pi();
  // Handle theta convention
  // 2nd quadrant
  if (x < 0 && y > 0) thetaDeg += 90;
  // 3rd quadrant
  if (x < 0 && y < 0) thetaDeg += 180;
  // 4th quadrant
  if (x > 0 && y < 0) thetaDeg += 270;
  // Remember, we overlaid a square grid on our disk, so some of the voxels
  // we need to ignore
  if (r >= m_diskRadius) return 0; 
  float angleXYandSiPMRad = ( (sipm - 1)*m_beta - thetaDeg )*( TMath::Pi()/180 );
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
 
  float relativeWeight = 1*cosAngleSiPMToXYandSiPM*TMath::Exp(-sipmToXY/m_attenuationLength)/sipmToXY;
  //std::cout << "Weight = " << relativeWeight << "  at sipm " << sipm << " from x = " << x << " y = " << y << "  sipmToXY = " << sipmToXY << " angleXYandSIPM = " << angleXYandSiPMRad*180/TMath::Pi() << " beta = " << m_beta << "  thetaDeg = " << thetaDeg <<std::endl;
  if (relativeWeight < 0) { std::cout << "UH OH! WEIGHT < 0!!\n"; std::exit(1); }
  return relativeWeight;
}

void Reconstructor::Reconstruct(SiPMToTriggerMap& sipmToTriggerMap, const SiPMInfoMap& sipmInfoMap, const unsigned& trigger)
{
  // Start the timer for this trigger
  clock_t start = clock();

  // First count number of photons in each hit, returning the sipm
  // with the largest count
  const auto maxSiPM_count = InitData(sipmToTriggerMap, sipmInfoMap, trigger);
  // Make our initial estimate for the source intensities
  FirstEstimate(maxSiPM_count);
  
  // Okay, here we go...
  std::cout << "Running MLE...\n";
  
  // We have our first estimate, now let's make our next estimate
  unsigned iterator(1);
  //MakeNextEstimate(iterator); 
}

void Reconstructor::FirstEstimate(const std::pair<unsigned, unsigned>& maxSiPM_counts)
{
  // We will make our first guesses:
  //
  //  R:     diskRadius/2
  //  Theta: (maxSiPM-1)*beta  ///< This is just the angular position of maxSiPM
  // 
  const float rGuess     = m_diskRadius/2;
  const float thetaGuess = m_beta*(maxSiPM_counts.first-1); 
  // What the x and y coordinates?
  const float xGuess     = rGuess*TMath::Cos(thetaGuess*TMath::Pi()/180);
  const float yGuess     = rGuess*TMath::Sin(thetaGuess*TMath::Pi()/180); 

 // std::cout << "RGuess = " << rGuess << "   thetaGuess = " << thetaGuess << "  maxSiPM = " << maxSiPM_counts.first << "  counts = " << maxSiPM_counts.second << "\n";

  // Let's make a gaussian distribution centered around this position
  // The amplitude at the mean will be the maxSiPM_counts.second + 1
  std::string gaus2D = "[0]*TMath::Exp( -((x-[1])*(x-[1]) + (y-[2])*(y-[2]))/(2*[3]*[3]) )";
  TF2 intensityDist("intensityDist", gaus2D.c_str(), -m_diskRadius, m_diskRadius, -m_diskRadius, m_diskRadius);
  intensityDist.SetParameters(maxSiPM_counts.second+1, xGuess, yGuess, m_diskRadius/4);

  // Now update our voxels
  for (auto& voxel : m_voxelList)
  {
    float intensityGuess = intensityDist.Eval(voxel.X(), voxel.Y());
    voxel.SetOldIntensity(intensityGuess);
  }
}

void Reconstructor::MakeNextEstimate(unsigned& iterator)
{
  // Our next estimate for the intensities is constructed from two terms:
  //
  //   1) Numerator:    Sum_i { c_ij*OldEstimate_j*Counts_i/Sum_k { c_ik*OldEstimate_k} }
  //   2) Denominator:  Sum_i { c_ij }
  // 
  // c_ij is the probability that photon leaving voxel j reaches sipm i.
  // 
  // The idea here is to construct our new intenstiy estimates and then compare to
  // our old estimates. Then we can check for our "epsilon convergence."

  std::cout << "Iteration = " << iterator << std::endl;

  // Helpful parameters
  std::vector<float> numDenom(m_nSiPMs, 0);

  // Compute sum in numerator's denominator
  for (const auto& voxel : m_voxelList) 
  {
    const auto probMap = voxel.ProbabilityMap();
    for (unsigned sipm = 1; sipm <= m_nSiPMs; sipm++) numDenom[sipm-1] += probMap.find(sipm)->second*voxel.OldIntensity(); 
  }
  // Now we can loop over each voxel
  for (auto& voxel : m_voxelList)
  {
    float num(0), denom(0);
    // Numerator
    const auto probMap = voxel.ProbabilityMap();
    for (unsigned sipm = 1; sipm <= m_nSiPMs; sipm++)
    {
      num += probMap.find(sipm)->second*voxel.OldIntensity()*m_data[sipm-1]/numDenom[sipm-1];
    }
    // Denominator
    for (const auto& sipmWeight : probMap) denom += sipmWeight.second;
    
    float newEstimate = num/denom;
    //std::cout << "Num = " << num << "  denom = " << denom << std::endl;
    voxel.SetNewIntensity(newEstimate);
  }

  // Check for convergence
  float eps(0);
  CheckConvergence(eps); 

  // If we're not happy with this, re-estimate
  std::cout << "Epsilon = " << eps << std::endl;
  if (iterator < 20/*m_epsilonConvergence*/) 
  { 
    iterator++;
    // Update our estimates 
    for (auto& voxel : m_voxelList)
    {
      //std::cout << "X = " << voxel.X() << "  Y = " << voxel.Y() << "  Old = " << voxel.OldIntensity() << "  New = " << voxel.NewIntensity() << std::endl;
      voxel.SetOldIntensity(voxel.NewIntensity());
    }
    MakeNextEstimate(iterator); 
  }
}

void Reconstructor::CheckConvergence(float& eps)
{
  // This idea is that if our new estimate doesn't 
  // change much, then we've converged
  float sum(0);
  unsigned counter(0);
  for (const auto& voxel : m_voxelList)
  {
    counter++;
    sum += (voxel.OldIntensity() - voxel.NewIntensity())*(voxel.OldIntensity() - voxel.NewIntensity());
    //std::cout << "old = " << voxel.OldIntensity() << "  new = " << voxel.NewIntensity() << std::endl;
  }
  // Our metric is the average difference
  std::cout << sum << "  " << counter << std::endl;
  eps = sum/counter;
}

void Reconstructor::SaveResults()
{
  // Need to think about plots here

  // First plot: 2D histograms showing intensities 
  unsigned nBins(std::sqrt(m_nVoxels));
  TH2F intensities("intensities", "intensities", nBins, -m_diskRadius, m_diskRadius, nBins, -m_diskRadius, m_diskRadius);
  for (const auto& voxel : m_voxelList)
  {
    auto xBin = intensities.GetXaxis()->FindBin(voxel.X()); 
    auto yBin = intensities.GetYaxis()->FindBin(voxel.Y());
    intensities.SetBinContent(xBin, yBin, voxel.OldIntensity());
  }

  TFile outputFile(m_recoOutputFile.c_str(), "UPDATE");
  intensities.Write();
  outputFile.Close();
}

std::pair<unsigned, unsigned> Reconstructor::InitData(SiPMToTriggerMap& sipmToTriggerMap, const SiPMInfoMap& sipmInfoMap, const unsigned& trigger)
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
