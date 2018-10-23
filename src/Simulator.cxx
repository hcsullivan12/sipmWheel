// 
// File: Simulator.cxx
//
// Author: Hunter Sullivan
//
// Discription: Structure to simulate photon propogation in sipmewheel.
//

#include "Simulator.h"
#include "TMath.h"

namespace wheel {

Photon::Photon()
{}

Photon::~Photon()
{}

SiPM::SiPM(const std::vector<float>& pos, const float& area)
 : m_totalWeight(0), m_position(pos), m_area(area)
{};

SiPM::~SiPM()
{}

Simulator::Simulator(const Configuration& config)
{
  // So we don't have to keep passing config
  m_simulateOutputPath     = config.simulateOutputPath;
  m_nSiPMs                 = config.nSiPMs;
  m_diskRadius             = config.diskRadius;
  m_diskThickness          = config.diskThickness;
  m_sourcePosition.push_back(config.sourcePosition[0]);
  m_sourcePosition.push_back(config.sourcePosition[1]);
  m_tpbEmissionPeak        = config.tpbEmissionPeak;
  m_indexRefractionDisk    = config.indexRefractionDisk;
  m_indexRefractionEnv     = config.indexRefractionEnv;
  m_nPhotonsToLaunch       = config.nPhotonsToLaunch;
  m_terminationThreshold   = config.terminationThreshold;
  m_sipmArea               = config.sipmArea;
  m_bulkAbsorption         = config.bulkAbsorption;
  m_surfaceAbsorptionCoeff = config.surfaceAbsorptionCoeff;
  // Set seed
  m_rg.SetSeed(0);
}

Simulator::~Simulator()
{}

void Simulator::Simulate()
{
  // **Steps for this algorithm**
  //
  //  1) Emit photons from source position according to some angular distribution.
  //     (Maybe take photon wavelength from TPB emission distribution as well!, but how to include into propogation?)
  //     This initializes directions for each photon.
  //  2) Get step size from exponontial distribution with bulk absorption length.
  //     if step size < distance to disk surface, we lost it!
  //  3) Check angle of incidence with critical angle
  //     (Maybe take angle of incidence as gaussian to model imperfect surface, quasispecular)
  //     if angle of incidence >= critical, R=1
  //     else calculate Reflection coeffiecent 
  //     REPEAT!
  //
  //     Check the positions of photons at each step
  //     if photon position = area of sipm, terminate
  //     if photon weight < threshold, we lost it!
  //
  //     The sum of the photon weights captured is the gives the probability of detection at the sipm
  //  

  // Initialize sipms
  Initialize();

  // Commence emission!
  for (unsigned photonCounter = 1; photonCounter <= m_nPhotonsToLaunch; photonCounter++)
  {
    // Initialize photon!
    Photon photon;
    // Emit!
    Emit(photon);
    // Step!
    if (!Step(photon)) continue;
    // Reflect!
    if (!Reflect(photon);

  }
  
}

void Simulator::Initialize()
{
  // Initialize our sipms
  for (unsigned sipmID = 1; sipmID <= m_nSiPMs; sipmID++)
  {
    // Position of this sipm
    float thetaDeg = (sipmID - 1)*360.0/m_nSiPMs;
    float r        = m_diskRadius;
    std::vector<float> pos = {r, thetaDeg};

    SiPM sipm(pos, m_sipmArea);
    m_sipmMap.emplace(sipmID, sipm);
  }
}

void Simulator::ConvertToPolar(float& r, float& thetaDeg, const float& x, const float& y)
{
  r        = std::sqrt(x*x + y*y);
  thetaDeg = TMath::ASin(std::abs(y/r))*180/TMath::Pi();
 
  // Handle theta convention
  if (x < 0 && y > 0) thetaDeg = 180 - thetaDeg;
  if (x < 0 && y < 0) thetaDeg = 180 + thetaDeg;
  if (x > 0 && y < 0) thetaDeg = 360 - thetaDeg;
}

void Simulator::ConvertToCartesian(float& x, float& y, const float& r, const float& thetaDeg)
{
  x = r*TMath::Cos(thetaDeg*TMath::Pi()/180);
  y = r*TMath::Sin(thetaDeg*TMath::Pi()/180);
}

void Simulator::Emit(Photon& photon)
{
  // We will assume uniform distribution for direction
  // Bottom of disk is the xy plane
  // **TODO: Pick the wavelength**
  
  // Only accept if z <= 0
  Double_t xTemp(0), yTemp(0), zTemp(1);
  while (zTemp > 0) m_rg.Sphere(xTemp, yTemp, zTemp, 1);
 
  // Exchange for floats (get rid of warning message)
  float x = xTemp; float y = yTemp; float z = zTemp;

  // Get XYZ position of our light source
  // Assuming light source is pointing in -z direction
  float lsX(0), lsY(0), lsZ(m_diskThickness);
  ConvertToCartesian(lsX, lsY, m_sourcePosition[0], m_sourcePosition[1]);

  std::vector<float> unitMomentumVec = {x,     y,   z};
  std::vector<float> lastPos         = {lsX, lsY, lsZ};
  std::vector<float> currentPos      = {lsX, lsY, lsZ}; 
  std::vector<float> nextPos         = currentPos;

  // Assuming no interaction, where is the next surface?
  CalculateNextPosition(nextPos, currentPos, unitMomentumVec); 

  photon.UpdateStatus(1.0, unitMomentumVec, lastPos, currentPos, nextPos);
}

void Simulator::CalculateNextPosition(std::vector<float>& nextPos, const std::vector<float>& currentPos, const std::vector<float>& unitMomentumVec)
{ 
  // Keep moving until we've crossed boundary
  float epsilon(0.001), r(0);

  while (r < m_diskRadius && nextPos[2] >= 0 && nextPos[2] <= m_diskThickness)
  {
    // Increment
    nextPos[0] = nextPos[0] + epsilon*unitMomentumVec[0];
    nextPos[1] = nextPos[1] + epsilon*unitMomentumVec[1];
    nextPos[2] = nextPos[2] + epsilon*unitMomentumVec[2];

    // Recalculate r
    r = std::sqrt(nextPos[0]*nextPos[0] + nextPos[1]*nextPos[1]);
  }
}

bool Simulator::Step(Photon& photon)
{
  // First: bulk propogation
  auto stepSize = m_rg.Exp(m_bulkAbsorption);
  auto dist     = CalculateDistance(photon.CurrentPosition(), photon.NextPosition());  

  if (dist >= stepSize) return false;
  return true;
}

float Simulator::CalculateDistance(const std::vector<float>& currentPos, const std::vector<float>& nextPos)
{
  float delta[3] = {currentPos[0]-nextPos[0], currentPos[1]-nextPos[1], currentPos[2]-nextPos[2]};
  return std::sqrt(delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2]);
}
}
