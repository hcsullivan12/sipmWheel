// 
// File: Simulator.cxx
//
// Author: Hunter Sullivan
//
// Discription: Structure to simulate photon propogation in sipmwheel.
//

#include "Simulator.h"
#include "TMath.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TF2.h"
#include "Reconstructor.h"

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
  m_simulateOutputPath     = config.recoOutputPath;
  m_recoAnaTreePath        = config.recoAnaTreePath;
  m_nSiPMs                 = config.nSiPMs;
  m_N0                     = config.N0;
  m_reconstruct            = config.reconstruct;
  m_diskRadius             = config.diskRadius;
  m_diskThickness          = config.diskThickness;
  m_sourcePosition.push_back(config.sourcePosition[0]);
  m_sourcePosition.push_back(config.sourcePosition[1]);
  m_sourceSigma            = config.sourceSigma;
  m_smearSigma             = config.smearSigma;
  m_tpbEmissionPeak        = config.tpbEmissionPeak;
  m_indexRefractionDisk    = config.indexRefractionDisk;
  m_indexRefractionEnv     = config.indexRefractionEnv;
  m_permittivityDisk       = config.permittivityDisk;
  m_permittivityEnv        = config.permittivityEnv;
  m_nPhotonsToLaunch       = config.nPhotonsToLaunch;
  m_terminationThreshold   = config.terminationThreshold;
  m_sipmArea               = config.sipmArea;
  m_bulkAttenuation        = config.bulkAttenuation;
  m_bulkAbsorption         = config.bulkAbsorption;
  m_surfaceAbsorptionCoeff = config.surfaceAbsorptionCoeff;
  // Set seed
  m_stepGenerator.SetSeed(0);
  m_rGenerator.SetSeed(0);
  m_phiGenerator.SetSeed(0);
  m_momentumGenerator.SetSeed(0);
  m_smearGenerator.SetSeed(0);
  m_rouletteGenerator.SetSeed(0); 
  // Initialize histograms
  m_stepHist.SetName("Step Size");
  m_stepHist.SetBins(100, 0, 200);
  m_intensityProfileHist.SetName("Intensity profile");
  m_intensityProfileHist.SetBins(500, -m_diskRadius-1, m_diskRadius+1, 500, -m_diskRadius-1, m_diskRadius+1);
  m_totalStepsHist.SetName("Total Steps");
  m_totalStepsHist.SetBins(500, 0, 1000);
  m_depositionHist.SetName("E Deposition");
  m_depositionHist.SetBins(102, -m_diskRadius-1, m_diskRadius+1, 102, -m_diskRadius-1, m_diskRadius+1); 
}

Simulator::~Simulator()
{}

void Simulator::Simulate(const Configuration& config)
{
  // **Steps for this algorithm**
  //
  //  1) Emit photons from source position! Assume Gaussian beam profile
  //     centered on source position, and tpb emission direction assumed
  //     uniform.
  //  2) Get step size from exponential distribution parameterized by
  //     bulk absorption length. If step size smaller than distance to
  //     to next boundary, absorb, then step again! Otherwise, reflect
  //     off next boundary.
  //  3) Check angle of incidence with critical angle. Assume rough surface.
  //     If angle of incidence >= critical, R=1. Otherwise, calculate 
  //     reflection coefficient. Compute outgoing momentum, then
  //     step!
  //  

  // Initialize sipms
  Initialize();

  std::vector<TGraph*> stepGraphs;

  // Commence emission!
  clock_t start = clock();
  for (unsigned photonCounter = 1; photonCounter <= m_nPhotonsToLaunch; photonCounter++)
  {
    if (photonCounter%1000 == 0) std::cout << "Photon " << photonCounter << std::endl;

    // Initialize photon!
    Photon photon;

    // Emit!
    Emit(photon);
    // If it missed the disk
    if(Terminate(photon)) continue;

    photon.XSteps().push_back(photon.CurrentPosition()[0]);
    photon.YSteps().push_back(photon.CurrentPosition()[1]);
    photon.ZSteps().push_back(photon.CurrentPosition()[2]);

    // Series of safety checks
    auto currentPos   = photon.CurrentPosition();
    auto nextBoundary = photon.NextBoundary(); 
    float rCurrent    = std::sqrt(currentPos[0]*currentPos[0] + currentPos[1]*currentPos[1]);
    float rBoundary   = std::sqrt(nextBoundary[0]*nextBoundary[0] + nextBoundary[1]*nextBoundary[1]);
    if (rCurrent        > (m_diskRadius+0.0001))    std::cout << "Error. Inital position greater than disk radius!\n";
    if (rBoundary       > (m_diskRadius+0.0001))    std::cout << "Error. Next boundary position r > disk radius!\n";
    if (nextBoundary[2] < (0-0.0001))               std::cout << "Error. Next boundary position z < 0!\n";
    if (nextBoundary[2] > (m_diskThickness+0.0001)) std::cout << "Error. Next boundary position z > diskThickness!\n";

    // Begin tracing!
    // Be carful here!
    unsigned stepNumber(0);
    Step(photon, stepNumber);
    m_totalStepsHist.Fill(stepNumber);
    
    auto stepsX = photon.XSteps();
    auto stepsY = photon.YSteps();
    auto stepsZ = photon.ZSteps();

    // Make a TPolyLine for only one photon, for visualization purposes
    if (photonCounter == 1)
    {
      TPolyLine3D polyLine(stepsX.size(), &stepsX[0], &stepsY[0], &stepsZ[0]);
      TFile f(m_simulateOutputPath.c_str(), "RECREATE");
      int LineBins = 100;
      TPolyLine3D *l1 = 0;
      TPolyLine3D *l2 = 0;
      l1 = new TPolyLine3D(LineBins);
      l2 = new TPolyLine3D(LineBins);
      for (int e = 0; e < LineBins; ++e) 
      {
        Double_t Angle = e*2*TMath::Pi()/LineBins;
        l1->SetPoint(e, m_diskRadius*TMath::Cos(Angle), m_diskRadius*TMath::Sin(Angle), 0);
        l2->SetPoint(e, m_diskRadius*TMath::Cos(Angle), m_diskRadius*TMath::Sin(Angle), m_diskThickness);
      }
      l1->SetLineColor(2);
      l2->SetLineColor(2);
      TCanvas c1("3D Photon Tracker", "3D Photon Tracker", 800, 800);
      l1->Draw("");
      l2->Draw("same");
      polyLine.Draw("l same");
      c1.Write();
      delete l1;
      delete l2;
    }

    if (photonCounter > 500) continue;
    stepGraphs.push_back(new TGraph(stepsX.size(), &stepsX[0], &stepsY[0]));
  }

  clock_t end = clock();
  double duration = ((double)(end - start))/CLOCKS_PER_SEC;

  std::cout << "\nSimulation run time of " << duration << " s\n" << std::endl;

  HandleSiPMInfo();
  MakePlots(stepGraphs);

  // Now reconstruction
  if (m_reconstruct)
  {
    std::cout << "\nInitializing reconstruction...\n";

    // First we need to fill the data
    std::map<unsigned, unsigned> data;
    unsigned max(0);
    for (unsigned sipm = 1; sipm <= m_nSiPMs; sipm++)
    {
      unsigned nPE = m_sipmMap.find(sipm)->second.NPhotons();
      if (nPE > max) max = nPE;
      data.emplace(sipm, nPE);
      std::cout << "SIPM " << sipm << " --> " << nPE << std::endl;
    }
    /*data.emplace(1,100);
    data.emplace(2,70);
    data.emplace(3,50);
    data.emplace(4,20);
    data.emplace(5,10);
    data.emplace(6,20);
    data.emplace(7,50);
    data.emplace(8,70);*/

    // Initialize the main work horse for reconstruction
    wheel::Reconstructor reconstructor;    
    reconstructor.Initialize(config);
    reconstructor.SetData(data);
    // Begin!
    reconstructor.Reconstruct(max); 
    // Plots
    reconstructor.MakePlot(1);
    // Add the true point
    UpdatePlot();
  }
}

void Simulator::Initialize()
{
  // This assumes that the sipms are seperated 
  // by the same angle and centered on the disk 
  // at each site

  // Initialize our sipms
  for (unsigned sipmID = 1; sipmID <= m_nSiPMs; sipmID++)
  {
    // Position of this sipm
    float thetaDeg = (sipmID - 1)*360.0/m_nSiPMs;
    float r        = m_diskRadius;
    // Remember, we're using cylindrical coordinates here
    std::vector<float> pos = {r, thetaDeg, m_diskThickness/2};

    SiPM sipm(pos, m_sipmArea);
    sipm.SetTotalWeight(0.0);
    m_sipmMap.emplace(sipmID, sipm);
  }
}

void Simulator::HandleSiPMInfo()
{
  TCanvas c("LightYield", "Light Yield", 800, 800);
  c.SetTheta(90);
  c.SetPhi(0);

  TH2D* h = new TH2D("h", "h", 2*8, 0, 360, 20, 0, 20);
  
  float sum(0);
  for (const auto& sipm : m_sipmMap) sum += sipm.second.TotalWeight();

    
  for (auto& sipm : m_sipmMap)
  {
    std::cout << "SiPM " << sipm.first << " captured " << (sipm.second.TotalWeight()/sum)*100 << "%!\n";
   
    unsigned nPE = std::round(m_N0*(sipm.second.TotalWeight()/sum));
    sipm.second.SetNPhotons(nPE);
    h->SetBinContent(2*(sipm.first -1) +1, m_diskRadius, sipm.second.TotalWeight()/sum);
  }

  h->Draw("lego2 pol");
  c.Update();
  TFile f(m_simulateOutputPath.c_str(), "UPDATE");
  c.Write();
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
  // We will assume our source has a gaussian intensity profile
  // We will also assume tpb has isotropic (Pz<0) emission 
 
  // Note: Bottom of disk is the xy plane
  // **TODO: Pick the wavelength**

  // First generate r and phi for this photon
  float r   = m_sourceSigma*std::sqrt(-std::log(m_rGenerator.Uniform(0,1)));
  float phi = 2*TMath::Pi()*m_phiGenerator.Uniform(0,1); 
  
  // What's x and y 
  float incidentX = r*TMath::Cos(phi);
  float incidentY = r*TMath::Sin(phi);  
  
  // Generate tpb emission direction
  // Only accept if z <= 0
  Double_t xTemp(0), yTemp(0), zTemp(1);
  while (zTemp > 0) m_momentumGenerator.Sphere(xTemp, yTemp, zTemp, 1);
  // To get rid of warnings
  float xDir = xTemp; float yDir = yTemp; float zDir = zTemp;

  // Get XYZ position of our light source
  // Assuming light source is pointing in -z direction
  float lsX(0), lsY(0), lsZ(m_diskThickness);
  ConvertToCartesian(lsX, lsY, m_sourcePosition[0], m_sourcePosition[1]);

  // Translate
  std::vector<float> unitMomentum    = {xDir, yDir, zDir};
  std::vector<float> lastPos         = {lsX+incidentX, lsY+incidentY, lsZ};
  std::vector<float> currentPos      = lastPos;

  // Assuming no interaction, where is the next surface?
  CalculateNextBoundary(photon, currentPos, unitMomentum); 

  // If this photon missed our disk, we missed it
  if (CalculateR(currentPos) > (m_diskRadius-0.01))
  {
    photon.SetWeight(0);
    photon.SetUnitMomentum(unitMomentum);
    photon.SetLastPosition(lastPos);
    photon.SetCurrentPosition(currentPos);
    return;
  }
 
  photon.SetWeight(1.0);
  photon.SetUnitMomentum(unitMomentum);
  photon.SetLastPosition(lastPos);
  photon.SetCurrentPosition(currentPos);
 
  // Update our histogram
  m_intensityProfileHist.Fill(currentPos[0], currentPos[1]);
}

void Simulator::CalculateNextBoundary(Photon& photon, const std::vector<float>& currentPos, const std::vector<float>& unitMomentum)
{ 
  // We have to check for the smallest epsilon > 0
  // Hitting bottom and top
  // We need to take the maximum of the z epsilons since
  // we may already be on the top or bottom
  float zEpsilonFloor   = (0 - currentPos[2])/unitMomentum[2]; 
  float zEpsilonCeiling = (m_diskThickness - currentPos[2])/unitMomentum[2];
  float zEpsilon = std::max(zEpsilonFloor, zEpsilonCeiling);
  // Hitting side
  float B = unitMomentum[0]*currentPos[0] + unitMomentum[1]*currentPos[1];
  float A = unitMomentum[0]*unitMomentum[0] + unitMomentum[1]*unitMomentum[1];
  float C = currentPos[0]*currentPos[0] + currentPos[1]*currentPos[1] - m_diskRadius*m_diskRadius;
  float rEpsilon = (-B + std::sqrt(B*B - A*C))/A;

  // Now sort the epsilons
  std::vector<float> epsilons = {zEpsilon, rEpsilon};
  std::sort(epsilons.begin(), epsilons.end(), [](const float& left, const float& right) { return left < right;});
  // Take the smallest > 0
  auto it = epsilons.begin();
  while (*it < 0) it++;
  float stepEpsilon = *it;

  // Update
  std::vector<float> nextBoundary(3,0);
  nextBoundary[0] = currentPos[0] + stepEpsilon*unitMomentum[0]; 
  nextBoundary[1] = currentPos[1] + stepEpsilon*unitMomentum[1]; 
  nextBoundary[2] = currentPos[2] + stepEpsilon*unitMomentum[2]; 

  // WE DON"T KNOW WHICH BOUNDARY YET!!!

  // Handle floating point errors
  float r = CalculateR(nextBoundary);
  float eps(0);
  while ((m_diskRadius-0.03) < r)
  {
    nextBoundary[0] = currentPos[0] + (stepEpsilon+eps)*unitMomentum[0];
    nextBoundary[1] = currentPos[1] + (stepEpsilon+eps)*unitMomentum[1];
    nextBoundary[2] = currentPos[2] + (stepEpsilon+eps)*unitMomentum[2];
    eps = eps - 0.01;
    r = CalculateR(nextBoundary);
  } 
  
  if (r > (m_diskRadius+0.01)) std::cout << "HOLDUP50\n"; 

  // Set the next boundary
  photon.SetNextBoundary(nextBoundary);

  // What boundary is this?
  if ((m_diskRadius-0.05) < r    && r < (m_diskRadius+0.05))    { photon.SetBoundary("side");   return; }
  if ((0-0.01) < nextBoundary[2] && nextBoundary[2] < (0+0.01)) { photon.SetBoundary("bottom"); return; }
  
  // Otherise
  photon.SetBoundary("top");
}

void Simulator::Step(Photon& photon, unsigned& stepNumber)
{
  // Have we hit an sipm yet?
  if (Captured(photon)) return;
  stepNumber++;

  // First bulk propogation
  auto stepSize = m_stepGenerator.Exp(m_bulkAttenuation);
  // What's the distance to the next boundary?
  auto distToBoundary = CalculateDistance(photon.CurrentPosition(), photon.NextBoundary());  
  m_stepHist.Fill(stepSize);

  float temp3 = CalculateR(photon.CurrentPosition());
  if (temp3 > (m_diskRadius+0.001)) std::cout << "HOLDUP30 " << temp3 << "\n";

  if (stepSize <= distToBoundary)
  {
    //std::cout << "\n\nHERE\n";
    // We're being absorbed!
    // Update the weight
    float weight   = photon.Weight();
    float lossFrac = weight*(m_bulkAbsorption/m_bulkAttenuation); 
    weight = weight - lossFrac;
    photon.SetWeight(weight);
    if (std::isnan(photon.Weight())) std::cout << weight << std::endl;

    // We have to step again!
    // Find next position
    std::vector<float> nextPos(3,0);
    nextPos[0] = photon.CurrentPosition()[0] + stepSize*photon.UnitMomentum()[0];
    nextPos[1] = photon.CurrentPosition()[1] + stepSize*photon.UnitMomentum()[1];
    nextPos[2] = photon.CurrentPosition()[2] + stepSize*photon.UnitMomentum()[2];
 
    // We're neglecting scattering, so we don't have to update the
    // next boundary or unit momentum
    photon.SetLastPosition(photon.CurrentPosition());
    photon.SetCurrentPosition(nextPos);
 
    // Update our deposition plot
    float temp = CalculateR(photon.CurrentPosition());
    if (temp > (m_diskRadius-0.001)) std::cout << "HOLDUP10 " << temp << "\n";
    auto xBin = m_depositionHist.GetXaxis()->FindBin(photon.CurrentPosition()[0]);
    auto yBin = m_depositionHist.GetYaxis()->FindBin(photon.CurrentPosition()[1]);
    auto cont = m_depositionHist.GetBinContent(xBin, yBin);
    m_depositionHist.SetBinContent(xBin, yBin, cont + lossFrac);
    
    photon.XSteps().push_back(photon.CurrentPosition()[0]);
    photon.YSteps().push_back(photon.CurrentPosition()[1]);
    photon.ZSteps().push_back(photon.CurrentPosition()[2]);

    // Step!
    Step(photon, stepNumber);
    return;
  }

  // We've hit a boundary!
  // Update current and last position
  photon.SetLastPosition(photon.CurrentPosition());
  photon.SetCurrentPosition(photon.NextBoundary());
  photon.XSteps().push_back(photon.CurrentPosition()[0]);
  photon.YSteps().push_back(photon.CurrentPosition()[1]);
  photon.ZSteps().push_back(photon.CurrentPosition()[2]);

  // We are assuming here that the sipms are 
  // optically coupled to the disk, this way 
  // we don't have to worry about transmission.
  // Handle reflection
  float temp = CalculateR(photon.CurrentPosition());
  if (temp > (m_diskRadius+0.001)) std::cout << "HOLDUP20 " << temp << "\n";
  if (std::isnan(photon.Weight())) std::cout << "her3e80 = " << std::endl;
  Reflect(photon);

  photon.XSteps().push_back(photon.CurrentPosition()[0]);
  photon.YSteps().push_back(photon.CurrentPosition()[1]);
  photon.ZSteps().push_back(photon.CurrentPosition()[2]);

  // Check the weight
  // Safety
  if (Terminate(photon)) return;

  // Recursive call
  Step(photon, stepNumber);
}

bool Simulator::Captured(Photon& photon)
{
  // We just need to check if our photon 
  // lies within the area of the sipm
  for (auto& sipmCount : m_sipmMap)
  {
    auto& sipm = sipmCount.second;
    // Define the area 
    float sideLength = std::sqrt(sipm.Area());
    float halfArc    = sideLength/2.0;
    float dAlphaDeg  = (halfArc/sipm.Position()[0])*180.0/TMath::Pi();

    float rEpsilon = 0.05;
    std::vector<float> rLimits        = {sipm.Position()[0]-rEpsilon,  sipm.Position()[0]+rEpsilon};
    std::vector<float> thetaLimitsDeg = {sipm.Position()[1]-dAlphaDeg, sipm.Position()[1]+dAlphaDeg};
       
    // Convert the photon's position to polar
    float pR(0), pThetaDeg(0);
    ConvertToPolar(pR, pThetaDeg, photon.CurrentPosition()[0], photon.CurrentPosition()[1]);
    // Subtlety here! Theta goes 0-360 > 0!
    // If sipm1, thetaLimitsDeg[0] < 0
    // Maybe not the best solution here, but we should handle sipm1 seperately
    if (sipmCount.first == 1) 
    {
      // Check if this sipm captured our photon
      if (rLimits[0] < pR && pR < rLimits[1])
      {
        thetaLimitsDeg[0] = thetaLimitsDeg[0] + 360;
        if (thetaLimitsDeg[0] < pThetaDeg || pThetaDeg < thetaLimitsDeg[1])
        {
          // Captured!
          float sipmWeight = sipm.TotalWeight();
          if (std::isnan(photon.Weight()) || std::isnan(sipmWeight)) std::cout << "SiPM " << sipmCount.first << "  " << sipmWeight << " photon weight " << photon.Weight() << std::endl;
          sipm.SetTotalWeight(sipmWeight + photon.Weight());
          photon.SetWeight(0.0);

          // Update our deposition plot
          float temp = CalculateR(photon.CurrentPosition());
          if (temp > (m_diskRadius+0.001)) std::cout << "HOLDUP2 " << temp << "\n";

          auto xBin = m_depositionHist.GetXaxis()->FindBin(photon.CurrentPosition()[0]);
          auto yBin = m_depositionHist.GetYaxis()->FindBin(photon.CurrentPosition()[1]);
          auto cont = m_depositionHist.GetBinContent(xBin, yBin);
          m_depositionHist.SetBinContent(xBin, yBin, sipm.TotalWeight());

          return true;
        }
      }
    }
    else 
    {
      if (thetaLimitsDeg[0] < 0) std::cout << "ERROR!\n";
      // Check if this sipm captured our photon
      if (rLimits[0] < pR && pR < rLimits[1]) 
      {
        if (thetaLimitsDeg[0] < pThetaDeg && pThetaDeg < thetaLimitsDeg[1])
        {
          // Captured!
          float sipmWeight = sipm.TotalWeight();
          if (std::isnan(photon.Weight()) || std::isnan(sipmWeight)) std::cout << "SiPM " << sipmCount.first << "  " << sipmWeight << "  photon weight " << photon.Weight() << std::endl;
          sipm.SetTotalWeight(sipmWeight + photon.Weight());
          photon.SetWeight(0.0);

          // Update our deposition plot
          float temp = CalculateR(photon.CurrentPosition());
          if (temp > (m_diskRadius+0.001)) std::cout << "HOLDUP3 " << temp << "\n";

          auto xBin = m_depositionHist.GetXaxis()->FindBin(photon.CurrentPosition()[0]);
          auto yBin = m_depositionHist.GetYaxis()->FindBin(photon.CurrentPosition()[1]);
          auto cont = m_depositionHist.GetBinContent(xBin, yBin);
          m_depositionHist.SetBinContent(xBin, yBin, sipm.TotalWeight());

          return true;
        }
      }
    }
  }
  // If it were captured, we wouldn't have made it here
  return false;
}

float Simulator::CalculateDistance(const std::vector<float>& currentPos, const std::vector<float>& nextPos)
{
  float delta[3] = {currentPos[0]-nextPos[0], currentPos[1]-nextPos[1], currentPos[2]-nextPos[2]};
  return std::sqrt(delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2]);
}

void Simulator::Reflect(Photon& photon)
{
  if (std::isnan(photon.Weight())) std::cout << "her3e110 = " << std::endl;

  // Have we hit an sipm yet?
  if (Captured(photon)) return;

  if (std::isnan(photon.Weight())) std::cout << "her3e100 = " << std::endl;

  // We have two mechanisms which attenuate our photons
  //  1) Loss of light from partial reflection (critical angle) 
  //  2) Absorption on surface 
  //

  // 1) 
  // First, let's calculate the angle of reflection
  if (m_indexRefractionEnv/m_indexRefractionDisk > 1) {std::cout << "Error. IndexRefractionEnv/indexRefractionDisk > 1!\n"; exit(1); }
  const float criticalAngleDeg  = TMath::ASin(m_indexRefractionEnv/m_indexRefractionDisk)*180/TMath::Pi();
  const auto  incidentPhotonDir = photon.UnitMomentum();

  // Now we have to know which surface we've hit
  // Top/bottom, Side
  std::vector<float> outgoingPhotonDir(3,0);
  std::vector<float> unitVec(3,0);

  // It's easier here if we smear the angle between zUnit and the surface norm
  if (photon.Boundary() == "top")    unitVec = {0, 0,-1};
  if (photon.Boundary() == "bottom") unitVec = {0, 0, 1};
  if (photon.Boundary() == "side")   unitVec = {-photon.CurrentPosition()[0], -photon.CurrentPosition()[1], 0};
  float temp1 = CalculateR(photon.CurrentPosition());
  if (temp1 > (m_diskRadius+0.001)) std::cout << "HOLDUP15 " << temp1 << "\n";

  // Normalize
  float mag = std::sqrt(unitVec[0]*unitVec[0] + unitVec[1]*unitVec[1] + unitVec[2]*unitVec[2]);
  unitVec[0] = unitVec[0]/mag;
  unitVec[1] = unitVec[1]/mag;
  unitVec[2] = unitVec[2]/mag;

  // Make a copy
  auto surfaceUnit = unitVec;
  // Smear the angle 
  const float smear1Deg = m_smearGenerator.Gaus(0, m_smearSigma);
  const float smear2Deg = m_smearGenerator.Gaus(0, m_smearSigma);
  const float smear3Deg = m_smearGenerator.Gaus(0, m_smearSigma);
  // Rotate surfaceUnit
  float sinSmear1 = TMath::Sin(smear1Deg*TMath::Pi()/180);
  float cosSmear1 = TMath::Cos(smear1Deg*TMath::Pi()/180);
  float sinSmear2 = TMath::Sin(smear2Deg*TMath::Pi()/180);
  float cosSmear2 = TMath::Cos(smear2Deg*TMath::Pi()/180);
  float sinSmear3 = TMath::Sin(smear3Deg*TMath::Pi()/180);
  float cosSmear3 = TMath::Cos(smear3Deg*TMath::Pi()/180);

  std::vector<float> surfaceUnitRot(3,0);
  surfaceUnitRot[0]   = cosSmear2*cosSmear3*surfaceUnit[0] - sinSmear3*cosSmear2*surfaceUnit[1] + sinSmear2*surfaceUnit[2];
  surfaceUnitRot[1]   = (sinSmear1*sinSmear2*cosSmear3 + cosSmear1*sinSmear3)*surfaceUnit[0] + 
                        (-sinSmear1*sinSmear2*sinSmear3 + cosSmear1*cosSmear3)*surfaceUnit[1] - sinSmear1*cosSmear2*surfaceUnit[2];
  surfaceUnitRot[2]   = (-sinSmear2*cosSmear1*cosSmear3 + sinSmear1*sinSmear3)*surfaceUnit[0] + 
                        (sinSmear2*cosSmear1*sinSmear3  + sinSmear1*cosSmear3)*surfaceUnit[1] + cosSmear1*cosSmear2*surfaceUnit[2];

  // Safety here, make sure our rotated normal is still a unit vector
  mag = std::sqrt(surfaceUnitRot[0]*surfaceUnitRot[0] + surfaceUnitRot[1]*surfaceUnitRot[1] + surfaceUnitRot[2]*surfaceUnitRot[2]);
  surfaceUnitRot[0] = surfaceUnitRot[0]/mag;
  surfaceUnitRot[1] = surfaceUnitRot[1]/mag;
  surfaceUnitRot[2] = surfaceUnitRot[2]/mag;

  // Compute the outgoing vector now
  // v and n are both unit vectors!
  float vDotn = incidentPhotonDir[0]*surfaceUnitRot[0] + incidentPhotonDir[1]*surfaceUnitRot[1] + incidentPhotonDir[2]*surfaceUnitRot[2];
  outgoingPhotonDir[0] = incidentPhotonDir[0] - 2*vDotn*surfaceUnitRot[0];
  outgoingPhotonDir[1] = incidentPhotonDir[1] - 2*vDotn*surfaceUnitRot[1];
  outgoingPhotonDir[2] = incidentPhotonDir[2] - 2*vDotn*surfaceUnitRot[2];
  // Update!
  photon.SetUnitMomentum(outgoingPhotonDir);

  // What's the angle of incidence?
  float cosIncidentAngle = std::abs(vDotn);
  float sinIncidentAngle = std::sqrt(1.001 - cosIncidentAngle*cosIncidentAngle);
  float incidentAngleDeg = TMath::ACos(cosIncidentAngle)*180/TMath::Pi();
      
  if (std::isnan(photon.Weight())) std::cout << "her3e150 = " << cosIncidentAngle << " " << incidentAngleDeg <<  std::endl;

  
  // If it's below the critical angle, attenuate
  // Be careful! 
  if (incidentAngleDeg < (criticalAngleDeg - 0.1))
  {
    // We need the reflection coeffecient here
    // Compute R for each polarization
    float radicalTerm = (m_indexRefractionDisk/m_indexRefractionEnv)*sinIncidentAngle;
    float b           = (m_indexRefractionDisk/m_indexRefractionEnv)*(m_indexRefractionDisk/m_indexRefractionEnv)*(m_permittivityEnv/m_permittivityDisk);
    
    float numPerpen = m_indexRefractionDisk*cosIncidentAngle - b*m_indexRefractionEnv*std::sqrt(1 - radicalTerm*radicalTerm);
    float denPerpen = m_indexRefractionDisk*cosIncidentAngle + b*m_indexRefractionEnv*std::sqrt(1 - radicalTerm*radicalTerm);

    float numPara   = b*m_indexRefractionEnv*m_indexRefractionEnv*cosIncidentAngle - m_indexRefractionEnv*m_indexRefractionDisk*std::sqrt(1 - radicalTerm*radicalTerm);
    float denPara   = b*m_indexRefractionEnv*m_indexRefractionEnv*cosIncidentAngle + m_indexRefractionEnv*m_indexRefractionDisk*std::sqrt(1 - radicalTerm*radicalTerm);

    float R_Perpen = (numPerpen/denPerpen)*(numPerpen/denPerpen);
    float R_Para   = (numPara/denPara)*(numPara/denPara);

    // Take the average
    float R = 0.5*(R_Perpen + R_Para);
    //std::cout << R << std::endl;

    float weight = photon.Weight();
    if (std::isnan(weight)) std::cout << "R = " << R << std::endl;
    photon.SetWeight(R*weight); 
 
    // E deposition
    auto xBin = m_depositionHist.GetXaxis()->FindBin(photon.CurrentPosition()[0]);
    auto yBin = m_depositionHist.GetYaxis()->FindBin(photon.CurrentPosition()[1]);
    auto cont = m_depositionHist.GetBinContent(xBin, yBin);
    m_depositionHist.SetBinContent(xBin, yBin, cont + (1-R)*weight);

    if (Terminate(photon)) return;
    if (std::isnan(photon.Weight())) std::cout << "her3e60 = " << std::endl;
  }

  // 2)
  // This occurs for every photon
  // Update our deposition plot
  float temp = CalculateR(photon.CurrentPosition()); 
  if (temp > (m_diskRadius+0.001)) std::cout << "HOLDUP1 " << temp << "\n";

  auto xBin = m_depositionHist.GetXaxis()->FindBin(photon.CurrentPosition()[0]);
  auto yBin = m_depositionHist.GetYaxis()->FindBin(photon.CurrentPosition()[1]);
  auto cont = m_depositionHist.GetBinContent(xBin, yBin);
  //if (temp > (m_diskRadius - 1))  std::cout << m_surfaceAbsorptionCoeff*photon.Weight() << std::endl;
  m_depositionHist.SetBinContent(xBin, yBin, cont + m_surfaceAbsorptionCoeff*photon.Weight());

  if (std::isnan(photon.Weight())) std::cout << "her1e " << photon.Weight() << std::endl;
  float weight = (1.0 - m_surfaceAbsorptionCoeff)*photon.Weight();
  photon.SetWeight(weight); 
  if (std::isnan(photon.Weight())) std::cout << "here " << weight << std::endl;
 
  // What is the next boundary?
  CalculateNextBoundary(photon, photon.CurrentPosition(), photon.UnitMomentum()); 
}

float Simulator::CalculateR(const std::vector<float>& pos)
{
  return std::sqrt(pos[0]*pos[0] + pos[1]*pos[1]);
}

bool Simulator::Terminate(Photon& photon)
{
  // Roulette procedure
  float rand = m_rouletteGenerator.Uniform(0,1);
  if (photon.Weight() < m_terminationThreshold) 
  {
    if (rand <= 0.1) 
    {
      float weight = photon.Weight();
      photon.SetWeight(weight/0.1);
      return false;
    }
    else { photon.SetWeight(0.0); return true; }
  } 
}

void Simulator::MakePlots(std::vector<TGraph*> stepGraphs)
{
  TFile f(m_simulateOutputPath.c_str(), "UPDATE");

  m_stepHist.Write();
  m_totalStepsHist.Write();

  TCanvas c("intensityProfile", "Intesnity Profile", 800, 800);
  TF2 diskf("diskf", "x*x + y*y", -m_diskRadius-1, m_diskRadius+1, -m_diskRadius-1, m_diskRadius+1);
  double contours[2] = {0, m_diskRadius*m_diskRadius};
  diskf.SetContour(2, contours);
  diskf.SetLineColor(1);
  diskf.Draw();
  m_intensityProfileHist.Draw("same");
  c.Write();
  TCanvas c2("stepHistory", "Step History", 800, 800);
  diskf.Draw();
  for (unsigned g = 0; g < stepGraphs.size(); g++)
  {
    stepGraphs[g]->SetMarkerStyle(7);
    stepGraphs[g]->SetMarkerSize(1);
    stepGraphs[g]->Draw("same");
  }
  c2.Write();

  m_depositionHist.Write();
}

void Simulator::UpdatePlot()
{
  // Open the file we created
  TFile f(m_simulateOutputPath.c_str(), "UPDATE");

  TCanvas *c1 = nullptr;
  std::string name = "logLikelihood_CL_" + std::to_string(1);
  f.GetObject(name.c_str(), c1);

  if (!c1) std::cout << "Error. Could not find object named " << name << " in " << m_simulateOutputPath.c_str() << std::endl;

  // Draw the true point
  float x(0), y(0);
  ConvertToCartesian(x, y, m_sourcePosition[0], m_sourcePosition[1]);
  TMarker trueXY(x, y, 8);
  trueXY.SetMarkerSize(2);
  trueXY.SetMarkerColor(2);

  c1->Draw();
  trueXY.Draw("same");

  c1->Write();
  delete c1;
}
}
