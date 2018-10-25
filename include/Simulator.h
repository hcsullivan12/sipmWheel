// 
// File: Simulator.h
//
// Author: Hunter Sullivan
//
// Discription: Structure to simulate photon propogation in sipmwheel.
//

#ifndef SIMULATOR_H
#define SIMULATOR_H

#include "Utilities.h"
#include <fstream>
#include <iostream>
#include "TRandom.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TPolyLine3D.h"

namespace wheel {

class Photon {

public:
  Photon();
  ~Photon();

  void UpdateStatus(const double&             weight, 
                    const std::vector<float>& unitMomentum, 
                    const std::vector<float>& lastPos,
                    const std::vector<float>& currentPos,
                    const std::vector<float>& nextBoundary) 
 {
   m_weight          = weight;
   m_unitMomentum    = unitMomentum;
   m_lastPosition    = lastPos;
   m_currentPosition = currentPos; 
   m_nextBoundary    = nextBoundary;
 };

  void SetWeight(const double& w)                        { m_weight          = w; };
  void SetCurrentPosition(const std::vector<float>& pos) { m_currentPosition = pos; };
  void SetLastPosition(const std::vector<float>& pos)    { m_lastPosition    = pos; };
  void SetNextBoundary(const std::vector<float>& pos)    { m_nextBoundary    = pos; };
  void SetUnitMomentum(const std::vector<float>& pos)    { m_unitMomentum    = pos; };
  void SetBoundary(const std::string& b)                 { m_boundary        = b; };

  std::vector<float> CurrentPosition() const { return m_currentPosition; };
  std::vector<float> NextBoundary()    const { return m_nextBoundary; };
  std::vector<float> UnitMomentum()    const { return m_unitMomentum; };
  double             Weight()          const { return m_weight; };
  std::string        Boundary()        const { return m_boundary; };
  std::vector<float>& XSteps()               { return m_xSteps; };
  std::vector<float>& YSteps()               { return m_ySteps; };
  std::vector<float>& ZSteps()               { return m_zSteps; };

private:
  
  double             m_weight;                 ///< Current weight 
  std::vector<float> m_unitMomentum;           ///< Momentum unit vector (cartesian)
  std::vector<float> m_lastPosition;           ///< Position at last step (x,y,z)
  std::vector<float> m_currentPosition;        ///< Position at current step (x,y,z)
  std::vector<float> m_nextBoundary;           ///< Next projected boundary position
  std::string        m_boundary;               ///< Current boundary (side, top, bottom)
  std::vector<float> m_xSteps;                 ///< History of steps X coordinate
  std::vector<float> m_ySteps;                 ///< History of steps Y coordinate
  std::vector<float> m_zSteps;                 ///< History of steps Z coordinate
};

class SiPM {

public:
  SiPM(const std::vector<float>& pos, const float& area);
  ~SiPM();

  const float              Area()        const { return m_area; };
  const std::vector<float> Position()    const { return m_position; }; 
  const double             TotalWeight() const { return m_totalWeight; };

  void SetTotalWeight(const double& w) { m_totalWeight = w; };


private:

  double             m_totalWeight;            ///< Sum of weights 
  std::vector<float> m_position;               ///< Position (r, thetaDeg, z)
  float              m_area;                   ///< Active area of sipm
};

class Simulator {

public:
  Simulator(const Configuration& config);
  ~Simulator();

  void Simulate();
   
private:

  void  Initialize();
  void  ConvertToPolar(float& r, float& thetaDeg, const float& x, const float& y);
  void  ConvertToCartesian(float& x, float& y, const float& r, const float& thetaDeg); 
  void  Emit(Photon& photon);
  void  CalculateNextBoundary(Photon& photon, const std::vector<float>& currentPos, const std::vector<float>& unitMomentumVec);
  float CalculateDistance(const std::vector<float>& currentPos, const std::vector<float>& nextPos);
  void  Step(Photon& photon, unsigned& stepNumber);
  void  Reflect(Photon& photon); 
  bool  Captured(Photon& photon);
  void  MakePlots(std::vector<TGraph*> stepGraphs);
  void  HandleSiPMInfo();
 
  std::string              m_simulateOutputPath;     ///< Output path for results
  unsigned                 m_nSiPMs;                 ///< Number of sipms
  float                    m_sipmArea;               ///< Active area of sipm
  unsigned                 m_nPhotonsToLaunch;       ///< Number of photons to launch
  float                    m_diskRadius;             ///< Disk radius
  float                    m_diskThickness;          ///< Disk thickness
  std::vector<float>       m_sourcePosition;         ///< Source position
  float                    m_sourceSigma;            ///< Spread in light source position
  float                    m_smearSigma;             ///< Smearing for incident angle
  float                    m_tpbEmissionPeak;        ///< TPB emission peak (nm)
  float                    m_indexRefractionDisk;    ///< Index of refraction of disk
  float                    m_indexRefractionEnv;     ///< Index of refraction of sorrounding environment
  double                   m_terminationThreshold;   ///< Termination threshold
  float                    m_surfaceAbsorptionCoeff; ///< Surface absorption coeffecient (0 < s < 1)
  float                    m_bulkAttenuation;        ///< Bulk attenuation length
  float                    m_bulkAbsorption;         ///< Bulk absorption length
  std::map<unsigned, SiPM> m_sipmMap;                ///< Map to sipm info
  TH1F                     m_stepHist;               ///< Histogram of steps
  TH2F                     m_intensityProfileHist;   ///< Histogram for intensity profile
  TH1F                     m_totalStepsHist;         ///< Hitsogram for total steps
  TRandom                  m_stepGenerator;          ///< Step generator
  TRandom                  m_rGenerator;             ///< Used for intensity profile
  TRandom                  m_phiGenerator;           ///< Used for intensity profile
  TRandom                  m_momentumGenerator;      ///< Generates tpb emission
  TRandom                  m_smearGenerator;         ///< Generates smearing
};
}

#endif
