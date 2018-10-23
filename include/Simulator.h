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

namespace wheel {

class Photon {

public:
  Photon();
  ~Photon();

  void UpdateStatus(const float&              weight, 
                    const std::vector<float>& unitMomentumVec, 
                    const std::vector<float>& lastPos,
                    const std::vector<float>& currentPos,
                    const std::vector<float>& nextPos) 
 {
   m_weight          = weight;
   m_unitMomentumVec = unitMomentumVec;
   m_lastPosition    = lastPos;
   m_currentPosition = currentPos; 
   m_nextPosition    = nextPos;
 };

  std::vector<float> CurrentPosition() const { return m_currentPosition; };
  std::vector<float> NextPosition()    const { return m_nextPosition; };

private:
  
  float              m_weight;                 ///< Current weight 
  std::vector<float> m_unitMomentumVec;        ///< Momentum unit vector (cartesian)
  std::vector<float> m_lastPosition;           ///< Position at last step (x,y,z)
  std::vector<float> m_currentPosition;        ///< Position at current step (x,y,z)
  std::vector<float> m_nextPosition;           ///< Next projected position
};

class SiPM {

public:
  SiPM(const std::vector<float>& pos, const float& area);
  ~SiPM();

private:

  float              m_totalWeight;            ///< Sum of weights 
  std::vector<float> m_position;               ///< Position (r, theta)
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
  void  CalculateNextPosition(std::vector<float>& nextPos, const std::vector<float>& currentPos, const std::vector<float>& unitMomentumVec);
  float CalculateDistance(const std::vector<float>& currentPos, const std::vector<float>& nextPos);
  bool  Step(Photon& photon);
  
  std::string              m_simulateOutputPath;     ///< Output path for results
  unsigned                 m_nSiPMs;                 ///< Number of sipms
  float                    m_sipmArea;               ///< Active area of sipm
  unsigned                 m_nPhotonsToLaunch;       ///< Number of photons to launch
  float                    m_diskRadius;             ///< Disk radius
  float                    m_diskThickness;          ///< Disk thickness
  std::vector<float>       m_sourcePosition;         ///< Source position
  float                    m_tpbEmissionPeak;        ///< TPB emission peak (nm)
  float                    m_indexRefractionDisk;    ///< Index of refraction of disk
  float                    m_indexRefractionEnv;     ///< Index of refraction of sorrounding environment
  float                    m_terminationThreshold;   ///< Termination threshold
  float                    m_surfaceAbsorptionCoeff; ///< Surface absorption coeffecient (0 < s < 1)
  float                    m_bulkAbsorption;         ///< Bulk absorption length
  std::map<unsigned, SiPM> m_sipmMap;                ///< Map to sipm info
  TRandom                  m_rg;                     ///< Random number generator
};
}

#endif
