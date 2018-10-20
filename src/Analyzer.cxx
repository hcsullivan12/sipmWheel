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
#include "TH2.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TBackCompFitter.h"

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
  std::cout << "\nInitializing reconstruction...\n";  
  Initialize(config);
  // Start reconstruction
  Reconstruct(sipmToTriggerMap, sipmInfoMap, trigger);
  // Make a useful plot
  MakePlot();
}

void Analyzer::Initialize(const Configuration& config)
{
  // Initialize
  m_maxIterations     = config.maxIterations;
  if (m_maxIterations == 0) { std::cout << "Error. Please specify maximum number of iterations!\n"; std::exit(1); }
  m_nVoxels           = config.nVoxels;  
  m_diskRadius        = config.diskRadius;
  if (m_diskRadius == 0) { std::cout << "Error. Please specify disk radius!\n"; std::exit(1); }
  m_attenuationLength = config.attenuationLength;
  if (m_attenuationLength == 0) { std::cout << "Error. Please specify an attenuation length!\n"; std::exit(1); }
  if (config.nSiPMs == 0) { std::cout << "Error. Please specify number of SiPMs!\n"; std::exit(1); }
  m_beta              = 360/config.nSiPMs;
  m_nSiPMs            = config.nSiPMs;
  m_mlLogLikelihood   = std::numeric_limits<double>::lowest();
  m_mlRadius = 0; m_mlTheta = 0; m_mlN0 = 0; 
  m_recoOutputFile    = config.recoOutputFile;
  m_data.clear();
  m_doCI = false;
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

  std::cout << "\nRunning MLE...\n";
  // Start main loop
  // We will cover the entire parameter space, at the expense of computation time
  // Once N0 is calculated, we'll use Newton-Raphson to better approximate x and y
  while ( N0 <= 300 ) 
  {
    Handle(N0);
    N0++;
  }
 
  // We should have the ml now
  std::cout << "Max log likelihood = " << m_mlLogLikelihood << std::endl
            << "X                  = " << m_mlX             << " cm\n"
            << "Y                  = " << m_mlY             << " cm\n"
            << "Radius             = " << m_mlRadius        << " cm\n"
            << "Theta              = " << m_mlTheta         << " deg\n" 
            << "N0                 = " << m_mlN0            << " photons\n";

  // Refine the estimation
  unsigned iterator(0);
  std::cout << "\nRefining estimate...\n";
  RefineEstimate(iterator);

  std::cout << "Max log likelihood = " << m_mlLogLikelihood << std::endl
            << "X                  = " << m_mlX             << " cm\n"
            << "Y                  = " << m_mlY             << " cm\n"
            << "Radius             = " << m_mlRadius        << " cm\n"
            << "Theta              = " << m_mlTheta         << " deg\n";

  // Estimate the confidence interval
  //ComputeConfidenceIntervals();
 
  clock_t end = clock();
  double duration = ((double) (end - start)) / CLOCKS_PER_SEC;
  std::cout << "Run time of " << duration << " s" << std::endl;
}

void Analyzer::Handle(const unsigned& N0)
{
  for (const auto& voxel : m_voxelList)
  {
    // Log likelihood for this parameter set
    double logLikelihood = ComputeLogLikelihood(voxel.X(), voxel.Y(), N0);
    if (logLikelihood > m_mlLogLikelihood) { m_mlLogLikelihood = logLikelihood; m_mlN0 = N0; m_mlX = voxel.X(); m_mlY = voxel.Y(); m_mlRadius = voxel.R(); m_mlTheta = voxel.Theta(); } 
  }
}

double Analyzer::ComputeLogLikelihood(const float& x, const float& y, const unsigned& N0)
{
  // Convert x and y to polar coordinates 
  float r        = 0;
  float thetaDeg = 0;
  ConvertToPolar(r, thetaDeg, x, y);  

  // Sum over terms ----> k_m*ln(lambda_m) - lambda_m - ln(k_m!)
  double sum = 0;
  for (int sipm = 1; sipm <= m_nSiPMs; sipm++) 
  {
    float lambda_m = ComputeLambda(r, thetaDeg, N0, sipm);
    //std::cout << "nPhotons: " << m_data.find(m)->second << "  lambda_m " << lambda_m << "   factorial " << TMath::Factorial(m_data.find(m)->second) << "  term ";
    double term = m_data.find(sipm)->second*log(lambda_m) - lambda_m - log(TMath::Factorial(m_data.find(sipm)->second));
    //std::cout << term << std::endl;
    sum = sum + term;
  }
  return sum;
}

void Analyzer::ConvertToPolar(float& r, float& thetaDeg, const float& x, const float& y)
{
  r        = std::sqrt(x*x + y*y);
  thetaDeg = TMath::ASin(std::abs(y/r))*180/TMath::Pi();

  // Handle theta convention
  if (x < 0 && y > 0) thetaDeg = 180 - thetaDeg;
  if (x < 0 && y < 0) thetaDeg = 180 + thetaDeg;
  if (x > 0 && y < 0) thetaDeg = 360 - thetaDeg;
}

float Analyzer::ComputeLambda(const float& r, const float& thetaDeg, const unsigned& N0, const unsigned& sipm)
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

void Analyzer::RefineEstimate(unsigned& iterator)
{
  // Make sure we haven't done this too many times
  iterator++;
  if (iterator > m_maxIterations) 
  {
    std::cout << "NOTE: Was not able to refine estimate in " << m_maxIterations << " iterations!" << std::endl;
    return;
  }

  // We will apply the Newton-Raphson method to
  // refine our estimate! We will take N0 as that
  // estimated from the "Handle" stage.
  //
  // We will use:
  //
  // V_n+1 = V_n - [l''(V_n)]^-1 * l'(V_n)
  //

  // We first need to calculate our derivatives
  // As h --> 0...
  const float h = 0.001;

  // dl/dx
  double l1   = ComputeLogLikelihood(m_mlX+h, m_mlY, m_mlN0);
  double dldx = (l1 - m_mlLogLikelihood)/h;

  // dl/dy
  double l2 = ComputeLogLikelihood(m_mlX, m_mlY+h, m_mlN0);
  double dldy = (l2 - m_mlLogLikelihood)/h;

  // d2l/dx2
  double l4 = ComputeLogLikelihood(m_mlX-h, m_mlY, m_mlN0);
  double l5 = ComputeLogLikelihood(m_mlX+h, m_mlY, m_mlN0);
  double d2ldx2 = (l4 - 2*m_mlLogLikelihood + l5)/(h*h);

  // d2l/dy2
  double l6 = ComputeLogLikelihood(m_mlX, m_mlY-h, m_mlN0);
  double l7 = ComputeLogLikelihood(m_mlX, m_mlY+h, m_mlN0);
  double d2ldy2 = (l6 - 2*m_mlLogLikelihood + l7)/(h*h);

  // d2l/dxdy
  double l10 = ComputeLogLikelihood(m_mlX+h, m_mlY+h, m_mlN0);
  double l11 = ComputeLogLikelihood(m_mlX+h, m_mlY-h, m_mlN0);
  double l12 = ComputeLogLikelihood(m_mlX-h, m_mlY+h, m_mlN0);
  double l13 = ComputeLogLikelihood(m_mlX-h, m_mlY-h, m_mlN0);
  double d2ldxdy = (l10 - l11 - l12 + l13)/(h*h);

  //std::cout << "2nd Derivatives => " << d2ldx2 << " " << d2ldy2 << " " << d2ldxdy << "\n";

  // Now we can form our Fisher matrix
  TArrayD data1(2);
  data1[0] = dldx; data1[1] = dldy;

  TArrayD data2(4);
  data2[0] = d2ldx2;
  data2[1] = d2ldxdy;
  data2[2] = d2ldxdy;
  data2[3] = d2ldy2;

  TMatrixD lPVector(2,1);
  lPVector.SetMatrixArray(data1.GetArray());
  TMatrixD fisherMatrix(2,2);
  fisherMatrix.SetMatrixArray(data2.GetArray());

  // Try to invert the matrix
  TDecompLU lu(fisherMatrix);
  TMatrixD inverseFisherMatrix = fisherMatrix;
  if (!lu.Decompose()) 
  {
    std::cout << "Decomposition failed, matrix singular!!" << std::endl;
    return;
  }
  else lu.Invert(inverseFisherMatrix);

  // Multiply Fisher and vector
  auto m = inverseFisherMatrix*lPVector;

  // Our convention --> x, y
  float newGuessX     = m_mlX  - m[0][0];
  float newGuessY     = m_mlY  - m[1][0];
  float newGuessR(0), newGuessT(0);
  ConvertToPolar(newGuessR, newGuessT, newGuessX, newGuessY);
  
  // Make sure r is not > disk radius
  if (newGuessR >= m_diskRadius) { std::cout << "Caution: Estimate went out of bounds!\n"; return; }
  double newGuessMLogL = ComputeLogLikelihood(newGuessX, newGuessY, m_mlN0);

  // Update
  m_mlX             = newGuessX;
  m_mlY             = newGuessY;
  ConvertToPolar(m_mlRadius, m_mlTheta, m_mlX, m_mlY);
  m_mlLogLikelihood = newGuessMLogL;
 
  // We want the eigenvalues and eigenvectors
  TVectorD eigenvalues(2);
  TMatrixD eigenvectors = fisherMatrix.EigenVectors(eigenvalues);
  // Clear old memory
  m_fisherEigenvalues.clear();
  m_fisherEigenvectors.clear();
  // Store eigenvalues
  m_fisherEigenvalues.push_back(eigenvalues[0]);
  m_fisherEigenvalues.push_back(eigenvalues[1]);
  // Store eigenvectors
  std::vector<float> e1 = {eigenvectors[0][0], eigenvectors[0][1]};
  std::vector<float> e2 = {eigenvectors[1][0], eigenvectors[1][1]};
  m_fisherEigenvectors.push_back(e1);
  m_fisherEigenvectors.push_back(e2);
 
  // Check the eigenvalues 
  if (eigenvalues[0] < 0 && eigenvalues[1] < 0) 
  { 
    m_doCI = true; 
    // Store the diagonal sigmas
    m_sigmaXDiag = 1.000/std::sqrt(-m_fisherEigenvalues[0]);
    m_sigmaYDiag = 1.000/std::sqrt(-m_fisherEigenvalues[1]);

    // We want to rotate back to our old coordinates now
    // Our rotation matrix is formed by the eigenvectors
    TArrayD rotData(4);
    rotData[0] = e1[0];
    rotData[1] = e2[0];
    rotData[2] = e1[1];
    rotData[3] = e2[1];
    // Our diagonal matrix has 1/m_sigmaXDiag^2 and 1/m_sigmaYDiag^2
    TArrayD diagData(4);
    diagData[0] = 1/(m_sigmaXDiag*m_sigmaXDiag);
    diagData[1] = 0;
    diagData[2] = 0;
    diagData[3] = 1/(m_sigmaYDiag*m_sigmaYDiag);

    TMatrixD rotationMatrix(2,2);
    rotationMatrix.SetMatrixArray(rotData.GetArray());
    TMatrixD diagMatrix(2,2);
    diagMatrix.SetMatrixArray(diagData.GetArray());

    TDecompLU lu(rotationMatrix);
    TMatrixD inverseRotationMatrix = rotationMatrix;
    if (!lu.Decompose()) 
    {
      std::cout << "Decomposition failed, matrix singular!!" << std::endl;
      return;
    }
    else lu.Invert(inverseRotationMatrix);

    // Now we have the matrices we need
    // Calculate:
    //
    // Sigma^-1 = R*Diag*R^-1
    //
   
    TMatrixD sigmaInverseMatrix = rotationMatrix*diagMatrix*inverseRotationMatrix;
    // Store these
    m_sigmaInverse.clear();
    m_sigmaInverse.push_back(sigmaInverseMatrix[0][0]);  m_sigmaInverse.push_back(sigmaInverseMatrix[0][1]);
    m_sigmaInverse.push_back(sigmaInverseMatrix[1][0]);  m_sigmaInverse.push_back(sigmaInverseMatrix[1][1]);

    return; 
  }
  // Otherwise, keep searching
  RefineEstimate(iterator);
}
/*
void Analyzer::ComputeConfidenceIntervals()
{
  // Some notes:
  //
  // [x-, x+]  --> x +- c*1/sqrt(I)
  //
  // Here: c*1/sqrt(I) = m_delta.find(percent)->second
  //
  // We will appoximate derivatives using standard numerical techniques.
  // We will compute the single CI, diagonalize to compute the joint CI,
  // where:
  //
  // CI = m_mlEstimate +- c*1/sqrt(-l'')
  //

  // We cannot form CI for events where the eigenvalues of Fisher are > 0!
  if (m_fisherEigenvalues[0] > 0 || m_fisherEigenvalues[1] > 0) return;

  // Store the sigmas
  m_sigmaX = 1.000/std::sqrt(-m_fisherEigenvalues[0]);
  m_sigmaY = 1.000/std::sqrt(-m_fisherEigenvalues[1]);


  // Form the margins for each CL
  m_deltaX.emplace(95, 1.960/std::sqrt(-m_fisherEigenvalues[0]));
  m_deltaX.emplace(90, 1.645/std::sqrt(-m_fisherEigenvalues[0]));
  m_deltaX.emplace(68, 1.000/std::sqrt(-m_fisherEigenvalues[0]));

  m_deltaY.emplace(95, 1.960/std::sqrt(-m_fisherEigenvalues[1]));
  m_deltaY.emplace(90, 1.645/std::sqrt(-m_fisherEigenvalues[1])); 
  m_deltaY.emplace(68, 1.000/std::sqrt(-m_fisherEigenvalues[1]));

  std::cout << m_deltaX.find(90)->second << "  " << m_deltaY.find(90)->second << std::endl;

}*/

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
/*
Double_t Analyzer::JointCLFunc(Double_t *x, Double_t *par)
{
  // Rotate our coordinate system
  TVectorD newXY(2);
  newXY[0] = m_fisherEigenvectors[0][0]*m_mlX + m_fisherEigenvectors[1][0]*m_mlY;
  newXY[1] = m_fisherEigenvectors[0][1]*m_mlX + m_fisherEigenvectors[1][1]*m_mlY;

  //std::cout << "New X = " << newXY[0] << "  new Y = " << newXY[1] << std::endl;

  // Now form the equations 
  Double_t a      = m_deltaX.find(par[0])->second;
  Double_t b      = m_deltaY.find(par[0])->second;
  Double_t xShift = m_fisherEigenvectors[0][0]*x[0] + m_fisherEigenvectors[1][0]*x[1] - newXY[0];
  Double_t yShift = m_fisherEigenvectors[0][1]*x[0] + m_fisherEigenvectors[1][1]*x[1] - newXY[1];

  Double_t eq = (xShift*xShift)/(a*a) + (yShift*yShift)/(b*b);

  return eq;
}*/

std::string Analyzer::JointCLFunc(const float& cl)
{
  // Rotate our coordinate system
  TVectorD newXY(2);
  newXY[0] = m_fisherEigenvectors[0][0]*m_mlX + m_fisherEigenvectors[1][0]*m_mlY;
  newXY[1] = m_fisherEigenvectors[0][1]*m_mlX + m_fisherEigenvectors[1][1]*m_mlY;

  //std::cout << "New X = " << newXY[0] << "  new Y = " << newXY[1] << std::endl;

  // Now form the equations 
  float a      = m_deltaX.find(cl)->second;
  float b      = m_deltaY.find(cl)->second;
  std::cout << "a = " << a << " b = " << b << std::endl;

  std::string xShift = "(" + std::to_string(m_fisherEigenvectors[0][0])+"*x +" + std::to_string(m_fisherEigenvectors[1][0]) + "*y - " + std::to_string(newXY[0]) + ")";
  std::string yShift = "(" + std::to_string(m_fisherEigenvectors[0][1]) + "*x +" + std::to_string(m_fisherEigenvectors[1][1]) + "*y - " + std::to_string(newXY[1]) + ")";

  std::string eq = std::to_string(TMath::Exp(m_mlLogLikelihood)) + "*TMath::Exp( -("+xShift+"*"+xShift+")/(2*"+std::to_string(a*a)+") )*TMath::Exp( -("+yShift+"*"+yShift+")/(2*"+std::to_string(b*b)+") )";
  return eq;
}

void Analyzer::MakePlot()
{
  TFile f(m_recoOutputFile.c_str(), "UPDATE");

  // Likelihood distribution for m_mlN0
  TH2D likelihoodDist("likelihoodDist", "Likelihood Distribution", std::sqrt(m_nVoxels), -m_diskRadius, m_diskRadius, std::sqrt(m_nVoxels), -m_diskRadius, m_diskRadius); 
  
  for (const auto& voxel : m_voxelList)
  {
    // Log likelihood for this parameter set
    double logLikelihood = ComputeLogLikelihood(voxel.X(), voxel.Y(), m_mlN0);

    unsigned xBin = likelihoodDist.GetXaxis()->FindBin(voxel.X());
    unsigned yBin = likelihoodDist.GetXaxis()->FindBin(voxel.Y());
    likelihoodDist.SetBinContent(xBin, yBin, TMath::Exp(logLikelihood));
  } 
    
  // Overlay the single confidence intervals only if this events "passed" the test
  if (!m_doCI) 
  {
    std::cout << "\nNOTE: EIGENVALUES > 0. NOT DRAWING CI!\n";
    return;
  }

  // Form our equation
  std::string delX2    = std::to_string(m_sigmaInverse[0])+"*(x - [1])*(x - [1])";
  std::string delY2    = std::to_string(m_sigmaInverse[3])+"*(y - [2])*(y - [2])";
  std::string delXdelY = "2*"+std::to_string(m_sigmaInverse[1])+"*(x - [1])*(y - [2])";
  std::string exp      = delX2 + " + " + delXdelY + " + " + delY2;
  std::string eq       = "[0]*TMath::Exp( -1*("+exp+") )";
 
  //TF2 *bigaus = new TF2("bigaus", eq.c_str(), -m_diskRadius, m_diskRadius, -m_diskRadius, m_diskRadius);
  //bigaus->SetParameters(1.0, m_mlX, m_mlY);


  /*// Compute the rotation angle
  float dotProduct  = m_fisherEigenvectors[0][0]*m_mlX + m_fisherEigenvectors[0][1]*m_mlY;
  float cosThetaRot = dotProduct/m_mlRadius; 
  float sinThetaRot = std::sqrt(1 - cosThetaRot*cosThetaRot);
  
  // Compute the rotated means
  TVectorD newXY(2);
  newXY[0] =  m_mlX*cosThetaRot + m_mlY*sinThetaRot;
  newXY[1] = -m_mlX*sinThetaRot + m_mlY*cosThetaRot; 
 
  std::string xP    = "x*" +std::to_string(cosThetaRot)+" + y*"+std::to_string(sinThetaRot);
  std::string yP    = "-x*"+std::to_string(sinThetaRot)+" + y*"+std::to_string(cosThetaRot);
  std::string delX2 = "(("+xP+"-[1])*("+xP+"-[1]))/(2*[2]*[2])";
  std::string delY2 = "(("+yP+"-[3])*("+yP+"-[3]))/(2*[4]*[4])";

  std::string eq = "[0]*TMath::Exp( -"+delX2+" - "+delY2+" )";*/

  TF2 bigaus("bigaus", "bigaus", m_mlX-1, m_mlX+1, m_mlY-1, m_mlY+1);
  bigaus.SetParameters(TMath::Exp(m_mlLogLikelihood), m_mlX, m_sigmaXDiag, m_mlY, m_sigmaYDiag, 0.1);
  likelihoodDist.Fit(&bigaus, "NR"); 

  /*// Contours
  TBackCompFitter *fitter = ((TBackCompFitter *)(TVirtualFitter::GetFitter()));
  
  int par1 = 1; // index of mean X
  int par2 = 3; // index of mean Y

  TGraph *gr1 = new TGraph( 80 );
  gr1->SetTitle(";top mass [GeV];#alpha_{S}^{}");
  gr1->SetFillColor(kGreen); gr1->SetLineColor(kBlack);
  double cl1 = 0.6827; // 1 sigma
  fitter->Contour(par1, par2, gr1, cl1);
  gr1->SetPoint(gr1->GetN(), gr1->GetX()[0], gr1->GetY()[0]); // "close" it
 
  TGraph *gr2 = new TGraph( 80 );
  gr2->SetTitle(";top mass [GeV];#alpha_{S}^{}");
  gr2->SetFillColor(kYellow); gr2->SetLineColor(kBlack);
  double cl2 = 0.9545; // 2 sigma
  fitter->Contour(par1, par2, gr2, cl2);
  gr2->SetPoint(gr2->GetN(), gr2->GetX()[0], gr2->GetY()[0]); // "close" it*/

  TCanvas c1("c1", "c1", 800, 800);
  likelihoodDist.Draw("colz");
  bigaus.Draw("same");
  /*gr2->Draw("AF");*/ //gr2->Draw("L");
  /*gr1->Draw("F");*/  //gr1->Draw("L"); 
  c1.Update();
  c1.Write();
  bigaus.Write();

//  delete bigaus;

  /*TLegend leg(0.1,0.6,0.3,0.7); 
  leg.AddEntry(&g1, "68% CL", "l");
  leg.AddEntry(&g2, "90% CL", "l"); 
  leg.AddEntry(&g3, "95% CL", "l"); */
 
  f.Close();
}
}
