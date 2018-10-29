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

  // First to cound number of photons in each hit
  // also N0 set to lower threshold to decrease computation time
  const auto maxSiPM_counts = InitData(sipmToTriggerMap, sipmInfoMap, trigger);
  unsigned N0 = maxSiPM_counts.second;

  // Start reconstruction
  Reconstruct(N0);
  // Make a useful plot
  MakePlot(trigger);
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
  m_recoOutputPath    = config.recoOutputPath;
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

void Analyzer::Reconstruct(unsigned& N0)
{
  // The idea here is to first get a rough estimate of the MLE parameters by
  // a simple grid search (this will be modified based on the application to
  // optimize computation time). Then we will try to apply the Newton-Raphson 
  // method to refine the estimate. 

  // Start the timer for this trigger
  clock_t start = clock();
 
  std::cout << "\nRunning MLE...\n";

  // Start main loop
  // We will cover the entire parameter space, at the expense of computation time
  // Once N0 is calculated, we'll use Newton-Raphson to better approximate x and y
  while (N0 <= 5000) 
  {
    Handle(N0);
    N0 = N0+5;
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
  // Set our initial guesses
  m_oldGuessX      = m_mlX;
  m_oldGuessY      = m_mlY;
  m_oldGuessMLogL  = m_mlLogLikelihood; 
  RefineEstimate(iterator);

  std::cout << "Max log likelihood = " << m_mlLogLikelihood << std::endl
            << "X                  = " << m_mlX             << " cm\n"
            << "Y                  = " << m_mlY             << " cm\n"
            << "Radius             = " << m_mlRadius        << " cm\n"
            << "Theta              = " << m_mlTheta         << " deg\n";
 
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
  //if (weight < 0) { std::cout << "UH OH! WEIGHT < 0!! " << r << " " <<  sipm << "  " <<  cosAngleSiPMToXYandSiPM << "\n"; std::exit(1); }
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
  const float h = 0.0001;

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

  // Multiply inverse Fisher and vector
  auto m = inverseFisherMatrix*lPVector;

  // Our convention --> x, y
  float newGuessX     = m_oldGuessX - m[0][0];
  float newGuessY     = m_oldGuessY - m[1][0];
  float newGuessR(0), newGuessT(0);
  ConvertToPolar(newGuessR, newGuessT, newGuessX, newGuessY);
  
  // Make sure r is not > disk radius
  if (newGuessR >= m_diskRadius) { std::cout << "\n!!!!Caution: Estimate went out of bounds!!!!\n"; return; }
  double newGuessMLogL = ComputeLogLikelihood(newGuessX, newGuessY, m_mlN0);
 
  // We want the eigenvalues and eigenvectors
  // Not sure why we're getting warnings here
  // These are symmetric matrices...
  // Nonetheless, this seems to be returning the correct values
  TVectorD eigenvalues(2);
  TMatrixD eigenvectors = inverseFisherMatrix.EigenVectors(eigenvalues); 

  // Check the eigenvalues and convergence 
  if (//eigenvalues[0] < 0 && eigenvalues[1] < 0 &&
      newGuessMLogL  < m_oldGuessMLogL) 
  { 
    // Update
    //m_mlX             = newGuessX;
    //m_mlY             = newGuessY;
    //ConvertToPolar(m_mlRadius, m_mlTheta, m_mlX, m_mlY);
    //m_mlLogLikelihood = newGuessMLogL;
    
    return; 
  }
  // Otherwise, keep searching
  m_oldGuessX     = newGuessX;
  m_oldGuessY     = newGuessY;
  m_oldGuessMLogL = newGuessMLogL;
  RefineEstimate(iterator);
}

std::pair<unsigned, unsigned> Analyzer::InitData(SiPMToTriggerMap& sipmToTriggerMap, const SiPMInfoMap& sipmInfoMap, const unsigned& trigger)
{ 
  // Which sipm saw the largest number of photons?
  unsigned maxSiPM;
  unsigned max(0);
  // Will will use this as a lower bound for N0
  unsigned total(0);

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
      hit.nPhotons = std::round(hit.hitAmplitude/thisGain);    //( thisGain*( hit.bias - thisBD ) ));   // hit.bias - thisBD
      sipmCounts += hit.nPhotons;
      //std::cout << hit.nPhotons << std::endl;
    }
    if (sipmCounts > max) { max = sipmCounts; maxSiPM = sipm.first; }

    // Store counts in data container
    m_data.emplace(sipm.first, sipmCounts);
    total += sipmCounts;
    std::cout << "SiPM " << sipm.first << " --> " << sipmCounts << " p.e.\n";
  }

  // For use later on
  m_maxCounts = max;

  return std::make_pair(maxSiPM, total);
}

void Analyzer::MakePlot(const unsigned& trigger)
{
  TFile f(m_recoOutputPath.c_str(), "UPDATE");

  // Log Likelihood distribution for m_mlN0
  TH2D logLikelihoodDist("logLikelihoodDist", "Likelihood Profile", 500, -m_diskRadius - 1, m_diskRadius + 1, 500, -m_diskRadius - 1, m_diskRadius + 1); 
  // Make a copy for contours
  auto contour68 = logLikelihoodDist;

  for (unsigned xBin = 1; xBin <= contour68.GetXaxis()->GetNbins(); xBin++)
  {
    for (unsigned yBin = 1; yBin <= contour68.GetYaxis()->GetNbins(); yBin++)
    {
      // Log likelihood for this parameter set
      float x(logLikelihoodDist.GetXaxis()->GetBinCenter(xBin));
      float y(logLikelihoodDist.GetYaxis()->GetBinCenter(yBin));
      double logLikelihood = ComputeLogLikelihood(x, y, m_mlN0);

      float r(0), theta(0);
      ConvertToPolar(r, theta, x, y);

      float diff = -2*(logLikelihood - m_mlLogLikelihood);

      // For plotting purposes
      if (diff > 20.0) diff = 20;  
      if (r > m_diskRadius) diff = 0;
      logLikelihoodDist.SetBinContent(xBin, yBin, diff);
      if (r > (m_diskRadius - 2)) diff = 20;
      contour68.SetBinContent(xBin, yBin, diff);
    }
  }  
  
  // Make copies to draw our contours
  auto contour90 = contour68;
  auto contour95 = contour68;

  // Set the confidence levels
  Double_t level68[1], level90[1], level95[1];
  level68[0] = TMath::ChisquareQuantile(0.68,3); // We had 3 d.o.f
  level90[0] = TMath::ChisquareQuantile(0.90,3);
  level95[0] = TMath::ChisquareQuantile(0.95,3);
  contour68.SetContour(1, level68);
  contour90.SetContour(1, level90);
  contour95.SetContour(1, level95);

  // Now draw the distribution and confidence regions
  std::string name1 = "logLikelihood_CL_" + std::to_string(trigger);
  TCanvas c1(name1.c_str(), name1.c_str(), 800, 800);
  logLikelihoodDist.Draw("colz");

  contour68.SetLineWidth(5);
  contour90.SetLineWidth(5);
  contour95.SetLineWidth(5);
  contour68.SetLineColor(4);
  contour90.SetLineColor(2);
  contour95.SetLineColor(3);
  contour68.Draw("cont3 same");
  contour90.Draw("cont3 same");
  contour95.Draw("cont3 same");

  // Add a marker for the MLE
  TMarker xy(m_mlX, m_mlY, 20);
  xy.SetMarkerSize(2);
  xy.SetMarkerColor(1);
  xy.Draw("same"); 

  // Add our legend
  TLegend leg1(0.1,0.6,0.3,0.7); 
  leg1.AddEntry(&contour68, "68% CL", "l");
  leg1.AddEntry(&contour90, "90% CL", "l"); 
  leg1.AddEntry(&contour95, "95% CL", "l"); 
  leg1.AddEntry(&xy, "MLE Position", "p");
  leg1.Draw("same");
  gStyle->SetPalette(53);
  TColor::InvertPalette();
  c1.Modified();

  c1.Update();
  c1.Write(); 

  // Make the comparison plot 
  std::string name2 = "data_mle_trigger" + std::to_string(trigger);
  TCanvas c2(name2.c_str(), name2.c_str(), 1000, 1000);
  // Get the counts lambda_m using the mlestimates for r, theta, and N0
  std::vector<unsigned> prediction;
  prediction.reserve(m_nSiPMs);
  // Maximum, used for plotting
  float max(m_maxCounts);
  for (int sipm = 1; sipm <= m_nSiPMs; sipm++) 
  {
    float lambda = ComputeLambda(m_mlRadius, m_mlTheta, m_mlN0, sipm); 
    prediction[sipm - 1] = static_cast<unsigned>(lambda);
    if (prediction[sipm - 1] > max) max = prediction[sipm - 1];
  }

  name2 = "pred_trigger" + std::to_string(trigger);
  TH1D pred(name2.c_str(), name2.c_str(), m_nSiPMs, 0, m_nSiPMs); 
  for (int posBin = 1; posBin <= m_nSiPMs; posBin++) pred.SetBinContent( posBin, prediction[posBin - 1] );
  
  pred.SetFillStyle(3001);
  pred.SetFillColor(kRed);
  pred.SetLineWidth(3);
  pred.SetLineColor(kRed);
  pred.GetXaxis()->SetTitle("SiPM Position");
  pred.GetYaxis()->SetTitle("p.e");
  pred.SetTitle("Estimator for SiPM Wheel");
  gStyle->SetOptStat(0);
  pred.SetMaximum(max+1);
  pred.SetMinimum(0);
  pred.Draw();

  name2 = "dataHisto_trigger" + std::to_string(trigger);
  TH1D dataHisto(name2.c_str(), name2.c_str(), m_nSiPMs, 0, m_nSiPMs);
  std::cout << "\nBin comparison:\n";
  for (int posBin = 1; posBin <= m_nSiPMs; posBin++) 
  {
    std::cout << "Data Bin " << posBin << " :  " << m_data.find(posBin)->second << " p.e." << "   Pred Bin " << posBin << " :  " << prediction[posBin - 1] << " p.e." << std::endl;
    dataHisto.SetBinContent(posBin, m_data.find(posBin)->second);
  }
  dataHisto.SetMarkerStyle(21);
  dataHisto.SetMarkerSize(2);
  dataHisto.Draw("same P");

  TLegend leg2(0.1,0.6,0.3,0.7);
  leg2.AddEntry(&pred, "Estimator", "f");
  leg2.AddEntry(&dataHisto, "Data", "p");

  c2.Write();
 
  f.Close();
}
}
