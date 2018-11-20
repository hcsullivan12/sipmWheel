//
// File: Reconstructor.cxx
//
// Author: Hunter Sullivan
//
// Description: Structure to run reconstruction on sipm data.
//

#include "Reconstructor.h"
#include "TMath.h"
#include "TF1.h"
#include "TF2.h"
#include "TH2.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TLegend.h"

namespace wheel {

Reconstructor::Reconstructor()
{}

Reconstructor::~Reconstructor()
{}

void Reconstructor::Reconstruct(SiPMToTriggerMap&    sipmToTriggerMap, 
                                const SiPMInfoMap&   sipmInfoMap, 
                                const Configuration& config, 
                                const unsigned&      trigger)
{
  // Initialize configuration
  std::cout << "\nInitializing reconstruction...\n";  
  Initialize(config); 

  // Count number of photons in each hit
  // also N0 set to lower threshold to decrease computation time
  InitData(sipmToTriggerMap, sipmInfoMap, trigger);
  m_logLikeFunc.Initialize(config, m_data);

  // Start reconstruction
  Reconstruct(m_totalCounts);
  // Make a useful plot
  MakePlot(trigger);
}

void Reconstructor::Initialize(const Configuration& config)
{
  // Initialize
  m_diskRadius        = config.diskRadius;
  if (m_diskRadius <= 0) { std::cout << "Error. Please specify disk radius!\n"; std::exit(1); }
  m_attenuationLength = config.attenuationLength;
  if (m_attenuationLength <= 0) { std::cout << "Error. Please specify an attenuation length!\n"; std::exit(1); }
  if (config.nSiPMs <= 0) { std::cout << "Error. Please specify number of SiPMs!\n"; std::exit(1); }
  m_beta              = 360/config.nSiPMs;
  m_nSiPMs            = config.nSiPMs;
  m_mlLogLikelihood   = std::numeric_limits<double>::lowest();
  m_mlRadius = 0; m_mlTheta = 0; m_mlN0 = 0; 
  m_recoOutputPath    = config.recoOutputPath;
  m_data.clear();
}

void Reconstructor::Reconstruct(unsigned& N0)
{
  // Start the timer for this trigger
  clock_t start = clock();

  if (N0 == 0) { std::cout << "No photons detected!\n"; return; }

  // Initialize minimizer
  ROOT::Math::Minimizer* m_minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
 
  // Set tolerance , etc...
  m_minimizer->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
  m_minimizer->SetMaxIterations(10000);  // for GSL
  m_minimizer->SetTolerance(0.001);
  m_minimizer->SetPrintLevel(2);

  // Create function wrapper for minimizer
  // a IMultiGenFunction type
  ROOT::Math::Functor f(&m_logLikeFunc, &LogLikelihoodFunction::Eval, 3);  
  m_minimizer->SetFunction(f);

  // Make our initial guess and define parameters
  // This will be:
  // r     = R/2
  // theta = angleOfSiPM
  // N0    = 10*m_maxCounts 
  double guessR        = m_diskRadius/2.0;
  double guessThetaDeg = (m_maxSiPM - 1)*m_beta;
  double guessN0       = 10*m_maxCounts; 
  m_minimizer->SetVariable(0, "N0",    guessN0,       0.5);
  m_minimizer->SetVariable(1, "R",     guessR,        0.01);
  m_minimizer->SetVariable(2, "Theta", guessThetaDeg, 0.01);
  m_minimizer->SetVariableLowerLimit(0, N0);
  m_minimizer->SetVariableLimits(1, 0, m_diskRadius);
  m_minimizer->SetVariableLimits(2, 0, 359.9);
 
  std::cout << "\nRunning MLE...\n";
  m_minimizer->Minimize();
 
  // We should have the ml now
  const double *params = m_minimizer->X();
  m_mlN0     = params[0];
  m_mlRadius = params[1];
  m_mlTheta  = params[2];
  m_mlX      = m_mlRadius*TMath::Cos(m_mlTheta*TMath::Pi()/180);
  m_mlY      = m_mlRadius*TMath::Sin(m_mlTheta*TMath::Pi()/180);
  m_mlLogLikelihood = -1*m_minimizer->MinValue();

  std::cout << "Max log likelihood = " << m_mlLogLikelihood << std::endl
            << "X                  = " << m_mlX             << " cm\n"
            << "Y                  = " << m_mlY             << " cm\n"
            << "Radius             = " << m_mlRadius        << " cm\n"
            << "Theta              = " << m_mlTheta         << " deg\n" 
            << "N0                 = " << m_mlN0            << " photons\n";
 
  clock_t end = clock();
  double duration = ((double) (end - start)) / CLOCKS_PER_SEC;
  std::cout << "\nRun time of " << duration << " s" << std::endl;
}

void Reconstructor::ConvertToPolar(float&       r, 
                                   float&       thetaDeg, 
                                   const float& x, 
                                   const float& y)
{
  r        = std::sqrt(x*x + y*y);
  thetaDeg = TMath::ASin(std::abs(y/r))*180/TMath::Pi();

  // Handle theta convention
  if (x < 0 && y > 0) thetaDeg = 180 - thetaDeg;
  if (x < 0 && y < 0) thetaDeg = 180 + thetaDeg;
  if (x > 0 && y < 0) thetaDeg = 360 - thetaDeg;
}

void Reconstructor::InitData(SiPMToTriggerMap&  sipmToTriggerMap, 
                             const SiPMInfoMap& sipmInfoMap, 
                             const unsigned&    trigger)
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
  m_maxCounts   = max;
  m_maxSiPM     = maxSiPM;
  m_totalCounts = total;
}

void Reconstructor::MakePlot(const unsigned& trigger)
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

      float r(0), theta(0);
      ConvertToPolar(r, theta, x, y);

      double logLikelihood = m_logLikeCalc.ComputeLogLikelihood(m_mlN0, r, theta);

      float diff = -2*(logLikelihood - m_mlLogLikelihood);

      // For plotting purposes
      if (diff > 10.0) diff = 10;  
      if (r > m_diskRadius) diff = 0;
      logLikelihoodDist.SetBinContent(xBin, yBin, diff);
      if (r > (m_diskRadius - 0.5)) diff = 10;
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

  TMarker trueXY(0.5, 10, 20);
  trueXY.SetMarkerSize(2);
  trueXY.SetMarkerColor(2);
  trueXY.Draw("same");
 
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
  float max(0);
  for (const auto& d : m_data) if (d.second > max) max = d.second;
  for (int sipm = 1; sipm <= m_nSiPMs; sipm++) 
  {
    float lambda = m_logLikeCalc.ComputeLambda(m_mlN0, m_mlRadius, m_mlTheta, sipm); 
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
