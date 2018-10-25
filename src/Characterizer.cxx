// 
// File: Characterizer.cxx
//
// Author: Hunter Sullivan
//
// Discription: Structure to characterize sipms.
//

#include "Characterizer.h"
#include "TSpectrum.h"
#include "TFile.h"
#include "TLatex.h"

namespace wheel {

Characterizer::Characterizer()
{}

Characterizer::~Characterizer()
{}

void Characterizer::Initialize(const Configuration& config)
{
  // Initialize our gain container 
  for (unsigned i = 1; i <= config.nSiPMs; i++)
  {    
    // This is the map from sipm to bias voltage to gain[0] and error[1]
    std::vector<std::vector<float>> b(config.nBiases,std::vector<float>(2,0));
    m_sipmGains.push_back(b);
  }
}

void Characterizer::Characterize(SiPMInfoMap& sipmInfoMap, const SiPMToTriggerMap& sipmToTriggerMap, const Configuration& config)
{    
  // Loop over the sipms and biases
  unsigned ampCounter(1);
  unsigned gainCounter(1);
  SiPMGains sipmGains;
  // There should only be one sipm here 
  for (const auto& sipm : sipmToTriggerMap)
  {
    // Make the amplitude distributions
    MakeHistograms(sipm.first, sipm.second, config); 
    std::cout << std::endl;
  }
  return;
}

void Characterizer::MakeHistograms(const unsigned& sipm, const std::vector<HitCandidateVec>& vecOfHitCandVec, const Configuration& config)
{
  // Reserve the space
  std::vector<TH1D> distVec;
  distVec.reserve(config.nBiases);
  m_ampDists.emplace(sipm, distVec);
  
  // m_ampDists ordered by bias
  for (const auto& bias : config.biases)
  {
    std::string name = "SiPM " + std::to_string(sipm) + " at " + std::to_string(bias) + " V";
    // We will rebin
    TH1D h(name.c_str(), name.c_str(), 1, 0, 1);
    h.GetXaxis()->SetTitle("Amplitude/V");
    m_ampDists.find(sipm)->second.emplace_back(h);
  }
  // Create a container for min and max
  std::vector<float> xMin(config.nBiases, std::numeric_limits<float>::max());
  std::vector<float> xMax(config.nBiases, 0);
  // Need to loop over hits initially to get min and max for binning
  for (const auto& hitVec : vecOfHitCandVec)
  {
    for (const auto& hit : hitVec)
    {
      // Get bias index for this hit 
      const unsigned index = std::distance(config.biases.begin(), config.biases.find(hit.bias));
      if (hit.hitAmplitude < xMin[index]) xMin[index] = hit.hitAmplitude;
      if (hit.hitAmplitude > xMax[index]) xMax[index] = hit.hitAmplitude;
    }
  }
  // Now reset the binning
  unsigned index = 0;
  for (auto& dist : m_ampDists.find(sipm)->second)
  {
    dist.SetBins(50, 0, 0.01); //xMax[0]);
    index++;
  }
  // Loop over the hits to fill histos
  for (const auto& hitVec : vecOfHitCandVec)
  {
    for (const auto& hit : hitVec)
    {
      // Get bias for this hit and get the hist that matches this
      const unsigned index = std::distance(config.biases.begin(), config.biases.find(hit.bias));
      auto& dist = m_ampDists.find(sipm)->second[index];
      dist.Fill(hit.hitAmplitude);
    }
  }
}

TGraphErrors Characterizer::FitGain(TH1D& hs, const unsigned& sipm, const unsigned& nBias, const Configuration& config)
{
  Float_t gain=0;
  TSpectrum s(3);

  Int_t nfound = s.Search(&hs, config.characterizeAmpSig, "", config.characterizeAmpThr);
  Int_t npeaks = s.GetNPeaks();
  printf("Found %d peaks to fit\n",npeaks);

  // Get the peaks 
  auto peaks = s.GetPositionX();

  Double_t y[npeaks],  ey[npeaks];
  Double_t gx[npeaks], gex[npeaks];
  Double_t gy[npeaks], gey[npeaks];
 
  // Print out the peaks found and store into a vec
  std::vector<Double_t> peaksVec(npeaks, 0);
  for (int p = 0; p < npeaks; p++) { /*std::cout << "Found peak at " << peaks[p] << "\n"; */peaksVec[p] = peaks[p]; }
  std::sort(peaksVec.begin(), peaksVec.end(), [](const Double_t& left, const Double_t& right){return left < right;});
 
  for (int peak = 0; peak < npeaks; peak++) 
  {
    TF1 gfit("gfit", "gaus", peaksVec[peak] - config.characterizeAmpFitRange/2, peaksVec[peak] + config.characterizeAmpFitRange/2);
    hs.Fit(&gfit,"QR+");    
    //std::cout << "before = " << peaksVec[peak] << "   after = " << gfit.GetParameter(1) << std::endl;
    gx[peak]  = peak + 1;
    gy[peak]  = gfit.GetParameter(1);  //Mean
    gey[peak] = gfit.GetParError(1);
  }

  TGraphErrors grpeaks(npeaks,gx,gy,0,gey);
  std::string name = "SiPM " + std::to_string(sipm) + " Gain from " + std::to_string(*std::next(config.biases.begin(), nBias));
  grpeaks.SetTitle(name.c_str());
  grpeaks.GetYaxis()->SetTitle("Area/a.u.");
  
  grpeaks.GetYaxis()->SetTitleOffset(1.4);
  grpeaks.GetXaxis()->SetTitle("Peak N");
  grpeaks.SetMarkerColor(4);
  grpeaks.SetMarkerStyle(20);
  grpeaks.SetFillColor(0);

  // Fit
  TF1 fit("fit","[0] + [1]*x",0.5,npeaks+1);
  fit.SetParName(1,"Gain");
  fit.SetParName(0, "Pedestal");
  grpeaks.Fit(&fit, "QR");
  gStyle->SetOptFit();

  // Save this gain
  std::vector<float> info;
  m_sipmGains[sipm-1][nBias][0] = fit.GetParameter(1);
  m_sipmGains[sipm-1][nBias][1] = fit.GetParError(1);

  return grpeaks;
}

void Characterizer::SaveCharacterizationPlots(const wheel::Configuration& config)
{
  // First make the gain plots
  std::cout << "Finding peaks for gain plots...\n";
  unsigned biasCounter = 0;
  // Containers for our graphs
  TGraphErrors* gainGraphs[config.nSiPMs][config.nBiases];
  TGraphErrors* sipmGainGraphs[config.nSiPMs];
  TGraphErrors* overvoltageGraphs[config.nSiPMs];
  TLatex*       bdLabels[config.nSiPMs];
  TLatex*       gainLabels[config.nSiPMs];

  // Useful containers for gains and biases
  float biases[config.nBiases], gy[config.nBiases];
  for (unsigned b = 0; b < config.nBiases; b++)
  {
    biases[b] = *std::next(config.biases.begin(), b); 
    // fix me!
    gy[b]     = 0.001;
  }

  // Add the gain plots
  for (unsigned sipm = 1; sipm <= config.nSiPMs; sipm++)
  {
    unsigned biasCounter(0);
    for (auto& ampDist : m_ampDists.find(sipm)->second)
    {
      auto plot = FitGain(ampDist, sipm, biasCounter, config);
      auto x    = plot.GetX();
      auto y    = plot.GetY();
      auto ex   = plot.GetEX();
      auto ey   = plot.GetEY();
      gainGraphs[sipm-1][biasCounter] = new TGraphErrors(plot.GetN(), x, y, ex, ey);

      biasCounter++;
    }

    float gains[config.nBiases], ge[config.nBiases];
    for (unsigned g = 0; g < config.nBiases; g++) 
    {
      gains[g] = m_sipmGains[sipm-1][g][0];
      ge[g]    = m_sipmGains[sipm-1][g][1];
    }
    sipmGainGraphs[sipm-1]    = new TGraphErrors(config.nBiases, biases, gains, 0, ge);  
    // Will change later
    overvoltageGraphs[sipm-1] = new TGraphErrors(config.nBiases, biases, gains, 0, ge);
  }

  std::cout << "\nMaking and saving plots...\n";

  // Canvases for out histos
  // Amp dist
  TCanvas masterAmpDist("masterAmpDist",  "All Amplitude Distributions", 1000, 1000);
  // Peak versus n p.e.
  TCanvas masterGainPlot("masterGainPlot", "All Gains", 1000, 1000); 
  // Gain versus bias
  TCanvas bdPlots("SiPM Breakdowns", "SiPM Breakdowns", 800, 800);
  // Gain versus overvoltage
  TCanvas overvoltagePlots("overvoltagPlots", "Overvoltage Plots", 800, 800);

  // Divide the canvases
  masterAmpDist.Divide(config.nBiases,config.nSiPMs);
  masterGainPlot.Divide(config.nBiases,config.nSiPMs);
  bdPlots.Divide(4,2); 
  overvoltagePlots.Divide(4,2); 

  // Loop over sipms
  unsigned ampCounter(1);
  unsigned gainCounter(1);
  for (unsigned sipm = 1; sipm <= config.nSiPMs; sipm++)
  {
    // Draw the amp dists 
    for (auto& dist : m_ampDists.find(sipm)->second)
    {
      // Amplitude dist
      masterAmpDist.cd(ampCounter);
      gStyle->SetOptStat(0);
      dist.GetXaxis()->SetTitle("Integral/a.u.");
      dist.Draw();
      masterAmpDist.Update();
      dist.GetXaxis()->SetTitle("Area/a.u."); 
      masterAmpDist.Modified();
      ampCounter++;
    }
    // Draw the gain plot for each sipm and bias
    for (unsigned bias = 0; bias < config.nBiases; bias++)
    {
      // Since we had to make copies, we need to refit
      std::string name = "SiPM " + std::to_string(sipm) + " Gain from " + std::to_string(*std::next(config.biases.begin(), bias));
      gainGraphs[sipm-1][bias]->SetTitle(name.c_str());
      gainGraphs[sipm-1][bias]->GetYaxis()->SetTitle("Area/a.u.");
      gainGraphs[sipm-1][bias]->GetYaxis()->SetTitleOffset(1.4);
      gainGraphs[sipm-1][bias]->GetXaxis()->SetTitle("Peak N");
      gainGraphs[sipm-1][bias]->SetMarkerColor(4);
      gainGraphs[sipm-1][bias]->SetMarkerStyle(20);
      gainGraphs[sipm-1][bias]->SetFillColor(0);

      // Fit
      TF1 fit("fit","[0] + [1]*x", 0.5, gainGraphs[sipm-1][bias]->GetN()+1);
      fit.SetParName(1,"Gain");
      fit.SetParName(0, "Pedestal");
      gainGraphs[sipm-1][bias]->Fit(&fit, "QR");
      gStyle->SetOptFit();

      masterGainPlot.cd(gainCounter);
      gainGraphs[sipm-1][bias]->Draw("AP");
      gainCounter++;
      gainGraphs[sipm-1][bias]->GetXaxis()->SetLimits(0, 5);
      gainGraphs[sipm-1][bias]->SetMinimum(0);      
      gainGraphs[sipm-1][bias]->SetMaximum(0.01);
    }
 
    // Now draw the breakdown plot and overvoltage plot
    bdPlots.cd(sipm);
//    overvoltagePlots.cd(sipm);

    // Plot versus bias
    sipmGainGraphs[sipm-1]->SetFillColor(0);
    std::string name = "Breakdown of SiPM " + sipm;
    sipmGainGraphs[sipm-1]->SetTitle(name.c_str());
    sipmGainGraphs[sipm-1]->GetXaxis()->SetTitle("Bias (V)");
    sipmGainGraphs[sipm-1]->GetYaxis()->SetTitle("Gain (V/p.e.)");
    sipmGainGraphs[sipm-1]->GetXaxis()->SetLimits(70,75);
    sipmGainGraphs[sipm-1]->SetMinimum(0);
    sipmGainGraphs[sipm-1]->SetMaximum( 0.01 );
    sipmGainGraphs[sipm-1]->SetMarkerStyle(20);
    sipmGainGraphs[sipm-1]->SetMarkerSize(2);
    sipmGainGraphs[sipm-1]->Draw("AP");

    // Fit for sipm gain plot to get the breakdown
    TF1 fit1("fit1", "[0] + [1]*x", *std::min_element(config.biases.begin(), config.biases.end()), *std::max_element(config.biases.begin(), config.biases.end()));

    fit1.SetParName(1,"Slope");
    fit1.SetParName(0, "YInt");
    sipmGainGraphs[sipm-1]->Fit(&fit1, "QR");
    gStyle->SetOptFit();

    const float breakdown = -fit1.GetParameter(0)/fit1.GetParameter(1);
    std::cout << "\nSiPM " << sipm << " breakdown is... " << breakdown << "\n";

    // Plot versus overvoltage
    float biasOvervoltage[config.nBiases];
    unsigned counter = 0;
    for (auto& overvoltage : biasOvervoltage) 
    {
      overvoltage = biases[counter] - breakdown;
      // Just so we can see if something went wrong
      if (overvoltage < 0) overvoltage = biases[counter] - 70.0;
      counter++;
    }

    // Now plot versus overvoltage
    overvoltagePlots.cd(sipm);

    // We need to change the x coordinates
    auto y = overvoltageGraphs[sipm-1]->GetY();
    for (unsigned p = 0; p < overvoltageGraphs[sipm-1]->GetN(); p++) overvoltageGraphs[sipm-1]->SetPoint(p, biasOvervoltage[p], y[p]);
    overvoltageGraphs[sipm-1]->SetTitle("SiPM Gain");
    overvoltageGraphs[sipm-1]->GetXaxis()->SetTitle("Over-voltage (#DeltaV)");
    overvoltageGraphs[sipm-1]->GetYaxis()->SetTitle("Gain (mV/p.e.)");
    overvoltageGraphs[sipm-1]->GetYaxis()->SetTitleOffset(1.4);
    overvoltageGraphs[sipm-1]->GetXaxis()->SetLimits(0,4);
    overvoltageGraphs[sipm-1]->SetMinimum(0);
    overvoltageGraphs[sipm-1]->SetMaximum(0.01);
    overvoltageGraphs[sipm-1]->SetMarkerStyle(24);
    overvoltageGraphs[sipm-1]->SetMarkerSize(2);
    overvoltageGraphs[sipm-1]->Draw("AP");

    ///Fit for overvoltage plot
    TF1 fit2("fit2", "[0] + [1]*x", 0, 4);
    fit2.SetParName(1,"Slope");
    fit2.SetParName(0, "YInt");
    overvoltageGraphs[sipm-1]->Fit(&fit2, "QR");
    gStyle->SetOptFit();

    std::cout << "Gain is... " << std::setprecision(3) << fit2.GetParameter(1) << " mV/p.e./O.V.\n";
    ///Place the quantities on plot
    std::string bd = std::to_string(breakdown);
    //double breakdown_error = breakdown*std::sqrt( std::pow(fit2->GetParError(0)/fit1->GetParameter(0), 2) + std::pow( fit2->GetParError(1)/fit2->GetParameter(1), 2) );
    //std::string bd_e = std::to_string(breakdown_error);
    std::string text1 = "V_{BD} = " + bd.substr(0,4) + " V";
    bdLabels[sipm] = new TLatex(1,0.4,text1.c_str());
    bdLabels[sipm]->SetTextSize(0.025);
    bdLabels[sipm]->Draw("same");

    std::string g = std::to_string(fit2.GetParameter(1));
    std::string e = std::to_string(fit2.GetParError(1));
    std::string text2 = "G = " + g.substr(0,5) + " #pm " + e.substr(0,5) + " mV/p.e./#DeltaV";
    gainLabels[sipm] = new TLatex(1,0.35,text2.c_str());
    gainLabels[sipm]->SetTextSize(0.025);
    gainLabels[sipm]->Draw("same");
  }

  // Write the plots to our output file
  TFile f(config.characterizeOutputPath.c_str(), "RECREATE");
  masterAmpDist.Write();
  masterGainPlot.Write();
  bdPlots.Write();
  overvoltagePlots.Write();
  f.Close();
}
}
