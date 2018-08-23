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

namespace wheel {

Characterizer::Characterizer()
{}

Characterizer::~Characterizer()
{}

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
    // Make the gain plot for each amplitude distribution
    // Reserve the space
    std::vector<TGraphErrors> g;
    g.reserve(ampDists.find(sipm.first)->second.size());
    ampPeaks.emplace(sipm.first, g);
    unsigned biasCounter = 0;
    // Add the gain plot
    for (auto& ampDist : ampDists.find(sipm.first)->second)
    {
      ampPeaks.find(sipm.first)->second.emplace_back(FitGain(ampDist, sipm.first, sipmGains, biasCounter, config));
      biasCounter++;
    }
  }
  // Ignore this for now
  return;

  // Now the breakdown plots
  TCanvas bdPlots("SiPM Breakdowns", "SiPM Breakdowns", 800, 800);
  bdPlots.Divide(4,4); 

  // Loop over sipms gains
  unsigned sipmCounter = 1;
  for (const auto& sipm : sipmGains)
  {
    // Plot versus bias
    // Temp vec for gains and baises
    std::vector<float> biases, gains;
    for (const auto& bias : config.biases) biases.emplace_back(bias);
    std::pair<SiPMGains::iterator, SiPMGains::iterator> lu = sipmGains.equal_range(sipm.first);
    for (auto& g = lu.first; g != lu.second; g++) gains.emplace_back(g->second);

    TGraph g1(config.nBiases, &biases[0], &gains[0]);
    std::string name = "Breakdown of SiPM " + sipm.first;
    g1.SetTitle(name.c_str());
    g1.GetXaxis()->SetTitle("Bias (V)");
    g1.GetYaxis()->SetTitle("Gain (mV/p.e.)");
    g1.GetXaxis()->SetLimits(70,75);
    g1.SetMinimum(0);
    g1.SetMaximum( 1.1*(*std::max_element(gains.begin(), gains.end())) );
    g1.SetMarkerStyle(20);
    g1.SetMarkerSize(2);
    //g1.Draw("AP");

    // Fit for bias plot
    TF1 fit1("fit1", "[0] + [1]*x", *std::min_element(biases.begin(), biases.end()), *std::max_element(biases.begin(), biases.end()));

    fit1.SetParName(1,"Slope");
    fit1.SetParName(0, "YInt");
    g1.Fit(&fit1, "QR");
    gStyle->SetOptFit();

    const float breakdown = -fit1.GetParameter(0)/fit1.GetParameter(1);
    std::cout << "\nSiPM " << sipm.first << " breakdown is... " << breakdown << "\n";

    // Plot versus overvoltage
    std::vector<float> biasOvervoltage(biases.size(), 0);
    unsigned counter = 0;
    for (auto& redBias : biasOvervoltage) {
      redBias = biases[counter] - breakdown;
      counter++;
    }
    TCanvas overvoltage_plots("overvoltage_plot", "Overvoltage Plot", 800, 800);
    TGraph g2(config.nBiases, &biasOvervoltage[0], &gains[0]);
    g2.SetTitle("SiPM Gain");
    g2.GetXaxis()->SetTitle("Over-voltage (#DeltaV)");
    g2.GetYaxis()->SetTitle("Gain (mV/p.e.)");
    g2.GetXaxis()->SetLimits(0,4);
    g2.SetMinimum(0);
    g2.SetMaximum(4);
    g2.SetMarkerStyle(24);
    g2.SetMarkerSize(2);
    //g2.Draw("AP");

    ///Fit for overvoltage plot
    TF1 fit2("fit2", "[0] + [1]*x", 0, 4);
    fit2.SetParName(1,"Slope");
    fit2.SetParName(0, "YInt");
    g2.Fit(&fit2, "QR");
    gStyle->SetOptFit();

    //std::cout << fit2->GetParError(0) << "   " << fit2->GetParError(1) <<  std::endl;
    //for (unsigned i = 0; i < nFiles; i++) std::cout << "Gain " << i + 1 << " = " << gains[i] << std::endl;
    std::cout << "Gain is... " << std::setprecision(3) << fit2.GetParameter(1) << " mV/p.e./O.V.\n";

   /* ///Place the quantities on plot
    std::string bd = std::to_string(breakdown);
    //double breakdown_error = breakdown*std::sqrt( std::pow(fit2->GetParError(0)/fit1->GetParameter(0), 2) + std::pow( fit2->GetParError(1)/fit2->GetParameter(1), 2) );
    //std::string bd_e = std::to_string(breakdown_error);
    std::string text1 = "V_{BD} = " + bd.substr(0,4) + " V";
    TLatex *bd_label = new TLatex(1,0.4,text1.c_str());
    bd_label->SetTextSize(0.025);
    bd_label->Draw("same");

    std::string g = std::to_string(fit2->GetParameter(1));
    std::string e = std::to_string(fit2->GetParError(1));
    std::string text2 = "G = " + g.substr(0,5) + " #pm " + e.substr(0,5) + " mV/p.e./#DeltaV";
    TLatex *gain_label = new TLatex(1,0.35,text2.c_str());
    gain_label->SetTextSize(0.025);
    gain_label->Draw("same");*/

   // Finished this sipm
    bdPlots.cd(sipmCounter);
    sipmCounter++;
  }

  /*TFile f(config.characterizeOutputFile.c_str(), "RECREATE");
  masterAmpDist.Write();
  masterGainPlot.Write();
  f.Close();*/
}

void Characterizer::MakeHistograms(const unsigned& sipm, const std::vector<HitCandidateVec>& vecOfHitCandVec, const Configuration& config)
{
  // Reserve the space
  std::vector<TH1D> distVec;
  distVec.reserve(config.nBiases);
  ampDists.emplace(sipm, distVec);
  
  // ampDists ordered by bias
  for (const auto& bias : config.biases)
  {
    std::string name = "SiPM " + std::to_string(sipm) + " at " + std::to_string(bias) + " V";
    // We will rebin
    TH1D h(name.c_str(), name.c_str(), 1, 0, 1);
    h.GetXaxis()->SetTitle("Integral/a.u.");
    ampDists.find(sipm)->second.emplace_back(h);
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
      if (hit.hitIntegral < xMin[index]) xMin[index] = hit.hitIntegral;
      if (hit.hitIntegral > xMax[index]) xMax[index] = hit.hitIntegral;
    }
  }
  // Now reset the binning
  unsigned index = 0;
  for (auto& dist : ampDists.find(sipm)->second)
  {
    dist.SetBins(100, 0, 100);
    index++;
  }
  // Loop over the hits to fill histos
  for (const auto& hitVec : vecOfHitCandVec)
  {
    for (const auto& hit : hitVec)
    {
      // Get bias for this hit and get the hist that matches this
      const unsigned index = std::distance(config.biases.begin(), config.biases.find(hit.bias));
      auto& dist = ampDists.find(sipm)->second[index];
      dist.Fill(hit.hitIntegral);
    }
  }
}

TGraphErrors Characterizer::FitGain(TH1D& hs, const unsigned& sipm, SiPMGains& sipmGains, const unsigned& nBias, const Configuration& config)
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
  for (int p = 0; p < npeaks; p++) { std::cout << "Found peak at " << peaks[p] << "\n"; peaksVec[p] = peaks[p]; }
  std::sort(peaksVec.begin(), peaksVec.end(), [](const Double_t& left, const Double_t& right){return left < right;});
 
  for (int peak = 0; peak < npeaks; peak++) 
  {
    TF1 gfit("gfit", "gaus", peaksVec[peak] - config.characterizeAmpFitRange/2, peaksVec[peak] + config.characterizeAmpFitRange/2);
    hs.Fit(&gfit,"QR+");    
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

  sipmGains.emplace(sipm, 1000*fit.GetParameter(1) ); //Converted to mV 

  return grpeaks;
}
}
