// 
// File: WaveformAlg.cxx
//
// Author: Hunter Sullivan
//
// Discription: Structure to perform various algos on the input waveform.
//

#include "WaveformAlg.h"
#include <algorithm>
#include <cmath>
#include "TH1.h"
#include "TF1.h"

namespace wheel {

WaveformAlg::WaveformAlg()
{}

WaveformAlg::~WaveformAlg()
{}

void WaveformAlg::SmoothWaveform(std::vector<float>& signal, const Configuration& config)
{
  // Going to perform a moving average here
  unsigned sample(1);
  for (auto& amp : signal)
  {
    // Ignore first few points
    if (sample <= config.smaRange) 
    {
      sample++;
      continue;
    }
    // Ignore last points
    if (sample >= signal.size() - config.smaRange) return;

    float leadSum(0), trailSum(0);
    for (unsigned offset = 1; offset <= config.smaRange; offset++)
    {
      leadSum  += signal[sample - offset];
      trailSum += signal[sample + offset];
    }
    float rangeAvg = (leadSum + trailSum + amp)/(2*config.smaRange + 1);

    // Reassign
    amp = rangeAvg;
    sample++;
  }
}

void WaveformAlg::SmoothWaveform2(std::vector<float>& signal, const Configuration& config)
{
  // Going to perform a moving average here
  unsigned sample(1);
  for (auto& amp : signal)
  {
    // Ignore first few points
    if (sample <= config.smaRange) 
    {
      sample++;
      continue;
    }
    // Ignore last points
    if (sample >= signal.size() - config.smaRange) return;

    float leadSum(0);
    for (unsigned offset = 1; offset <= config.smaRange; offset++)
    {
      leadSum += signal[sample - offset];
    }
    float rangeAvg = (leadSum + amp)/(config.smaRange + 1);

    // Reassign
    amp = rangeAvg;
    sample++;
  }
}

void WaveformAlg::FindHitCandidates(std::vector<float>&  waveform,
                                    size_t               roiStartTick,
                                    size_t               channel,
                                    const float&         bias,
		                                HitCandidateVec&     hitCandVec,
                                    const Configuration& config)
{
  // Find noise parameters first
  float baselineInter(0);
  const auto noiseParameters = ComputeNoise(waveform, baselineInter, config);

  // Use the recursive version to find the candidate hits
  // We only need to search a limited range 
  FindHitCandidates(waveform.begin(), waveform.end(), noiseParameters, bias, roiStartTick, hitCandVec, config);  
  return;
}

void WaveformAlg::FindHitCandidates(std::vector<float>::const_iterator startItr,
                                    std::vector<float>::const_iterator stopItr,
                                    const std::pair<float, float>&     noiseParameters,
                                    const float&                       bias,
                                    size_t                             roiStartTick,
		                                HitCandidateVec&                   hitCandVec,
                                    const Configuration&               config)
{
  unsigned configOffset = config.hitFinderSearch;

  // Need a minimum number of ticks to do any work here
  if (std::distance(startItr,stopItr) > configOffset)
  {
    // Find the highest peak in the range given
    auto maxItr = std::max_element(startItr, stopItr);

    float maxValue = *maxItr;
    int   maxTime  = std::distance(startItr,maxItr);

    float hitThreshold;
    float hitThreshold = std::max( noiseParameters.first + config.hitSigma*noiseParameters.second, config.minimumHitAmp ); 

    //std::cout << noiseParameters.first + hitSigma*noiseParameters.second << "  " << minimumHitAmp << "  " << maxValue << std::endl;

    if (maxValue > hitThreshold)
    {
      // backwards to find first bin for this candidate hit
      auto firstItr = std::distance(startItr,maxItr) > configOffset ? maxItr - (configOffset - 1) : startItr;

      while(firstItr != startItr)
      {
        // Check for pathology where waveform goes too negative
        if (*firstItr < -hitThreshold) break;

        // Check both sides of firstItr and look for min/inflection point
        if (*firstItr < *(firstItr+configOffset) && *firstItr <= *(firstItr-configOffset)) break;

        firstItr-=configOffset;
      }

      float min     = *firstItr;
      int firstTime = std::distance(startItr,firstItr);

      // Recursive call to find all candidate hits earlier than this peak
      FindHitCandidates(startItr, firstItr + 1, noiseParameters, bias, roiStartTick, hitCandVec, config);

      // forwards to find last bin for this candidate hit
      auto lastItr = std::distance(maxItr,stopItr) > configOffset ? maxItr + (configOffset - 1) : stopItr - 1;

      while(lastItr != stopItr - 1)
      {
        // Check for pathology where waveform goes too negative
        if (*lastItr < -hitThreshold) break;

        // Check both sides of firstItr and look for min/inflection point
        if (*lastItr <= *(lastItr+configOffset) && *lastItr < *(lastItr-configOffset)) break;

        lastItr+=configOffset;
      }

      int lastTime = std::distance(startItr,lastItr);

      // Now save this candidate's start and max time info
      HitCandidate hitCandidate;
      hitCandidate.startTick     = roiStartTick + firstTime;
      hitCandidate.stopTick      = roiStartTick + lastTime;
      hitCandidate.hitBase       = min;
      hitCandidate.hitCenter     = roiStartTick + maxTime;
      hitCandidate.hitHeight     = maxValue;
      hitCandidate.bias          = bias;

      hitCandVec.push_back(hitCandidate);

      // Recursive call to find all candidate hits later than this peak
      FindHitCandidates(lastItr + 1, stopItr, noiseParameters, bias, roiStartTick + std::distance(startItr,lastItr + 1), hitCandVec, config);
    }
  }

  return;
}

std::pair<float, float> WaveformAlg::ComputeNoise(std::vector<float>& signal, float& baselineInter, const Configuration& config)
{
  // Find min and max of this signal for number of bins
  auto min_max = std::minmax_element( signal.begin(), signal.end() );
  int minSample = static_cast<int>( min_max.first - signal.begin() );
  int maxSample = static_cast<int>( min_max.second - signal.begin() );

  float min = signal.at(minSample);
  float max = signal.at(maxSample);

  unsigned nBins = static_cast<unsigned>(std::abs(max - min)*10000);

  // Make a histogram to compute mean noise
  TH1S amplitudeDist("amplitudeDist", "amplitudeDist", nBins, min, max);
  for (unsigned sample = 0; sample < signal.size(); ++sample) {
    amplitudeDist.Fill( signal.at(sample) );
  }

  // Fit the amplitude distribution
  TF1 gauss("gauss", "gaus");
  amplitudeDist.GetXaxis()->SetRangeUser(amplitudeDist.GetMean() - 3 * amplitudeDist.GetStdDev(), amplitudeDist.GetMean() + 3 * amplitudeDist.GetStdDev());
  amplitudeDist.Fit(&gauss, "QN"); 
  // Extrapolate
  TF1 line("line", "[0] + [1]*(x-[2])", 0, max);
  float x     = gauss.GetParameter(1) + 0.4*gauss.GetParameter(2);
  float y     = gauss.Eval(x);
  float slope = gauss.Derivative(x);
  line.SetParameters(y,slope,x);
  TF1 fint("fint", "TMath::Abs(0-line)", gauss.GetParameter(1), max);

  baselineInter = fint.GetMinimumX();
  return std::pair<double, double>(gauss.GetParameter(1), gauss.GetParameter(2));
}

}
