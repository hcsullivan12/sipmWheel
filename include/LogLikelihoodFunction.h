// 
// File: LogLikelihoodFunction.h
//
// Author: Hunter Sullivan
//
// Description: Structure for log likelihood function. 
//

#ifndef LOGLIKELIHOODFUNCTION_H
#define LOGLIKELIHOODFUNCTION_H

#include "LogLikelihoodCalculator.h"

namespace wheel {

class LogLikelihoodFunction 
{
public:
  double operator() (const double* par) 
  {
    return -1*m_logLikeCalc.ComputeLogLikelihood(par[0], par[1], par[3]); 
  }

  void Initialize(const Configuration&                config, 
                  const std::map<unsigned, unsigned>& data)
  {
    m_logLikeCalc.Initialize(config, data);
  }

  double Eval(const double* par)  
  {
    return -1*m_logLikeCalc.ComputeLogLikelihood(par[0], par[1], par[2]); 
  }

private:
  LogLikelihoodCalculator m_logLikeCalc;
};

}

#endif
