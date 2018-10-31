// 
// File: ConfigReader.h
//
// Author: Hunter Sullivan
//
// Discription: Structure to read the configuration file.
//

#ifndef CONFIGREADER_H
#define CONFIGREADER_H

#include "Utilities.h"
#include <fstream>
#include <iostream>

namespace wheel {

class ConfigReader {

public:
  ConfigReader();
  ~ConfigReader();
  
  void ReadFile(Configuration& config);
  void OutputConfigInfo(const Configuration& config);

private:

  void RecordBiases(Configuration& config, const std::string& value);
  void RecordGains(std::map<unsigned, float>& map, const Configuration& config, const std::string& value);
  void RecordPosition(Configuration& config, const std::string& value);
 
};
}

#endif
