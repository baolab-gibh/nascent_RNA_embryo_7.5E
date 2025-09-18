#pragma once
#include <string>
#include "FastqQuality.h"
#include "DNASequence.h"

class FastqRecord {
  std::string name;
  DNASequence sequence;
  FastqQuality quality;

public:
  FastqRecord(const std::string& n, const std::string& s, const std::string& q, const std::string& d);
  FastqRecord(const std::string& n, const std::string& s, const std::string& q);
};
