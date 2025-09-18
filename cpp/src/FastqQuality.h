#pragma once
#include <iostream>
#include <string>
#include <memory>
#include <algorithm>
#include "htslib/sam.h"

class FastqQuality {
  std::string quality;
  std::string description;
  bool valid;

public:
  FastqQuality() = delete;
  FastqQuality(const bam1_t *rec);
  FastqQuality(const bam1_t *rec, int base);
  FastqQuality(const std::string&);
  ~FastqQuality();

  // Getters
  const std::string& get_quality() const;

  // Estimations
  bool is_valid() const;
  bool is_valid_deep() const;

  // Operations
  const std::string& to_string() const;

  // Stream
  friend std::ostream& operator<<(std::ostream& o, const FastqQuality& fq);
  friend std::ostream& operator<<(std::ostream& o, const std::shared_ptr<FastqQuality> fq);
};

typedef std::shared_ptr<FastqQuality> FastqQualitySPtr;
