#pragma once
#include <string>
#include <memory>
#include "htslib/sam.h"

class SamHeader {
  std::string filename;
  samFile *fp;
  sam_hdr_t *header;

public:
  SamHeader() = delete;
  SamHeader(const std::string& fn);
  SamHeader(std::string fn, sam_hdr_t *header);
  ~SamHeader();
};
