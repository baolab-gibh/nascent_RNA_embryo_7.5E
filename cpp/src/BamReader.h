#pragma once
#include <memory>
#include <string>
#include <iostream>
#include <algorithm>

#include "Typedefs.h"
#include "SamRecord.h"
#include "SamAuxTag.h"
#include "CigarString.h"
#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include "htslib/hts_endian.h"

class BamReader {
  size_t count{0};

  std::string filename;
  std::string version;

  samFile *fp;
  sam_hdr_t *header;
  bam1_t *bam_rec;

  const char* lut{"=ACMGRSVTWYHKDBN"};

public:
  BamReader() = delete;
  BamReader(const std::string& filename);
  ~BamReader();

  std::shared_ptr<SamRecord> next();
  std::shared_ptr<SamAuxTag<std::string>> get_tag(const uint8_t*, const uint8_t* end);
  size_t current_line();

  typedef std::shared_ptr<BamReader> BamReaderSPtr;
};
