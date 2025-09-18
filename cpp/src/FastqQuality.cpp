#include "FastqQuality.h"

FastqQuality::FastqQuality(const bam1_t *rec, int base) : quality(""), description("") {
  const uint8_t* q_str = bam_get_qual(rec);
  const auto *bam_core = &rec->core;

  if (q_str[0] != 0xff) {
    this->quality = std::string((char*)q_str, bam_core->l_qseq);
    std::transform(this->quality.begin(), this->quality.end(), this->quality.begin(), [base, this](char c) -> char {
      if (c < 33 || c > 126) { this->valid = false; }
      return c + base;
    });
  }
}

FastqQuality::FastqQuality(const bam1_t *rec): FastqQuality(rec, 33) { }
FastqQuality::FastqQuality(const std::string& s): quality(s) { }
FastqQuality::~FastqQuality() = default;


bool FastqQuality::is_valid() const { return this->quality != ""; }
bool FastqQuality::is_valid_deep() const {
  for (const auto x : this->quality) {
    if (x < 33 || x > 126) { return false; }
  }
  return true;
}

const std::string& FastqQuality::to_string() const { return this->quality; }
const std::string& FastqQuality::get_quality() const { return this->quality; }

// Stream
std::ostream& operator<<(std::ostream& o, const FastqQuality& fq) { o << fq.quality; return o; }
std::ostream& operator<<(std::ostream& o, const std::shared_ptr<FastqQuality> fq) { o << fq->quality; return o; }
