#include "CigarString.h"

CigarString::CigarString(const std::string &str) { }
CigarString::CigarString(const std::vector<std::pair<size_t, CigarOperation>> &op) : op(std::move(op)) {}
CigarString::CigarString(const bam1_t *rec) {
  size_t ops_len{0};
  CigarOperation ops{CigarOperation::MATCH};
  const bam1_core_t *bam_core = &rec->core;

  const uint32_t *_cigar = bam_get_cigar(rec);
  for (uint32_t i = 0; i < bam_core->n_cigar; ++i) {
    ops_len = bam_cigar_oplen(_cigar[i]);
    ops = static_cast<CigarOperation>(bam_cigar_opchr(_cigar[i]));
    this->op.push_back(std::make_pair(ops_len, ops));
  }
}

CigarString::~CigarString() {}

std::shared_ptr<std::string> CigarString::to_string() const {
  std::shared_ptr<std::string> str(new std::string(""));
  for (auto x : this->op) { *str += std::to_string(x.first) + char(x.second); }
  return str;
}

std::vector<std::pair<size_t, CigarOperation>> CigarString::get_op() const {
  return this->op;
}

size_t CigarString::infer_sequence_length() const {
  for (auto x : this->op) {
    switch (x.first) {
      case CigarOperation::MATCH:
      case CigarOperation::SEQUENCE_MATCH:
      case CigarOperation::SEQUENCE_MISMATCH:
        return x.second;
      case CigarOperation::INSERTION:
      case CigarOperation::DELETION:
      case CigarOperation::REFERENCE_SKIP:
      case CigarOperation::SOFT_CLIPPING:
      case CigarOperation::HARD_CLIPPING:
      case CigarOperation::PADDING:
        break;
      default:
        break;
    }
  }

  return 0;
}


// Streams
std::ostream& operator<<(std::ostream& o, const CigarString cs)     { o << *cs.to_string(); return o; }
std::ostream& operator<<(std::ostream& o, const CigarStringSPtr cs) { return o << *cs;   }
