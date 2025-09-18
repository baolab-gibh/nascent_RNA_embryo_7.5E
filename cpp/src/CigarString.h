#pragma once
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include "htslib/sam.h"
#include "htslib/kstring.h"

enum CigarOperation {
  MATCH             = 'M', // M
  INSERTION         = 'I', // I
  DELETION          = 'D', // D
  REFERENCE_SKIP    = 'N', // N
  SOFT_CLIPPING     = 'S', // S
  HARD_CLIPPING     = 'H', // H
  PADDING           = 'P', // P
  SEQUENCE_MATCH    = '=', // =
  SEQUENCE_MISMATCH = 'X'  // X
};

class CigarString {
  std::vector<std::pair<size_t, CigarOperation>> op;

public:
  CigarString();
  CigarString(const bam1_t *rec);
  CigarString(const std::string&);
  CigarString(const std::vector<std::pair<size_t, CigarOperation>> &);
  CigarString(const bam1_t*, const bam1_core_t*);
  ~CigarString();


  // Getters
  std::vector<std::pair<size_t, CigarOperation>> get_op() const;
  std::shared_ptr<std::string> to_string() const;

  // Inferred results
  size_t infer_sequence_length() const;

  // Stream
  friend std::ostream& operator<<(std::ostream& o, const std::shared_ptr<CigarString> cs); 
  friend std::ostream& operator<<(std::ostream& o, const CigarString cs); 

private:
  void _parse_cigar_string(const std::string &);
};

typedef std::shared_ptr<CigarString> CigarStringSPtr;
