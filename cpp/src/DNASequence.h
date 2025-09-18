#pragma once
#include <string>
#include <memory>
#include "htslib/sam.h"

enum STRAND { FORWARD = 0, REVERSE = 1 };
enum VALID_BASES { A = 0, C = 1, G = 2, T = 3, N = 4, a = 5, c = 6, g = 7, t = 8, n = 9 };

class DNASequence {
  std::string sequence;
  std::string name;
  std::string description;
  STRAND strand;

  const char* lut{"=ACMGRSVTWYHKDBN"};

public:
  DNASequence() = delete;
  DNASequence(const bam1_t *rec);
  DNASequence(const std::string&);
  DNASequence(const std::string&, const std::string&, const std::string&, const std::string&);

  ~DNASequence();

  // Getters
  const std::string& get_name() const;
  const std::string& get_sequence() const;
  const std::string& get_description() const;

  size_t size() const;

  // Estimation
  bool is_valid() const;
  bool is_valid_deep() const;

  // Operations
  std::shared_ptr<std::string> reverse_complement();
  void substract();
  const std::string& to_string() const;
  const char& at(size_t) const;

  // Stream
  friend std::ostream& operator<<(std::ostream& o, const DNASequence& s);
  friend std::ostream& operator<<(std::ostream& o, const std::shared_ptr<DNASequence> s);

  // Match
  bool match(const std::string& s) const;
  bool fuzzy_match(const std::string& s) const;

  // Mutation operations
  void addMutations();
};

typedef std::shared_ptr<DNASequence> DNASequenceSPtr;
