#include "DNASequence.h"

DNASequence::DNASequence(
  const bam1_t *rec
) : sequence(""), name(""), description("")
{
  const uint8_t *str = bam_get_seq(rec);
  const bam1_core_t * bam_core = &rec->core;

  for (auto i = 0; i < bam_core->l_qseq / 2; ++i) {
    this->sequence += this->lut[str[i] >> 4 & 0x0F];
    this->sequence += this->lut[str[i] &  0x0F];
  }

  if (bam_core->l_qseq % 2) {
    this->sequence += this->lut[str[bam_core->l_qseq / 2] >> 4 & 0x0F];
  }
};

DNASequence::DNASequence(
  const std::string& s,
  const std::string& n,
  const std::string& d,
  const std::string& q
) : sequence(s), name(n), description(d)
{ }

DNASequence::DNASequence(const std::string& s) : DNASequence(s, "", "", "") { }

DNASequence::~DNASequence() = default;

// Getters
const std::string& DNASequence::get_name() const { return this->name; }
const std::string& DNASequence::get_sequence() const { return this->sequence; }
const std::string& DNASequence::get_description() const { return this->description; }

size_t DNASequence::size() const { return this->sequence.size(); }

std::shared_ptr<std::string> DNASequence::reverse_complement() {
  std::shared_ptr<std::string> rev_comp(new std::string(""));

  for (auto it = this->sequence.rbegin(); it != this->sequence.rend(); ++it) {
    switch (*it) {
      case 'A': rev_comp->append("T"); break;
      case 'C': rev_comp->append("G"); break;
      case 'G': rev_comp->append("C"); break;
      case 'T': rev_comp->append("A"); break;
      case 'a': rev_comp->append("t"); break;
      case 'c': rev_comp->append("g"); break;
      case 'g': rev_comp->append("c"); break;
      case 't': rev_comp->append("a"); break;
      case 'N': rev_comp->append("N"); break;
      case 'n': rev_comp->append("n"); break;
      default : rev_comp->append("N"); break;
    }
  }

  return rev_comp;
}

const char& DNASequence::at(size_t i) const { return this->sequence.at(i); }

// Estimation
bool DNASequence::is_valid() const { return this->sequence != ""; }
bool DNASequence::is_valid_deep() const {
  for (const auto c : this->sequence) {
    switch(c) {
      case 'A': case 'C': case 'G': case 'T': case 'N': break;
      case 'a': case 'c': case 'g': case 't': case 'n': break;
      default:
        return false;
    }
  }
  return true;
}

const std::string& DNASequence::to_string() const { return this->sequence; }

std::ostream& operator<<(std::ostream& o, const DNASequence& s) { o << s.sequence; return o; }
std::ostream& operator<<(std::ostream& o, const std::shared_ptr<DNASequence> s) { o << s->sequence; return o; }
