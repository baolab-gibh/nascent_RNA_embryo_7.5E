#pragma once
#include <iostream>
#include <string>
#include <memory>
#include <stack>

#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include "htslib/hts_endian.h"

#include "Typedefs.h"
#include "SamAuxTag.h"
#include "CigarString.h"
#include "DNASequence.h"
#include "FastqQuality.h"


class SamRecord {
  // General fields (n >= 11)
  std::string      qname; // Field 1
  uint32_t          flag; // Field 2
  std::string      rname; // Field 3
  uint32_t           pos; // Field 4
  uint32_t          mapq; // Field 5
  CigarStringSPtr  cigar; // Field 6
  std::string      rnext; // Field 7
  uint32_t         pnext; // Field 8
  uint32_t          tlen; // Field 9
  DNASequenceSPtr    seq; // Field 10
  FastqQualitySPtr  qual; // Field 11
  std::vector<SamAuxTagSPtr> tags; // Field 12 -
  
  // Mutations
  std::vector<std::pair<std::string, size_t>> _seq_ops;
  std::vector<std::shared_ptr<std::pair<size_t, unsigned int>>> _mismatches; // TODO move Typedefs.h
  std::vector<std::shared_ptr<std::pair<size_t, std::string>>>  _insertions; // TODO move Typedefs.h
  std::vector<std::shared_ptr<std::pair<size_t, std::string>>>   _deletions; // TODO move Typedefs.h

  // Count mutations
  void collectSeqOps();
  void countMutations();

public:
  SamRecord() = delete;
  SamRecord(SamRecord&& o);
  SamRecord(const SamRecord& o);
  SamRecord(
    const std::string&        qname,
    const uint32_t             flag,
    const std::string&        rname,
    const uint32_t              pos,
    const uint32_t             mapq,
    const CigarStringSPtr     cigar,
    const std::string&        rnext,
    const uint32_t            pnext,
    const uint32_t             tlen,
    const DNASequenceSPtr       seq,
    const FastqQualitySPtr     qual,
    std::vector<SamAuxTagSPtr> tags
  );
  ~SamRecord() = default;

  // Streams
  friend std::ostream& operator<<(std::ostream&, const SamRecord&);
  friend std::ostream& operator<<(std::ostream&, const std::shared_ptr<SamRecord>);

  // Comparisons
  bool operator==(const SamRecord&) const;
  bool operator!=(const SamRecord&) const;
  bool operator> (const SamRecord&) const;
  bool operator>=(const SamRecord&) const;
  bool operator< (const SamRecord&) const;
  bool operator<=(const SamRecord&) const;

  // Getters
  const std::string  get_qname() const;
  const uint32_t      get_flag() const;
  const std::string  get_rname() const;
  const uint32_t       get_pos() const;
  const uint32_t      get_mapq() const;
  const CigarString  get_cigar() const;
  const std::string  get_rnext() const;
  const uint32_t     get_pnext() const;
  const uint32_t      get_tlen() const;
  const DNASequence    get_seq() const;
  const FastqQuality  get_qual() const;

  const std::vector<SamAuxTagSPtr> get_tags(std::string) const;
  const SamAuxTagSPtr get_tag(std::string) const;

  const std::string& get_cell_barcode() const;
  const std::string& get_unique_molecule_id() const;
  const std::string& get_md() const;
  const int get_nm() const;

  // Mapping traits
  bool is_paired() const;
  bool is_proper_pair() const;
  bool is_unmapped() const;
  bool is_mate_unmapped() const;
  bool is_reverse() const;
  bool is_mate_reverse() const;
  bool is_read1() const;
  bool is_read2() const;
  bool is_secondary() const;
  bool is_qc_failed() const;
  bool is_duplicates() const;
  bool is_supplementary() const;

  // Sequence operations
  StrSPtr to_string();
  StrSPtr reverse_complement();
  std::string get_reference_sequence();

  std::vector<std::shared_ptr<std::pair<size_t, unsigned int>>> get_mismatches() const;
  std::vector<std::shared_ptr<std::pair<size_t, unsigned int>>> get_mismatches(int) const;
  std::vector<std::shared_ptr<std::pair<size_t, unsigned int>>> get_mismatches(int, int) const;

  // Nascent RNA
  bool skipRead(uint32_t, uint32_t, uint8_t, bool) const;
  bool hasT2C(uint8_t max_mm, uint8_t trim_head, uint8_t trim_tail) const;
  typedef std::shared_ptr<SamRecord> SamRecordSPtr;
};

typedef std::pair<size_t, std::string> Indel;
typedef std::pair<size_t, unsigned int> Mismatch;
typedef std::vector<std::shared_ptr<std::pair<size_t, unsigned int>>> MismatchVec;
typedef std::vector<std::shared_ptr<std::pair<size_t, std::string>>>  IndelVec;


/*
 * Inline functions to convert mutation codes
 *
***/
inline unsigned int encode_mutations(char r, char q) {
  unsigned int d = 0;
  if (q == 'T') { d |=  1; } else if (q == 'G') { d |=  2; } else if (q == 'C') { d |=  4; } else if (q == 'A') { d |=   8; }
  if (r == 'T') { d |= 16; } else if (r == 'G') { d |= 32; } else if (r == 'C') { d |= 64; } else if (r == 'A') { d |= 128; }
  return d;
}

inline unsigned int encode_mutations(std::array<char, 2> m) { return encode_mutations(m.at(0), m.at(1)); }
inline unsigned int encode_mutations(std::shared_ptr<std::array<char, 2>> m) { return encode_mutations(*m); }

inline std::shared_ptr<std::array<char, 2>> decode_mutations(unsigned int d) {
  char q{'N'}, r{'N'};
  if (d &  1) { q = 'T'; } else if (d &  2) { q = 'G'; } else if (d &  4) { q = 'C'; } else if (d & 8) { q = 'A'; }
  d >>= 4;
  if (d &  1) { r = 'T'; } else if (d &  2) { r = 'G'; } else if (d &  4) { r = 'C'; } else if (d & 8) { r = 'A'; }
  return std::shared_ptr<std::array<char, 2>>(new std::array<char, 2>({r, q}));
}
