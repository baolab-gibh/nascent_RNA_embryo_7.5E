#include "SamRecord.h"
#include <bitset>

/*
* SamRecord
* */
SamRecord::SamRecord(
  const SamRecord& o
): qname(o.qname),   flag(o.flag), rname(o.rname), pos(o.pos), mapq(o.mapq), cigar(o.cigar),
   rnext(o.rnext), pnext(o.pnext),   tlen(o.tlen), seq(o.seq), qual(o.qual),   tags(o.tags)
{
  this->collectSeqOps();
  this->countMutations();
}

SamRecord::SamRecord(
  SamRecord&& o
) : qname(std::move(o.qname)),              flag(o.flag), rname(std::move(o.rname)),              pos(o.pos),
                 mapq(o.mapq), cigar(std::move(o.cigar)), rnext(std::move(o.rnext)),          pnext(o.pnext),
                 tlen(o.tlen),     seq(std::move(o.seq)),   qual(std::move(o.qual)), tags(std::move(o.tags))
{
  this->collectSeqOps();
  this->countMutations();
}

SamRecord::SamRecord(
  const std::string& qname, const uint32_t         flag, const std::string&    rname,     const uint32_t          pos,
  const uint32_t      mapq, const CigarStringSPtr cigar, const std::string&    rnext,     const uint32_t        pnext,
  const uint32_t      tlen, const DNASequenceSPtr   seq, const FastqQualitySPtr qual, std::vector<SamAuxTagSPtr> tags
): qname(qname),   flag(flag), rname(rname),   pos(pos), mapq(mapq),          cigar(cigar),
   rnext(rnext), pnext(pnext),   tlen(tlen),   seq(seq), qual(qual), tags(std::move(tags))
{
  this->collectSeqOps();
  this->countMutations();
}

// Streams
std::ostream& operator<<(std::ostream& o, const SamRecord& br)
{
  o << br.qname << "\t" << br.flag  << "\t" << br.rname << "\t" << br.pos << "\t"
    << br.mapq  << "\t" << br.cigar << "\t" << br.rnext << "\t" << br.pnext << "\t"
    << br.tlen  << "\t" << *br.seq  << "\t" << *br.qual;
  for (const auto& tag : br.tags) { o << "\t" << tag; }
  return o;
}

std::ostream& operator<<(std::ostream& o, const std::shared_ptr<SamRecord> br) { return o << *br; }

// Comparisons
bool SamRecord::operator==(const SamRecord& o) const { return this->rname == o.rname && this->pos == o.pos && this->tlen == o.tlen; }
bool SamRecord::operator!=(const SamRecord& o) const { return !(*this == o); }
bool SamRecord::operator> (const SamRecord& o) const { return this->rname == o.rname ? this->pos > o.pos : this->rname > o.rname; }
bool SamRecord::operator>=(const SamRecord& o) const { return *this > o || *this == o; }
bool SamRecord::operator< (const SamRecord& o) const { return this->rname == o.rname ? this->pos < o.pos : this->rname < o.rname; }
bool SamRecord::operator<=(const SamRecord& o) const { return *this < o || *this == o; }

// Getters
const std::string  SamRecord::get_qname() const { return this->qname;  }
const uint32_t     SamRecord::get_flag()  const { return this->flag;   }
const std::string  SamRecord::get_rname() const { return this->rname;  }
const uint32_t     SamRecord::get_pos()   const { return this->pos;    }
const uint32_t     SamRecord::get_mapq()  const { return this->mapq;   }
const CigarString  SamRecord::get_cigar() const { return *this->cigar; }
const std::string  SamRecord::get_rnext() const { return this->rnext;  }
const uint32_t     SamRecord::get_pnext() const { return this->pnext;  }
const uint32_t     SamRecord::get_tlen()  const { return this->tlen;   }
const DNASequence  SamRecord::get_seq()   const { return *this->seq;   }
const FastqQuality SamRecord::get_qual()  const { return *this->qual;  }

const std::vector<SamAuxTagSPtr>  SamRecord::get_tags(std::string tag) const { return this->tags; }

// Get given tag by tag name.
const SamAuxTagSPtr SamRecord::get_tag(std::string tag) const {
  for (const auto& t : this->tags) { if (t->get_name() == tag) { return t; } }
  throw std::runtime_error("Tag " + tag + " not found!");
}

// Get the value of NM tag.
const int SamRecord::get_nm() const { return this->get_tag("NM")->get_val_int(); }

// Get the value of MD tag.
const std::string& SamRecord::get_md() const { return this->get_tag("MD")->get_value(); }

// Get the value of CR tag.
const std::string& SamRecord::get_cell_barcode() const { return this->get_tag("CR")->get_value(); }

// Get the value of UR tag.
const std::string& SamRecord::get_unique_molecule_id() const { return this->get_tag("UR")->get_value(); }

// collection sequence operations, including CIGAR and MD tag.
void SamRecord::collectSeqOps() {
  auto cigar_ops = this->get_cigar();
  auto md_ops = this->get_md();

  for (const auto& op : cigar_ops.get_op()) {
    if (op.second == CigarOperation::INSERTION) {
      this->_seq_ops.push_back(std::make_pair("I", op.first));
    } else if (op.second == CigarOperation::DELETION) {
      this->_seq_ops.push_back(std::make_pair("D", op.first));
    } else if (op.second == CigarOperation::MATCH) {
      this->_seq_ops.push_back(std::make_pair("M", op.first));
    } else if (op.second == CigarOperation::REFERENCE_SKIP) {
      this->_seq_ops.push_back(std::make_pair("N", op.first));
    }
  }
}

// Count mutations
void SamRecord::countMutations() {
  size_t pos = 0;
  std::string* match_len = new std::string("");

  // Count mutations by MD tag, including mismatches and deletions.
  std::cout << "------------------------------\n";
  std::cout << this->qname << "\t" << this->get_md() << "\t" << this->get_cigar() << "\t" << this->seq << std::endl;

  // Count mutations by CIGAR, including insertions.
  for (const auto& op : this->get_cigar().get_op()) {
    if (op.second == CigarOperation::INSERTION) {
      StrSPtr insertion_str(new std::string(""));
      for (size_t i = 0; i < op.first; i++) {
        insertion_str->push_back(this->seq->at(pos + i));
      }
      this->_insertions.push_back(std::make_shared<Indel>(pos, *insertion_str));
    }
    pos += op.first;
  }

  pos = 0;
  for (auto it = this->get_md().begin(); it != this->get_md().end(); ++it) {
    if ('0' <= *it && *it <= '9') {
      match_len->push_back(*it);
    } else {
      if (match_len->size() != 0) {
        pos += std::stoi(*match_len);
        match_len->clear();
      }

      pos++;
      if (*it == '^') {
        std::shared_ptr<std::string> deletion_str(new std::string(""));
        pos--;
        it++;
        while (it != this->get_md().end() && !('0' <= *it && *it <= '9')) {
          deletion_str->push_back(*it);
          it++;
        }

        this->_deletions.push_back(std::make_shared<Indel>(pos, *deletion_str));
        it--;
        continue;
      }

      std::cout << "running position: " << pos << " " << *it << "->" << this->seq->at(pos-1) << std::endl;

      auto mmp = std::make_shared<Mismatch>(pos, encode_mutations(*it, this->seq->at(pos - 1))); // reference -> query (read)
      this->_mismatches.push_back(mmp);
    }
  }
  delete match_len;
}

MismatchVec SamRecord::get_mismatches(int start, int end) const {
  if (start < 0) { start = 0; }
  if (end > this->seq->size()) { end = this->seq->size(); }

  auto sub_mmp = MismatchVec();
  if (this->_mismatches.size() != 0) {
    for (const auto& mmp : this->_mismatches) {
      if (start <= mmp->first && mmp->first < end) {
        sub_mmp.push_back(mmp);
      }
    }
  }
  return sub_mmp;
}
MismatchVec SamRecord::get_mismatches(int start) const { return this->get_mismatches(start, this->seq->size()); }
MismatchVec SamRecord::get_mismatches() const { return this->get_mismatches(0, this->seq->size()); }

// Mapping traits
bool SamRecord::is_paired()        const { return this->flag &    1; }
bool SamRecord::is_proper_pair()   const { return this->flag &    2; }
bool SamRecord::is_unmapped()      const { return this->flag &    4; }
bool SamRecord::is_mate_unmapped() const { return this->flag &    8; }
bool SamRecord::is_reverse()       const { return this->flag &   16; }
bool SamRecord::is_mate_reverse()  const { return this->flag &   32; }
bool SamRecord::is_read1()         const { return this->flag &   64; }
bool SamRecord::is_read2()         const { return this->flag &  128; }
bool SamRecord::is_secondary()     const { return this->flag &  256; }
bool SamRecord::is_qc_failed()     const { return this->flag &  512; }
bool SamRecord::is_duplicates()    const { return this->flag & 1024; }
bool SamRecord::is_supplementary() const { return this->flag & 2048; }

// Sequence operations
StrSPtr SamRecord::reverse_complement() { return this->seq->reverse_complement(); }
std::string get_reference_sequence() { return ""; }

StrSPtr SamRecord::to_string() {
  StrSPtr bam_str = std::make_shared<std::string>("");

  bam_str->append(this->qname);                  bam_str->append("\t");
  bam_str->append(std::to_string(this->flag));   bam_str->append("\t");
  bam_str->append(this->rname);                  bam_str->append("\t");
  bam_str->append(std::to_string(this->pos));    bam_str->append("\t");
  bam_str->append(std::to_string(this->mapq));   bam_str->append("\t");
  bam_str->append(*(this->cigar->to_string()));  bam_str->append("\t");
  bam_str->append(this->rnext);                  bam_str->append("\t");
  bam_str->append(std::to_string(this->pnext));  bam_str->append("\t");
  bam_str->append(std::to_string(this->tlen));   bam_str->append("\t");
  bam_str->append(this->get_seq().to_string());  bam_str->append("\t");
  bam_str->append(this->get_qual().to_string()); bam_str->append("\t");
  for (const auto& tag : this->tags) {
    bam_str->append(tag->to_string());
    bam_str->append("\t");
  }
  if (bam_str->back() == '\t') { bam_str->pop_back(); }

  return bam_str;
}

bool SamRecord::skipRead(uint32_t excl_flag, uint32_t incl_flag, uint8_t max_mm, bool excl_unpaired) const {
  uint8_t read_skip_code = 0;
  int nm_counts = this->get_tag("nM")->get_val_int();

  if ((this->flag & excl_flag) != 0) read_skip_code |= 1; // Skip due to containing non-required flag
  if ((this->flag > 0 && (this->flag & incl_flag) == 0)) read_skip_code |= 2; // Skip due to missing required flags
  if (excl_unpaired && this->is_paired() == 0) { read_skip_code |= 4; } // Skip due to unpaired reads
  if (nm_counts < 1) read_skip_code |= 8; // Remove reads without mismatches.
  if (nm_counts > max_mm) read_skip_code |= 16; // Remove reads without mismatches.
  if (!this->seq->is_valid()) read_skip_code |= 32; // Remove reads without base sequence if this->seq->is_valid();
  if (!this->qual->is_valid()) read_skip_code |= 64; // Remove reads without base quality
  
  return read_skip_code != 0 ? true : false;
}

bool SamRecord::hasT2C(uint8_t max_mm, uint8_t trim_head, uint8_t trim_tail) const {
  uint8_t n_t2c = 0;
  trim_tail = this->seq->size() - trim_tail;
  for (auto pm : this->get_mismatches(trim_head, trim_tail)) {
    auto x = decode_mutations(pm->second);
    std::cout << "Selected mutations: " << x->at(0) << "->" << x->at(1) << std::endl;
    if ((pm->second & (16 | 1)) == 17) { n_t2c++; }
  }

  return n_t2c > 0;
}
