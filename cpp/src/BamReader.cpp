#include "BamReader.h"

BamReader::BamReader(const std::string& filename): filename(filename)
{
  const htsFormat *fmt;

  this->fp = sam_open(filename.c_str(), "r");
  fmt = hts_get_format(this->fp);
  if (fmt->format != htsExactFormat::bam) { std::cerr << "Not a BAM file" << std::endl; }
  this->header = sam_hdr_read(this->fp);
  this->bam_rec = bam_init1();
}

BamReader::~BamReader() { hts_close(this->fp); }

// Read one record
std::shared_ptr<SamRecord> BamReader::next()
{
  int        ret{0};
  uint32_t  flag{0};
  uint32_t   pos{0};
  uint32_t pnext{0};
  uint32_t  mapq{0};
  uint32_t  tlen{0};
  std::vector<SamAuxTagSPtr>  tags;

  std::string    qname{"*"};
  std::string    rname{"*"};
  CigarStringSPtr     cigar;
  std::string    rnext{"*"};
  DNASequenceSPtr       seq;
  FastqQualitySPtr     qual;

  const bam1_core_t *bam_core{&this->bam_rec->core};

  ret = sam_read1(this->fp, this->header, this->bam_rec);
  if (ret < 0) {
    throw std::runtime_error("Failed to read BAM files! Exit 1.");
    exit(1);
  }

  qname = std::string((char*)this->bam_rec->data, bam_core->l_qname); // qname
  if (qname.back() == '\0') { qname.pop_back(); }
  flag = bam_core->flag; // flag
  rname = std::string(this->header->target_name[bam_core->tid]); // chrom
  if (rname.back() == '\0') { rname.pop_back(); }
  pos = bam_core->pos + 1; // pos
  mapq = bam_core->qual; // qual
  cigar = std::make_shared<CigarString>(this->bam_rec); // cigar

  if (bam_core->mtid == bam_core->tid) { // mate chr
    rnext = "=";
  } else if (bam_core->mtid >= 0) {
    rnext = std::string(this->header->target_name[bam_core->mtid]);
  }

  pnext = bam_core->mpos + 1; // mate pos
  tlen = bam_core->isize; // template len
  seq = std::make_shared<DNASequence>(this->bam_rec); // sequence
  qual = std::make_shared<FastqQuality>(this->bam_rec); // qual

  // aux
  uint8_t* end = this->bam_rec->data + this->bam_rec->l_data;
  for (auto aux = bam_aux_first(this->bam_rec); aux; aux = bam_aux_next(this->bam_rec, aux)) {
    tags.push_back(this->get_tag(aux, end));
  }
  this->count++;

  SamRecord::SamRecordSPtr br(new SamRecord(qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, tags));
  return br;
}

// Obtain auxiliary tags.
SamAuxTagSPtr BamReader::get_tag(const uint8_t* aux, const uint8_t *end) {
  auto tag_name = std::string{""}, tag_type = std::string{""}, tag_value = std::string{""};

  tag_name = std::string{(char*)(aux - 2), 2};

  uint8_t _type = *aux;
  if (_type == 'C') {
    tag_type = "i";
    tag_value = std::to_string(*(aux + 1));
  } else if (_type == 'c') {
    tag_type = "i";
    tag_value = std::to_string(le_to_i8(aux));
  } else if (_type == 'S') {
    tag_type = "i";
    tag_value = std::to_string(le_to_u16(aux));
  } else if (_type == 's') {
    tag_type = "i";
    tag_value = std::to_string(le_to_i16(aux));
  } else if (_type == 'I') {
    tag_type = "i";
    tag_value = std::to_string(le_to_u32(aux));
  } else if (_type == 'i') {
    tag_type = "i";
    tag_value = std::to_string(le_to_i32(aux));
  } else if (_type == 'A') {
    tag_type = "A";
    tag_value = (char*)(aux + 1);
  } else if (_type == 'f') {
    tag_type = "f";
    tag_value = std::to_string((double)le_to_float(aux));
  } else if (_type == 'd') {
    tag_type = "d";
    tag_value = std::to_string(le_to_double(aux));
  } else if (_type == 'Z' || _type == 'H') {
    tag_type = "Z";
    tag_value = (char*)(aux + 1);
  } else if (_type == 'B') {
    uint32_t n;
    uint8_t sub_type = *(aux++);
    unsigned sub_type_size;

    // or externalise sam.c's aux_type2size function?
    switch (sub_type) {
      case 'A': case 'c': case 'C':
        sub_type_size = 1;
        break;
      case 's': case 'S':
        sub_type_size = 2;
        break;
      case 'i': case 'I': case 'f':
        sub_type_size = 4;
        break;
      default:
        sub_type_size = 0;
        break;
    }

    if (sub_type_size == 0 || end - aux < 4)
      throw std::runtime_error("Truncated aux data");

    n = le_to_u32(aux);
    aux += 4; // now points to the start of the array
    if ((size_t)(end - aux) / sub_type_size < n)
      throw std::runtime_error("Truncated aux data");

    tag_type = "B:" + std::string(sub_type, 1);

    switch (sub_type) {
      case 'c':
        for (uint32_t i = 0; i < n; ++i) {
          tag_value += ",";
          tag_value += *(int8_t*)aux;
          ++aux;
        }
        break;
      case 'C':
        for (uint32_t i = 0; i < n; ++i) {
          tag_value += ",";
          tag_value += *(uint8_t*)aux;
          ++aux;
        }
        break;
      case 's':
        for (uint32_t i = 0; i < n; ++i) {
          tag_value += ",";
          tag_value += std::to_string(le_to_i16(aux));
          aux += 2;
        }
        break;
      case 'S':
        for (uint32_t i = 0; i < n; ++i) {
          tag_value += ",";
          tag_value += std::to_string(le_to_u16(aux));
          aux += 2;
        }
        break;
      case 'i':
        for (uint32_t i = 0; i < n; ++i) {
          tag_value += ",";
          tag_value += std::to_string(le_to_i32(aux));
          aux += 4;
        }
        break;
      case 'I':
        for (uint32_t i = 0; i < n; ++i) {
          tag_value += ",";
          tag_value += std::to_string(le_to_u32(aux));
          aux += 4;
        }
        break;
      case 'f':
        for (uint32_t i = 0; i < n; ++i) {
          tag_value += ",";
          tag_value += std::to_string((double)le_to_float(aux));
          aux += 4;
        }
        break;
      default:
        throw std::runtime_error("Unknown aux type");
    }
  } else { // Unknown type
    throw std::runtime_error("Unknown aux type");
  }

  SamAuxTagSPtr sam_tag(new SamAuxTag<std::string>(tag_name, tag_type, tag_value));

  return sam_tag;
}

// Return current line in SAM format.
size_t BamReader::current_line() { return this->count; }
