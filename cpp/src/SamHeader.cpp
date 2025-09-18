#include "SamHeader.h"
SamHeader::SamHeader(const std::string& fn): filename(fn) {
    samFile *fp = sam_open(fn.c_str(), "r");
    this->header = sam_hdr_read(fp);
  }

SamHeader::SamHeader(std::string fn, sam_hdr_t *header): filename(fn), header(header) { this->fp = nullptr; }
SamHeader::~SamHeader() {
  if (this->header != nullptr) {
    if (this->fp == nullptr) {
      sam_hdr_destroy(this->header);
    } else {
      sam_close(this->fp);
    }
  }
}
